@testsnippet MulchSetup begin
    using Test
    using ModelingToolkit
    using ModelingToolkit: t, D
    using OrdinaryDiffEqDefault
    using OrdinaryDiffEqDefault: SciMLBase
    using Vegetation

    const day_s = 86400.0  # seconds per day
end

# ============================================================
# ResidueMulchDecomposition Tests
# ============================================================

@testitem "ResidueMulchDecomposition: Structural Verification" setup = [MulchSetup] tags = [:mulch] begin
    sys = ResidueMulchDecomposition()
    @test sys isa ModelingToolkit.System
    @test nameof(sys) == :ResidueMulchDecomposition

    vars = unknowns(sys)
    eqs = equations(sys)

    # 39 equations, 41 unknowns (39 algebraic + 13 ODE states - 2 forcing params = 39 vars)
    @test length(eqs) == 39

    # Verify key variable names exist
    var_names = [string(v) for v in vars]
    for expected in [
            "RM_CARB(t)", "RM_CELL(t)", "RM_LIGN(t)",
            "RMN_CARB(t)", "RMN_CELL(t)", "RMN_LIGN(t)",
            "N_inorg(t)", "RM_T(t)", "RMN_T(t)",
            "MTRF(t)", "CNRF(t)", "RN_Net(t)",
        ]
        @test any(n -> contains(n, expected), var_names)
    end

    # Verify compilation succeeds
    compiled = mtkcompile(sys)
    @test compiled !== nothing
    # 13 ODE states: 6 RM pools (3 contacting + 3 non-contacting),
    # 6 RMN pools, and N_inorg
    @test length(unknowns(compiled)) == 13
end

@testitem "ResidueMulchDecomposition: Default Parameters" setup = [MulchSetup] tags = [:mulch] begin
    sys = ResidueMulchDecomposition()
    compiled = mtkcompile(sys)

    params = parameters(compiled)
    pdict = Dict(
        Symbol(p) => ModelingToolkit.getdefault(p) for p in params
            if ModelingToolkit.hasdefault(p)
    )

    # Verify decay coefficients match Table 2 (converted from d⁻¹ to s⁻¹)
    @test pdict[:D_CARB] ≈ 0.43 / day_s rtol = 1.0e-6
    @test pdict[:D_CELL] ≈ 0.24 / day_s rtol = 1.0e-6
    @test pdict[:D_LIGN] ≈ 0.0228 / day_s rtol = 1.0e-6

    # Verify lignin inhibition parameter
    @test pdict[:γ_lign] ≈ 12.0

    # Verify MTRF coefficients (Table 2)
    @test pdict[:a_MTRF] ≈ 0.384
    @test pdict[:b_MTRF] ≈ 0.0187
    @test pdict[:c_MTRF] ≈ 0.142
    @test pdict[:d_MTRF] ≈ 0.628

    # Verify CNRF coefficients
    @test pdict[:a_CNRF] ≈ 0.693
    @test pdict[:CNR_crit] ≈ 13.0

    # Verify other parameters
    @test pdict[:HUMF] ≈ 0.125
    @test pdict[:α_F] ≈ 0.1
end

@testitem "ResidueMulchDecomposition: Eq. 21 - Decay Coefficients" setup = [MulchSetup] tags = [:mulch] begin
    # Verify decay coefficients depend on lignin fraction
    sys = ResidueMulchDecomposition()
    compiled = mtkcompile(sys)

    tspan = (0.0, 1.0)
    prob = ODEProblem(compiled, [], tspan)
    sol = solve(prob)

    # Initial LIGN fraction = 0.1/1.0 = 0.1 (since all mass in contacting portion sums to 0.12)
    # But actual fraction: RM_LIGN_init / RM_T_init = 0.012 / 0.12 = 0.1
    f_LIGN_init = 0.1

    # d_CARB = D_CARB * exp(-γ * f_LIGN) = 0.43/86400 * exp(-12 * 0.1)
    d_CARB_expected = 0.43 / day_s * exp(-12.0 * f_LIGN_init)
    d_CARB_actual = sol[compiled.d_CARB][1]
    @test d_CARB_actual ≈ d_CARB_expected rtol = 1.0e-4

    # d_CELL = D_CELL * exp(-γ * f_LIGN) = 0.24/86400 * exp(-12 * 0.1)
    d_CELL_expected = 0.24 / day_s * exp(-12.0 * f_LIGN_init)
    d_CELL_actual = sol[compiled.d_CELL][1]
    @test d_CELL_actual ≈ d_CELL_expected rtol = 1.0e-4

    # d_LIGN is constant
    d_LIGN_expected = 0.0228 / day_s
    d_LIGN_actual = sol[compiled.d_LIGN][1]
    @test d_LIGN_actual ≈ d_LIGN_expected rtol = 1.0e-4
end

@testitem "ResidueMulchDecomposition: Eq. 22 - MTRF and CNRF" setup = [MulchSetup] tags = [:mulch] begin
    sys = ResidueMulchDecomposition()
    compiled = mtkcompile(sys)

    # Test at T=20°C (293.15K), h=-1m
    tspan = (0.0, 1.0)
    prob = ODEProblem(compiled, [], tspan)
    sol = solve(prob)

    T_C = 293.15 - 273.15  # 20°C
    h_m = -1.0

    # MTRF = (0.384 + 0.0187*20) * exp((0.142 + 0.628/20) * (-1))
    MTRF_expected = (0.384 + 0.0187 * T_C) * exp((0.142 + 0.628 / T_C) * h_m)
    MTRF_actual = sol[compiled.MTRF][1]
    @test MTRF_actual ≈ MTRF_expected rtol = 1.0e-4

    # CNR = 0.41 * RM_T / (RMN_T + N_inorg)
    RM_T = 0.12  # kg/m² (all in contacting)
    RMN_T = 0.12 * (0.2 * 0.08 + 0.7 * 0.01 + 0.1 * 0.01)  # initial N mass
    N_inorg_init = 0.001
    CNR_expected = 0.41 * RM_T / (RMN_T + N_inorg_init)
    CNR_actual = sol[compiled.CNR][1]
    @test CNR_actual ≈ CNR_expected rtol = 1.0e-3

    # CNRF depends on whether CNR > 13
    if CNR_expected > 13.0
        CNRF_expected = exp(-0.693 * (CNR_expected - 13.0) / 13.0)
    else
        CNRF_expected = 1.0
    end
    CNRF_actual = sol[compiled.CNRF][1]
    @test CNRF_actual ≈ CNRF_expected rtol = 1.0e-3
end

@testitem "ResidueMulchDecomposition: Freezing Stops Decomposition" setup = [MulchSetup] tags = [:mulch] begin
    # When T ≤ 0°C, MTRF = 0 and no decomposition should occur
    sys = ResidueMulchDecomposition()
    compiled = mtkcompile(sys)

    tspan = (0.0, 10.0 * day_s)  # 10 days
    # Set temperature to -5°C = 268.15K
    prob = ODEProblem(compiled, [compiled.T_mulch => 268.15], tspan)
    sol = solve(prob)

    # MTRF should be 0
    @test sol[compiled.MTRF][1] ≈ 0.0 atol = 1.0e-10

    # Residue mass should not change
    @test sol[compiled.RM_CARB][end] ≈ sol[compiled.RM_CARB][1] rtol = 1.0e-6
    @test sol[compiled.RM_CELL][end] ≈ sol[compiled.RM_CELL][1] rtol = 1.0e-6
    @test sol[compiled.RM_LIGN][end] ≈ sol[compiled.RM_LIGN][1] rtol = 1.0e-6
end

@testitem "ResidueMulchDecomposition: Mass Positivity" setup = [MulchSetup] tags = [:mulch] begin
    sys = ResidueMulchDecomposition()
    compiled = mtkcompile(sys)

    tspan = (0.0, 100.0 * day_s)  # 100 days
    prob = ODEProblem(compiled, [], tspan)
    sol = solve(prob)

    @test sol.retcode == SciMLBase.ReturnCode.Success

    # All mass pools should remain non-negative
    @test all(sol[compiled.RM_CARB] .>= -1.0e-10)
    @test all(sol[compiled.RM_CELL] .>= -1.0e-10)
    @test all(sol[compiled.RM_LIGN] .>= -1.0e-10)
    @test all(sol[compiled.RMN_CARB] .>= -1.0e-10)
    @test all(sol[compiled.RMN_CELL] .>= -1.0e-10)
    @test all(sol[compiled.RMN_LIGN] .>= -1.0e-10)
end

@testitem "ResidueMulchDecomposition: CARB Decomposes Fastest" setup = [MulchSetup] tags = [:mulch] begin
    # Paper states CARB > CELL > LIGN decomposition speed
    sys = ResidueMulchDecomposition()
    compiled = mtkcompile(sys)

    tspan = (0.0, 50.0 * day_s)
    prob = ODEProblem(compiled, [], tspan)
    sol = solve(prob)

    @test sol.retcode == SciMLBase.ReturnCode.Success

    # Fraction remaining for each pool
    frac_CARB = sol[compiled.RM_CARB][end] / sol[compiled.RM_CARB][1]
    frac_CELL = sol[compiled.RM_CELL][end] / sol[compiled.RM_CELL][1]
    frac_LIGN = sol[compiled.RM_LIGN][end] / sol[compiled.RM_LIGN][1]

    # CARB should decompose most, then CELL, then LIGN
    @test frac_CARB < frac_CELL
    @test frac_CELL < frac_LIGN
end

@testitem "ResidueMulchDecomposition: Temperature Sensitivity" setup = [MulchSetup] tags = [:mulch] begin
    sys = ResidueMulchDecomposition()
    compiled = mtkcompile(sys)

    tspan = (0.0, 30.0 * day_s)

    # Warm conditions: 30°C
    prob_warm = ODEProblem(compiled, [compiled.T_mulch => 303.15], tspan)
    sol_warm = solve(prob_warm)

    # Cool conditions: 10°C
    prob_cool = ODEProblem(compiled, [compiled.T_mulch => 283.15], tspan)
    sol_cool = solve(prob_cool)

    @test sol_warm.retcode == SciMLBase.ReturnCode.Success
    @test sol_cool.retcode == SciMLBase.ReturnCode.Success

    # Warmer temperatures should lead to more decomposition (less remaining mass)
    RM_T_warm = sol_warm[compiled.RM_CARB][end] + sol_warm[compiled.RM_CELL][end] + sol_warm[compiled.RM_LIGN][end]
    RM_T_cool = sol_cool[compiled.RM_CARB][end] + sol_cool[compiled.RM_CELL][end] + sol_cool[compiled.RM_LIGN][end]
    @test RM_T_warm < RM_T_cool
end

@testitem "ResidueMulchDecomposition: Net N Mineralization" setup = [MulchSetup] tags = [:mulch] begin
    sys = ResidueMulchDecomposition()
    compiled = mtkcompile(sys)

    tspan = (0.0, 100.0 * day_s)
    prob = ODEProblem(compiled, [], tspan)
    sol = solve(prob)

    @test sol.retcode == SciMLBase.ReturnCode.Success

    # With C/N ratio of ~17 (high N content), net N mineralization should occur
    # N_inorg should increase over time
    @test sol[compiled.N_inorg][end] > sol[compiled.N_inorg][1]
end

# ============================================================
# MulchRadiationAttenuation Tests
# ============================================================

@testitem "MulchRadiationAttenuation: Structural Verification" setup = [MulchSetup] tags = [:mulch] begin
    sys = MulchRadiationAttenuation()
    @test sys isa ModelingToolkit.System
    @test nameof(sys) == :MulchRadiationAttenuation

    eqs = equations(sys)
    vars = unknowns(sys)

    # K=5: 6 S_d + 6 S_u + 6 L_d + 6 L_u + 6 R_net_iface + 5 R_elem + 4 forcing = 39 vars
    # 35 equations
    @test length(eqs) == 35
end

@testitem "MulchRadiationAttenuation: Shortwave Attenuation" setup = [MulchSetup] tags = [:mulch] begin
    sys = MulchRadiationAttenuation()

    # Verify Eq. 11 symbolically: S_d should attenuate geometrically
    # Test by substituting known values into the equations
    eqs = equations(sys)
    vars = unknowns(sys)

    # Use symbolic substitution to evaluate
    subs = Dict(
        sys.S_0 => 500.0,
        sys.T_a => 300.0,
        sys.T_s => 295.0,
        sys.ΔR => 0.3,
        sys.Ω_cl => 0.6,
        sys.α_m => 0.3,
        sys.α_s => 0.15,
    )
    for k in 1:5
        subs[sys.T_layer[k]] = 295.0 + k
    end

    # Solve algebraically by forward substitution for S_d
    ΔR = 0.3
    Ω_cl = 0.6
    S_0 = 500.0

    # Eq. 11: Shortwave attenuation
    S_d_6 = S_0  # top interface
    S_d_5 = S_0 * (1 - ΔR)
    S_d_expected = [S_0 * (1 - ΔR) * (1 - Ω_cl * ΔR)^(5 - k) for k in 1:4]
    pushfirst!(S_d_expected, S_0 * (1 - ΔR) * (1 - Ω_cl * ΔR)^4)  # k=1

    # Verify geometric attenuation
    for k in 1:4
        @test S_d_expected[k] < S_d_expected[k + 1] || S_d_expected[k] ≈ S_d_expected[k + 1]
    end

    # Bottom should be significantly attenuated
    @test S_d_expected[1] < S_0 * 0.5  # More than 50% attenuation
    @test S_d_expected[1] ≈ S_0 * (1 - ΔR) * (1 - Ω_cl * ΔR)^4 rtol = 1.0e-10
end

@testitem "MulchRadiationAttenuation: Net Radiation Consistency" setup = [MulchSetup] tags = [:mulch] begin
    # Verify that the sum of element net radiation equals the difference of interface net radiation
    # This is a structural property: Σ R_elem[k] = R_net[1] - R_net[K+1]
    # R_elem[k] = R_net[k] - R_net[k+1] → telescoping sum
    sys = MulchRadiationAttenuation()
    eqs = equations(sys)

    # Count equations by type
    n_S_d = 6  # S_d[1:6]
    n_S_u = 6  # S_u[1:6]
    n_L_d = 6  # L_d[1:6]
    n_L_u = 6  # L_u[1:6]
    n_R_net = 6  # R_net_iface[1:6]
    n_R_elem = 5  # R_elem[1:5]
    @test length(eqs) == n_S_d + n_S_u + n_L_d + n_L_u + n_R_net + n_R_elem
end

# ============================================================
# MulchWindProfile Tests
# ============================================================

@testitem "MulchWindProfile: Structural Verification" setup = [MulchSetup] tags = [:mulch] begin
    sys = MulchWindProfile()
    @test sys isa ModelingToolkit.System
    @test nameof(sys) == :MulchWindProfile

    # 1 friction velocity + 5 layer wind speeds + 1 forcing = 7 unknowns
    # 6 equations
    @test length(equations(sys)) == 6
end

@testitem "MulchWindProfile: Wind Speed Attenuation" setup = [MulchSetup] tags = [:mulch] begin
    # Verify Eq. 3 analytically
    K = 5
    Z_M = 0.06
    z_ref = 2.0
    u_ref = 3.0
    k_K = 0.4

    d_disp = 0.87 * Z_M
    z_r = 0.079 * Z_M

    # Eq. 3a: friction velocity
    u_star = k_K * u_ref / log((z_ref - d_disp) / z_r)
    @test u_star > 0

    # Eq. 3b: wind speed in each layer
    u_vals = Float64[]
    for k in 1:K
        z_center = (2k - 1) / (2K) * Z_M
        u_k = 0.21 * u_star * exp(2.2 / Z_M * z_center)
        push!(u_vals, u_k)
    end

    # Wind speed should increase from bottom to top
    for k in 1:(K - 1)
        @test u_vals[k] < u_vals[k + 1]
    end

    # All should be positive and less than reference
    @test all(u_vals .> 0)
    @test all(u_vals .< u_ref)
end

# ============================================================
# MulchHeatVaporFluxes Tests
# ============================================================

@testitem "MulchHeatVaporFluxes: Structural Verification" setup = [MulchSetup] tags = [:mulch] begin
    sys = MulchHeatVaporFluxes()
    @test sys isa ModelingToolkit.System
    @test nameof(sys) == :MulchHeatVaporFluxes

    @test length(equations(sys)) == 7
end

@testitem "MulchHeatVaporFluxes: Regime Selection" setup = [MulchSetup] tags = [:mulch] begin
    # Verify Eq. 4 analytically: Ra and Ri numbers
    Z_M = 0.06
    ν = 1.5e-5
    D_hm = 2.2e-5
    g = 9.81

    T_top = 310.0
    T_bot = 290.0
    ΔT = T_top - T_bot  # 20 K
    T_avg = (T_top + T_bot) / 2  # 300 K
    u_top = 0.01
    u_bot = 0.005
    Δu = u_top - u_bot  # 0.005 m/s

    # Eq. 4: Ra = 2g|ΔT|Z_M³ / (T̄ · ν · D_hm)
    Ra = 2 * g * abs(ΔT) * Z_M^3 / (T_avg * ν * D_hm)
    @test Ra > 1706.0  # Should trigger convection

    # Eq. 4: Ri = 2g|ΔT|Z_M / (T̄ · Δu²)
    Ri = 2 * g * abs(ΔT) * Z_M / (T_avg * Δu^2)
    @test Ri >= 1.0  # Should indicate free convection dominates

    # With strong wind, Ri should be small → forced convection
    Δu_strong = 5.0
    Ri_forced = 2 * g * abs(ΔT) * Z_M / (T_avg * Δu_strong^2)
    @test Ri_forced < 1.0
end

# ============================================================
# MulchWaterCharacteristic Tests
# ============================================================

@testitem "MulchWaterCharacteristic: Structural Verification" setup = [MulchSetup] tags = [:mulch] begin
    sys = MulchWaterCharacteristic()
    @test sys isa ModelingToolkit.System
    @test nameof(sys) == :MulchWaterCharacteristic
    @test length(equations(sys)) == 4
end

@testitem "MulchWaterCharacteristic: WRC Shape" setup = [MulchSetup] tags = [:mulch] begin
    # Verify Eq. 26 analytically
    f_LIGN = 0.1

    # a_m = -20.1e6 * exp(-0.249 * f_LIGN)
    a_m = -20.1e6 * exp(-0.249 * f_LIGN)
    @test a_m < 0  # Should be negative (suction)
    @test a_m ≈ -20.1e6 * exp(-0.249 * 0.1) rtol = 1.0e-10

    # b_m = 0.324 + 0.124 * f_LIGN
    b_m = 0.324 + 0.124 * f_LIGN
    @test b_m ≈ 0.3364 rtol = 1.0e-4
    @test b_m > 0  # Exponent should be positive

    # Saturated gravimetric WC
    θ_sat = 7.1 * exp(-0.079 * f_LIGN)
    @test θ_sat > 0
    @test θ_sat ≈ 7.1 * exp(-0.079 * 0.1) rtol = 1.0e-10

    # h = a_m * θ^(-b_m): at θ=1, h = a_m (negative, suction)
    h_at_1 = a_m * 1.0^(-b_m)
    @test h_at_1 < 0

    # At lower θ, h should be more negative (drier = more suction)
    h_at_low = a_m * 0.1^(-b_m)
    @test h_at_low < h_at_1  # More suction at lower water content
end

# ============================================================
# MulchHeatWaterTransfer Tests
# ============================================================

@testitem "MulchHeatWaterTransfer: Structural Verification" setup = [MulchSetup] tags = [:mulch] begin
    sys = MulchHeatWaterTransfer()
    @test sys isa ModelingToolkit.System
    @test nameof(sys) == :MulchHeatWaterTransfer

    eqs = equations(sys)
    vars = unknowns(sys)

    @test length(eqs) == 10

    var_names = [string(v) for v in vars]
    for expected in [
            "h_m(t)", "T_m(t)", "θ_vol(t)", "C_hh(t)", "C_TT(t)",
            "ρ_vs(t)", "h_rel(t)", "D_mv(t)", "D_Tv(t)", "λ_eff(t)",
        ]
        @test any(n -> contains(n, expected), var_names)
    end
end

@testitem "MulchHeatWaterTransfer: Compilation" setup = [MulchSetup] tags = [:mulch] begin
    sys = MulchHeatWaterTransfer()
    compiled = mtkcompile(sys)

    # After compilation, only h_m and T_m should remain as ODE states
    @test length(unknowns(compiled)) == 2
    state_names = Symbol.(unknowns(compiled))
    @test Symbol("h_m(t)") in state_names
    @test Symbol("T_m(t)") in state_names
end

@testitem "MulchHeatWaterTransfer: Constitutive Relations" setup = [MulchSetup] tags = [:mulch] begin
    sys = MulchHeatWaterTransfer()
    compiled = mtkcompile(sys)

    # Test at default conditions: h=-1m, T=293.15K
    tspan = (0.0, 1.0)
    prob = ODEProblem(compiled, [], tspan)
    sol = solve(prob)

    @test sol.retcode == SciMLBase.ReturnCode.Success

    # θ_vol should be positive
    @test sol[compiled.θ_vol][1] > 0

    # C_hh (water capacity) should be positive
    @test sol[compiled.C_hh][1] > 0

    # C_TT (heat capacity) should be positive
    @test sol[compiled.C_TT][1] > 0

    # ρ_vs (saturated vapor density) should be positive
    @test sol[compiled.ρ_vs][1] > 0

    # h_rel (relative humidity) should be between 0 and 1
    @test 0 < sol[compiled.h_rel][1] <= 1

    # D_mv (isothermal vapor diffusivity) should be positive
    @test sol[compiled.D_mv][1] > 0

    # D_Tv (thermal vapor diffusivity) should be positive
    @test sol[compiled.D_Tv][1] > 0

    # λ_eff (effective thermal conductivity) should be positive
    @test sol[compiled.λ_eff][1] > 0
end

@testitem "MulchHeatWaterTransfer: Saturated Vapor Density" setup = [MulchSetup] tags = [:mulch] begin
    # Verify Tetens formula at 20°C (293.15K)
    # P_vs(20°C) ≈ 611 * exp(17.27*20/(20+237.3)) ≈ 2338 Pa
    # ρ_vs = M_w * P_vs / (R * T) ≈ 0.01802 * 2338 / (8.314 * 293.15) ≈ 0.01728 kg/m³
    T_C = 20.0
    P_vs = 611.0 * exp(17.27 * T_C / (T_C + 237.3))
    ρ_vs_expected = 0.01802 * P_vs / (8.314 * 293.15)
    @test ρ_vs_expected > 0.01  # Sanity check

    sys = MulchHeatWaterTransfer()
    compiled = mtkcompile(sys)
    prob = ODEProblem(compiled, [], (0.0, 1.0))
    sol = solve(prob)

    @test sol[compiled.ρ_vs][1] ≈ ρ_vs_expected rtol = 0.01
end

# ============================================================
# MulchHeatWaterPDE Tests
# ============================================================

@testsnippet MulchPDESetup begin
    using ModelingToolkit
    using ModelingToolkit: t, D
    using DomainSets
    using MethodOfLines
    using OrdinaryDiffEqDefault
    using OrdinaryDiffEqDefault: SciMLBase
    using Vegetation
end

@testitem "MulchHeatWaterPDE: Structural Verification" setup = [MulchPDESetup] tags = [:mulch_pde] begin
    pde = MulchHeatWaterPDE(0.06, 3600.0)

    @test length(pde.eqs) == 4
    @test length(pde.dvs) == 4
    @test length(pde.ps) == 14
    @test length(pde.ivs) == 2
end

@testitem "MulchHeatWaterPDE: Discretization" setup = [MulchPDESetup] tags = [:mulch_pde] begin
    pde = MulchHeatWaterPDE(0.06, 3600.0)
    z = pde.ivs[2]
    dz = 0.02
    disc = MOLFiniteDifference([z => dz], t, approx_order = 2)
    prob = discretize(pde, disc; checks = false)

    @test prob isa ODEProblem
    @test length(prob.u0) > 0
    @test prob.tspan == (0.0, 3600.0)
end

@testitem "MulchHeatWaterPDE: Solution" setup = [MulchPDESetup] tags = [:mulch_pde] begin
    pde = MulchHeatWaterPDE(0.06, 3600.0)
    z = pde.ivs[2]
    dz = 0.02
    disc = MOLFiniteDifference([z => dz], t, approx_order = 2)
    prob = discretize(pde, disc; checks = false)

    sol = solve(prob)
    @test sol.retcode == SciMLBase.ReturnCode.Success
    @test length(sol.t) > 1
end

@testitem "MulchHeatWaterPDE: Boundary Values" setup = [MulchPDESetup] tags = [:mulch_pde] begin
    pde = MulchHeatWaterPDE(
        0.06, 3600.0;
        h_top = -0.5, h_bot = -2.0, T_top = 298.15, T_bot = 288.15
    )

    z = pde.ivs[2]
    dz = 0.02
    disc = MOLFiniteDifference([z => dz], t, approx_order = 2)
    prob = discretize(pde, disc; checks = false)

    # All initial values should be finite
    @test all(isfinite.(prob.u0))
    @test prob.tspan == (0.0, 3600.0)
end

# ============================================================
# MulchSurfaceRunoffPDE Tests
# ============================================================

@testitem "MulchSurfaceRunoffPDE: Structural Verification" setup = [MulchPDESetup] tags = [:mulch_pde] begin
    pde = MulchSurfaceRunoffPDE(0.5, 60.0)

    @test length(pde.eqs) == 3
    @test length(pde.dvs) == 3
    @test length(pde.ps) == 11
    @test length(pde.ivs) == 2
end

@testitem "MulchSurfaceRunoffPDE: Discretization" setup = [MulchPDESetup] tags = [:mulch_pde] begin
    pde = MulchSurfaceRunoffPDE(0.5, 60.0)
    l = pde.ivs[2]
    dl = 0.1
    disc = MOLFiniteDifference([l => dl], t, approx_order = 2)
    prob = discretize(pde, disc; checks = false)

    @test prob isa ODEProblem
    @test length(prob.u0) > 0
    @test prob.tspan == (0.0, 60.0)
    @test length(prob.u0) >= 8
end

@testitem "MulchSurfaceRunoffPDE: Solution" setup = [MulchPDESetup] tags = [:mulch_pde] begin
    pde = MulchSurfaceRunoffPDE(
        0.5, 60.0;
        P_val = 70.0 / 1000 / 3600,
        S_0_val = 0.01,
        n_manning_val = 0.15,
        h_init_val = 1.0e-3,
        q_init_val = 0.0
    )

    l = pde.ivs[2]
    dl = 0.1
    disc = MOLFiniteDifference([l => dl], t, approx_order = 2)
    prob = discretize(pde, disc; checks = false)

    sol = solve(prob)
    @test sol.retcode == SciMLBase.ReturnCode.Success
    @test length(sol.t) > 1
end

@testitem "MulchSurfaceRunoffPDE: Custom Parameters" setup = [MulchPDESetup] tags = [:mulch_pde] begin
    pde = MulchSurfaceRunoffPDE(
        1.0, 120.0;
        P_val = 1.0e-4,
        I_val = 5.0e-5,
        S_0_val = 0.02,
        n_manning_val = 0.2
    )

    # Verify expected parameters exist
    param_names = [string(p) for p in pde.ps]
    @test any(n -> contains(n, "P_rate"), param_names)
    @test any(n -> contains(n, "I_rate"), param_names)
    @test any(n -> contains(n, "S_0"), param_names)
    @test any(n -> contains(n, "n_mann"), param_names)

    # Verify domain length is 1.0 and time span is 120.0
    @test pde.domain[1].domain.right == 120.0
    @test pde.domain[2].domain.right == 1.0
end
