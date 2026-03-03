@testsnippet MaizeRootSetup begin
    using Test
    using ModelingToolkit
    using ModelingToolkit: t, D
    using OrdinaryDiffEqDefault
    using OrdinaryDiffEqDefault: SciMLBase
    using Vegetation

    # Unit conversion constants
    const one_day = 86400.0       # seconds
    const one_bar = 1.0e5           # Pa
    const one_cm2_day = 1.0e-4 / 86400.0  # m²/s
end

@testitem "MaizeRootGrowth: Structural Verification" setup = [MaizeRootSetup] tags = [:maize] begin
    sys = MaizeRootGrowth()
    @test sys isa ModelingToolkit.System
    @test nameof(sys) == :MaizeRootGrowth

    vars = unknowns(sys)
    eqs = equations(sys)

    # 2 state variables (Y, M) + 10 algebraic variables = 12 unknowns / equations
    @test length(vars) == 12
    @test length(eqs) == 12

    # Verify key variable names exist
    var_names = [string(v) for v in vars]
    for expected in [
            "Y(t)", "M(t)", "f1(t)", "f2(t)", "f3(t)", "f4(t)",
            "R_bar(t)", "D_eff_xx(t)", "D_eff_zz(t)",
        ]
        @test any(n -> contains(n, expected), var_names)
    end

    # Verify compilation succeeds
    compiled = mtkcompile(sys)
    @test compiled !== nothing
    @test length(unknowns(compiled)) == 2  # Y and M are the only ODE states
end

@testitem "MaizeRootGrowth: Parameter Defaults" setup = [MaizeRootSetup] tags = [:maize] begin
    sys = MaizeRootGrowth()
    compiled = mtkcompile(sys)

    params = parameters(compiled)
    pdict = Dict(
        Symbol(p) => ModelingToolkit.getdefault(p) for p in params
            if ModelingToolkit.hasdefault(p)
    )

    # Verify A_growth = 0.55/day in SI (s⁻¹)
    @test pdict[:A_growth] ≈ 0.55 / one_day rtol = 1.0e-6

    # Verify D0_xx = 50 cm²/day in SI (m²/s) — horizontal diffusivity
    @test pdict[:D0_xx] ≈ 50.0 * one_cm2_day rtol = 1.0e-6

    # Verify D0_zz = 3 cm²/day in SI (m²/s) — vertical diffusivity
    # Note: D0_zz is defined but not used in the ODE form (only needed for PDE),
    # so check it on the uncompiled system's parameters instead.
    sys_params = parameters(MaizeRootGrowth())
    sys_pdict = Dict(
        Symbol(p) => ModelingToolkit.getdefault(p) for p in sys_params
            if ModelingToolkit.hasdefault(p)
    )
    @test sys_pdict[:D0_zz] ≈ 3.0 * one_cm2_day rtol = 1.0e-6

    # Verify soil temperature default = 298 K
    @test pdict[:T_soil] ≈ 298.0 rtol = 1.0e-6

    # Verify bulk density default = 1380 kg/m³
    @test pdict[:ρ_b] ≈ 1380.0 rtol = 1.0e-6
end

@testitem "MaizeRootGrowth: Equation Verification - f₁" setup = [MaizeRootSetup] tags = [:maize] begin
    # Verify f₁ against hand-computed reference values from Eq. 1
    # f₁ = (1/2)(ψ_trd - 5.4|ψ|^0.25·exp(-10.58(1.7-ρ_b))) - (1/4)(ψ_trd - ψ)
    # where ψ_trd, ψ in bar and ρ_b in Mg/m³

    sys = MaizeRootGrowth()
    compiled = mtkcompile(sys)

    # Test case: default conditions (ψ_rtd=5bar, ψ_soil=-0.3bar, ρ_b=1.38 Mg/m³)
    tspan = (0.0, 1.0)  # very short, just evaluate
    prob = ODEProblem(compiled, Dict(), tspan)
    sol = solve(prob)

    # Hand-compute f₁ with defaults:
    # ψ_trd_bar = 5.0e5/1.0e5 = 5.0
    # ψ_soil_bar = -3.0e4/1.0e5 = -0.3
    # ρ_b_Mg = 1380/1000 = 1.38
    psi_trd = 5.0
    psi_s = -0.3
    rho = 1.38
    f1_expected = 0.5 * (psi_trd - 5.4 * abs(psi_s)^0.25 * exp(-10.58 * (1.7 - rho))) -
        0.25 * (psi_trd - psi_s)
    f1_expected = clamp(f1_expected, 0.0, 1.0)

    @test sol[compiled.f1][1] ≈ f1_expected rtol = 1.0e-6
end

@testitem "MaizeRootGrowth: Equation Verification - f₂" setup = [MaizeRootSetup] tags = [:maize] begin
    # Verify f₂ temperature favorability
    sys = MaizeRootGrowth()
    compiled = mtkcompile(sys)
    tspan = (0.0, 1.0)

    # At 25°C (298K): 18 ≤ 25 < 33, so f₂ = 1.0
    prob = ODEProblem(compiled, Dict(compiled.T_soil => 298.0), tspan)
    sol = solve(prob)
    @test sol[compiled.f2][1] ≈ 1.0 rtol = 1.0e-6

    # At 10°C (283.15K): f₂ = (10/18)^1.66
    prob_cold = ODEProblem(compiled, Dict(compiled.T_soil => 283.15), tspan)
    sol_cold = solve(prob_cold)
    f2_expected = (10.0 / 18.0)^1.66
    @test sol_cold[compiled.f2][1] ≈ f2_expected rtol = 1.0e-3

    # At 40°C (313.15K): f₂ = (40/33)^(-1.66)
    prob_hot = ODEProblem(compiled, Dict(compiled.T_soil => 313.15), tspan)
    sol_hot = solve(prob_hot)
    f2_expected_hot = (40.0 / 33.0)^(-1.66)
    @test sol_hot[compiled.f2][1] ≈ f2_expected_hot rtol = 1.0e-3
end

@testitem "MaizeRootGrowth: Equation Verification - f₃" setup = [MaizeRootSetup] tags = [:maize] begin
    # Verify f₃ aeration favorability
    # f₃ = ([O₂] - 0.02)^7.14 where [O₂] in mol/L
    sys = MaizeRootGrowth()
    compiled = mtkcompile(sys)
    tspan = (0.0, 1.0)

    # Default: O2_soil = 1020 mol/m³ = 1.02 mol/L → f₃ = (1.0)^7.14 = 1.0
    prob = ODEProblem(compiled, Dict(), tspan)
    sol = solve(prob)
    @test sol[compiled.f3][1] ≈ 1.0 rtol = 1.0e-6

    # Low O₂: 100 mol/m³ = 0.1 mol/L → f₃ = (0.08)^7.14
    prob_low = ODEProblem(compiled, Dict(compiled.O2_soil => 100.0), tspan)
    sol_low = solve(prob_low)
    f3_expected = (0.1 - 0.02)^7.14
    @test sol_low[compiled.f3][1] ≈ f3_expected rtol = 1.0e-3
end

@testitem "MaizeRootGrowth: Equation Verification - f₄" setup = [MaizeRootSetup] tags = [:maize] begin
    # Verify f₄ root density favorability
    # f₄ = 1 - min(1, (M+Y)/0.03)
    sys = MaizeRootGrowth()
    compiled = mtkcompile(sys)
    tspan = (0.0, 1.0)

    # Default: Y=0.001, M=0.0, total=0.001 → f₄ = 1 - 0.001/0.03 ≈ 0.9667
    prob = ODEProblem(compiled, Dict(), tspan)
    sol = solve(prob)
    f4_expected = 1.0 - 0.001 / 0.03
    @test sol[compiled.f4][1] ≈ f4_expected rtol = 1.0e-4

    # High density: Y=0.001, M=0.028, total=0.029 → f₄ = 1 - 0.029/0.03 ≈ 0.0333
    prob_high = ODEProblem(compiled, Dict(compiled.Y => 0.001, compiled.M => 0.028), tspan)
    sol_high = solve(prob_high)
    f4_high = 1.0 - 0.029 / 0.03
    @test sol_high[compiled.f4][1] ≈ f4_high rtol = 1.0e-3

    # Above threshold: Y=0.001, M=0.03, total=0.031 → f₄ = 0
    prob_over = ODEProblem(compiled, Dict(compiled.Y => 0.001, compiled.M => 0.03), tspan)
    sol_over = solve(prob_over)
    @test sol_over[compiled.f4][1] ≈ 0.0 atol = 1.0e-10
end

@testitem "MaizeRootGrowth: Equation Verification - f̃₁(ψ)" setup = [MaizeRootSetup] tags = [:maize] begin
    # Verify f̃₁ water potential factor for diffusion (Eq. 4)
    # f̃₁(ψ) = 0.5*sin(π*(ψ-(ψ_s+ψ_r)/2)/(ψ_s-ψ_r)) + 0.5
    # ψ_s = -150 cm, ψ_r = -500 cm, midpoint = -325 cm
    sys = MaizeRootGrowth()
    compiled = mtkcompile(sys)
    tspan = (0.0, 1.0)

    # At ψ_soil corresponding to -325 cm head (midpoint):
    # -325 cm × 98.0665 Pa/cm = -31871.6 Pa
    psi_mid_Pa = -325.0 * 98.0665
    prob_mid = ODEProblem(compiled, Dict(compiled.ψ_soil => psi_mid_Pa), tspan)
    sol_mid = solve(prob_mid)
    @test sol_mid[compiled.f_tilde_psi][1] ≈ 0.5 rtol = 1.0e-3

    # At ψ_soil corresponding to -150 cm head (wet limit, ψ_s):
    psi_wet_Pa = -150.0 * 98.0665
    prob_wet = ODEProblem(compiled, Dict(compiled.ψ_soil => psi_wet_Pa), tspan)
    sol_wet = solve(prob_wet)
    @test sol_wet[compiled.f_tilde_psi][1] ≈ 1.0 rtol = 1.0e-3

    # At ψ_soil corresponding to -500 cm head (dry limit, ψ_r):
    psi_dry_Pa = -500.0 * 98.0665
    prob_dry = ODEProblem(compiled, Dict(compiled.ψ_soil => psi_dry_Pa), tspan)
    sol_dry = solve(prob_dry)
    @test sol_dry[compiled.f_tilde_psi][1] ≈ 0.0 atol = 1.0e-3
end

@testitem "MaizeRootGrowth: Equation Verification - R̄ (Eq. 2)" setup = [MaizeRootSetup] tags = [:maize] begin
    # Verify R̄ = (M + Y) × A × min{f₁, f₂, f₃, f₄}
    sys = MaizeRootGrowth()
    compiled = mtkcompile(sys)
    tspan = (0.0, 1.0)

    # Under default optimal conditions (f₁≈1, f₂=1, f₃≈1, f₄≈0.967):
    # f_min ≈ f₄ ≈ 0.967
    # R̄ = (0.001 + 0.0) × (0.55/86400) × 0.967 ≈ 6.155e-9 kg/m³/s
    prob = ODEProblem(compiled, Dict(), tspan)
    sol = solve(prob)

    A_val = 0.55 / one_day
    f4_val = 1.0 - 0.001 / 0.03
    R_bar_expected = 0.001 * A_val * f4_val
    @test sol[compiled.R_bar][1] ≈ R_bar_expected rtol = 1.0e-3
end

@testitem "MaizeRootGrowth: Basic Integration" setup = [MaizeRootSetup] tags = [:maize] begin
    sys = MaizeRootGrowth()
    compiled = mtkcompile(sys)

    # 30-day simulation
    tspan = (0.0, 30.0 * one_day)
    prob = ODEProblem(compiled, Dict(), tspan)
    sol = solve(prob)

    @test sol.retcode == SciMLBase.ReturnCode.Success

    # Y and M should be non-negative throughout
    @test all(sol[compiled.Y] .>= -1.0e-15)
    @test all(sol[compiled.M] .>= -1.0e-15)
end

@testitem "MaizeRootGrowth: Maturation Conservation" setup = [MaizeRootSetup] tags = [:maize] begin
    # Test that total root carbon (Y + M) increases at the rate of carbon input R
    # and that the maturation term transfers from Y to M without loss
    sys = MaizeRootGrowth()
    compiled = mtkcompile(sys)

    tspan = (0.0, 10.0 * one_day)
    prob = ODEProblem(compiled, Dict(), tspan)
    sol = solve(prob)

    @test sol.retcode == SciMLBase.ReturnCode.Success

    # Total root density should increase monotonically (since R > 0 and no death)
    total = sol[compiled.Y] .+ sol[compiled.M]
    for i in 2:length(total)
        @test total[i] >= total[i - 1] - 1.0e-12
    end

    # At the end, total roots should exceed initial
    @test total[end] > total[1]
end

@testitem "MaizeRootGrowth: Steady State Y" setup = [MaizeRootSetup] tags = [:maize] begin
    # With constant R_total input, Y should approach a steady state where
    # R_input = T_YM * Y, so Y_ss = R_input / T_YM
    sys = MaizeRootGrowth()
    compiled = mtkcompile(sys)

    # Run for a long time so Y approaches steady state
    # Use a small R_total so f4 doesn't saturate too fast
    R_val = 1.0e-4 / one_day  # small carbon input
    T_YM_val = 0.1 / one_day
    tspan = (0.0, 100.0 * one_day)
    prob = ODEProblem(
        compiled,
        Dict(compiled.R_total => R_val, compiled.T_YM => T_YM_val),
        tspan
    )
    sol = solve(prob)

    @test sol.retcode == SciMLBase.ReturnCode.Success

    # Y should stabilize; check that dY/dt is small relative to Y at the end
    Y_end = sol[compiled.Y][end]
    Y_prev = sol[compiled.Y][end - 1]
    dt = sol.t[end] - sol.t[end - 1]
    dYdt = (Y_end - Y_prev) / dt

    # dY/dt should be small relative to T_YM * Y
    @test abs(dYdt) < 0.1 * T_YM_val * Y_end || Y_end < 1.0e-10
end

@testitem "MaizeRootGrowth: Temperature Favorability" setup = [MaizeRootSetup] tags = [:maize] begin
    # Test f₂ behavior:
    # - f₂ = 1 for T in [18, 33] °C
    # - f₂ < 1 for T < 18 or T > 33
    sys = MaizeRootGrowth()
    compiled = mtkcompile(sys)

    tspan = (0.0, 1.0 * one_day)

    # At 25°C (298K), f2 should be 1.0
    prob_optimal = ODEProblem(
        compiled,
        Dict(compiled.T_soil => 298.0),
        tspan
    )
    sol_optimal = solve(prob_optimal)
    @test sol_optimal.retcode == SciMLBase.ReturnCode.Success

    # At 10°C (283K), f2 should be < 1
    prob_cold = ODEProblem(
        compiled,
        Dict(compiled.T_soil => 283.0),
        tspan
    )
    sol_cold = solve(prob_cold)
    @test sol_cold.retcode == SciMLBase.ReturnCode.Success

    # Root growth should be faster at optimal temperature
    total_optimal = sol_optimal[compiled.Y][end] + sol_optimal[compiled.M][end]
    total_cold = sol_cold[compiled.Y][end] + sol_cold[compiled.M][end]
    @test total_optimal > total_cold
end

@testitem "MaizeRootGrowth: High Root Density Limits Growth" setup = [MaizeRootSetup] tags = [:maize] begin
    # When M + Y approaches 0.03 kg/m³ threshold, f4 → 0 and growth should slow
    sys = MaizeRootGrowth()
    compiled = mtkcompile(sys)

    # Start with high initial root density near the threshold
    tspan = (0.0, 10.0 * one_day)
    prob_high = ODEProblem(
        compiled,
        Dict(compiled.Y => 0.001, compiled.M => 0.028),
        tspan
    )
    sol_high = solve(prob_high)

    # Start with low initial root density
    prob_low = ODEProblem(
        compiled,
        Dict(compiled.Y => 0.001, compiled.M => 0.0),
        tspan
    )
    sol_low = solve(prob_low)

    @test sol_high.retcode == SciMLBase.ReturnCode.Success
    @test sol_low.retcode == SciMLBase.ReturnCode.Success

    # Growth rate should be slower when starting near threshold
    growth_high = (sol_high[compiled.Y][end] + sol_high[compiled.M][end]) -
        (sol_high[compiled.Y][1] + sol_high[compiled.M][1])
    growth_low = (sol_low[compiled.Y][end] + sol_low[compiled.M][end]) -
        (sol_low[compiled.Y][1] + sol_low[compiled.M][1])
    @test growth_low > growth_high
end

@testitem "MaizeRootGrowth: Positivity Preservation" setup = [MaizeRootSetup] tags = [:maize] begin
    # Y and M should remain non-negative for various initial conditions
    sys = MaizeRootGrowth()
    compiled = mtkcompile(sys)

    tspan = (0.0, 60.0 * one_day)

    # Test with very small initial Y
    prob = ODEProblem(
        compiled,
        Dict(compiled.Y => 1.0e-8, compiled.M => 0.0),
        tspan
    )
    sol = solve(prob)

    @test sol.retcode == SciMLBase.ReturnCode.Success
    @test all(sol[compiled.Y] .>= -1.0e-15)
    @test all(sol[compiled.M] .>= -1.0e-15)
end

@testitem "MaizeRootGrowth: Diffusion Factors (Eq. 4)" setup = [MaizeRootSetup] tags = [:maize] begin
    # Verify diffusion factors f̃₁ and f̃₂ from Eq. 4
    # D_eff = D⁰ × min(f̃₁, f̃₂) — D_eff itself is only relevant in the PDE form,
    # but the factors are observable from the compiled ODE system.
    sys = MaizeRootGrowth()
    compiled = mtkcompile(sys)
    tspan = (0.0, 1.0)

    # Under default conditions, f̃₂ ≈ 1.0 (temperature factor is essentially always 1)
    prob = ODEProblem(compiled, Dict(), tspan)
    sol = solve(prob)
    @test sol[compiled.f_tilde_T][1] ≈ 1.0 rtol = 1.0e-6

    # Under default ψ_soil = -3e4 Pa = ~-306 cm head (between ψ_s and ψ_r),
    # f̃₁ should be between 0 and 1
    @test 0.0 < sol[compiled.f_tilde_psi][1] < 1.0

    # At the dry limit (ψ_r = -500 cm), f̃₁ ≈ 0
    psi_dry_Pa = -500.0 * 98.0665
    prob_dry = ODEProblem(compiled, Dict(compiled.ψ_soil => psi_dry_Pa), tspan)
    sol_dry = solve(prob_dry)
    @test sol_dry[compiled.f_tilde_psi][1] ≈ 0.0 atol = 1.0e-6

    # At the wet limit (ψ_s = -150 cm), f̃₁ ≈ 1
    psi_wet_Pa = -150.0 * 98.0665
    prob_wet = ODEProblem(compiled, Dict(compiled.ψ_soil => psi_wet_Pa), tspan)
    sol_wet = solve(prob_wet)
    @test sol_wet[compiled.f_tilde_psi][1] ≈ 1.0 atol = 1.0e-6
end
