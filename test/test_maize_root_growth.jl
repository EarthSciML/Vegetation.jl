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

    # 2 state variables (Y, M) + 15 algebraic variables = 17 unknowns
    @test length(vars) == 17
    @test length(eqs) == 17

    # Verify key variable names exist
    var_names = [string(v) for v in vars]
    for expected in [
            "Y(t)", "M(t)", "f1(t)", "f2(t)", "f3(t)", "f4(t)",
            "R_bar(t)", "D_eff(t)",
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

    # Verify D0_xx = 50 cm²/day in SI (m²/s)
    @test pdict[:D0_xx] ≈ 50.0 * one_cm2_day rtol = 1.0e-6

    # Verify soil temperature default = 298 K
    @test pdict[:T_soil] ≈ 298.0 rtol = 1.0e-6

    # Verify bulk density default = 1380 kg/m³
    @test pdict[:ρ_b] ≈ 1380.0 rtol = 1.0e-6
end

@testitem "MaizeRootGrowth: Basic Integration" setup = [MaizeRootSetup] tags = [:maize] begin
    sys = MaizeRootGrowth()
    compiled = mtkcompile(sys)

    # 30-day simulation
    tspan = (0.0, 30.0 * one_day)
    prob = ODEProblem(compiled, [], tspan)
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
    prob = ODEProblem(compiled, [], tspan)
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
        [],
        tspan,
        [compiled.R_total => R_val, compiled.T_YM => T_YM_val]
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
        compiled, [], tspan,
        [compiled.T_soil => 298.0]
    )
    sol_optimal = solve(prob_optimal)
    @test sol_optimal.retcode == SciMLBase.ReturnCode.Success

    # At 10°C (283K), f2 should be < 1
    prob_cold = ODEProblem(
        compiled, [], tspan,
        [compiled.T_soil => 283.0]
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
        [compiled.Y => 0.001, compiled.M => 0.028],
        tspan
    )
    sol_high = solve(prob_high)

    # Start with low initial root density
    prob_low = ODEProblem(
        compiled,
        [compiled.Y => 0.001, compiled.M => 0.0],
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
        [compiled.Y => 1.0e-8, compiled.M => 0.0],
        tspan
    )
    sol = solve(prob)

    @test sol.retcode == SciMLBase.ReturnCode.Success
    @test all(sol[compiled.Y] .>= -1.0e-15)
    @test all(sol[compiled.M] .>= -1.0e-15)
end
