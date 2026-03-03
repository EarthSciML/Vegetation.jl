# Maize Diffusive Root Growth

## Overview

This module implements a diffusive root growth model for maize (*Zea mays* L.) based on
the MAIZSIM crop simulator. The model treats root carbon expansion as a diffusion process
in a 2D soil domain, with young roots diffusing spatially and maturing into non-diffusing
mature roots.

The model consists of three coupled processes:
1. **Carbon Allocation** — partitions photosynthetic carbon between shoot and root growth
2. **Carbon Assignment** — distributes root carbon to spatial locations based on favorability indices
3. **Root Diffusion** — solves the spatial distribution of young and mature root density

Four favorability indices control root growth at each location: penetration resistance (f₁),
temperature (f₂), aeration (f₃), and root density (f₄).

**Reference**: Wang, Z., Timlin, D., Li, S., Fleisher, D., Dathe, A., Luo, C., Dong, L.,
Reddy, V.R., and Tully, K. (2021). A diffusive model of maize root growth in MAIZSIM
and its applications in Ridge-Furrow Rainfall Harvesting. *Agricultural Water Management*,
254, 106966. [doi:10.1016/j.agwat.2021.106966](https://doi.org/10.1016/j.agwat.2021.106966)

```@docs
MaizeRootGrowth
```

## Implementation

The component implements Equations 1–4 from Wang et al. (2021). The non-spatial (ODE)
form is provided here; the full 2D PDE form (Eq. 3a) with spatial diffusion can be
constructed using MethodOfLines.jl.

```@example maize_root
using ModelingToolkit, Vegetation, DynamicQuantities, Symbolics, DataFrames

sys = MaizeRootGrowth()
nothing # hide
```

### State Variables

```@example maize_root
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example maize_root
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape=false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example maize_root
eqs = equations(sys)
```

## Analysis

### Root Growth Dynamics Under Default Conditions

The following example simulates root growth for 90 days under default soil conditions
(well-aerated, moderate bulk density, optimal temperature range). Young roots (Y)
accumulate from the carbon source, then mature into non-diffusing mature roots (M).

```@example maize_root
using OrdinaryDiffEqDefault, Plots

compiled = mtkcompile(sys)
tspan = (0.0, 90.0 * 86400.0)  # 90 days in seconds
prob = ODEProblem(compiled, [], tspan)
sol = solve(prob)

days = sol.t ./ 86400.0
p = plot(days, sol[compiled.Y] .* 1000, label="Young roots (Y)", xlabel="Time (days)",
    ylabel="Root density (g/m³)", title="Root Growth Dynamics (Wang et al. 2021)",
    linewidth=2)
plot!(p, days, sol[compiled.M] .* 1000, label="Mature roots (M)", linewidth=2)
plot!(p, days, (sol[compiled.Y] .+ sol[compiled.M]) .* 1000,
    label="Total (Y + M)", linewidth=2, linestyle=:dash)
p
```

### Effect of Temperature on Root Growth

Temperature favorability (f₂) follows a piecewise function: growth increases from 0 to 18°C,
is optimal between 18–33°C, and decreases above 33°C (Eq. 1).

```@example maize_root
T_vals = [283.0, 288.0, 293.0, 298.0, 303.0, 308.0, 313.0]  # K
T_labels = round.(T_vals .- 273.15, digits=0)

tspan = (0.0, 60.0 * 86400.0)
p = plot(xlabel="Time (days)", ylabel="Total root density (g/m³)",
    title="Temperature Effect on Root Growth", legend=:topleft)

for (T, lab) in zip(T_vals, T_labels)
    prob_T = ODEProblem(compiled, [], tspan, [compiled.T_soil => T])
    sol_T = solve(prob_T)
    days = sol_T.t ./ 86400.0
    total = (sol_T[compiled.Y] .+ sol_T[compiled.M]) .* 1000
    plot!(p, days, total, label="$(Int(lab))°C", linewidth=2)
end
p
```

### Effect of Soil Bulk Density on Root Growth

Higher bulk density increases penetration resistance, reducing f₁ and slowing root growth.
The model captures the observed effect of compacted soil layers inhibiting root expansion
(Laboski et al. 1998).

```@example maize_root
rho_vals = [1200.0, 1380.0, 1500.0, 1600.0]  # kg/m³

tspan = (0.0, 60.0 * 86400.0)
p = plot(xlabel="Time (days)", ylabel="Total root density (g/m³)",
    title="Bulk Density Effect on Root Growth", legend=:topleft)

for rho in rho_vals
    prob_rho = ODEProblem(compiled, [], tspan, [compiled.ρ_b => rho])
    sol_rho = solve(prob_rho)
    days = sol_rho.t ./ 86400.0
    total = (sol_rho[compiled.Y] .+ sol_rho[compiled.M]) .* 1000
    plot!(p, days, total, label="ρ_b = $(Int(rho)) kg/m³", linewidth=2)
end
p
```
