# Residue Mulch Water and Thermal Dynamics

## Overview

This module implements the residue mulch decomposition, radiation attenuation, wind profile,
heat/vapor flux, and water characteristic models from Wang et al. (2021). The model simulates
water and thermal dynamics at soil surfaces with residue mulch and surface runoff.

Residue mulch influences soil surface processes by partially blocking solar radiation,
conserving soil water, and mitigating soil temperature variations. The model tracks three
residue mass pools (carbohydrate CARB, cellulose CELL, lignin LIGN) and their nitrogen
dynamics, including decomposition, N mineralization, immobilization, and humification.

**Reference**: Wang, Z., Thapa, R., Timlin, D., Li, S., Sun, W., Beegum, S., et al.
(2021). Simulations of water and thermal dynamics for soil surfaces with residue mulch
and surface runoff. *Water Resources Research*, 57, e2021WR030431.
[https://doi.org/10.1029/2021WR030431](https://doi.org/10.1029/2021WR030431)

```@docs
ResidueMulchDecomposition
MulchRadiationAttenuation
MulchWindProfile
MulchHeatVaporFluxes
MulchWaterCharacteristic
MulchHeatWaterTransfer
MulchHeatWaterPDE
MulchSurfaceRunoffPDE
```

## Implementation

The implementation consists of eight ModelingToolkit components:

### ODE Components
1. **`ResidueMulchDecomposition`**: Tracks residue mass and N pools with decomposition kinetics (Eq. 18-25)
2. **`MulchRadiationAttenuation`**: Computes shortwave and longwave radiation through mulch layers (Eq. 11-15)
3. **`MulchWindProfile`**: Models wind speed attenuation through the mulch (Eq. 3)
4. **`MulchHeatVaporFluxes`**: Parameterizes diffusive and convective heat/vapor fluxes (Eq. 4-8)
5. **`MulchWaterCharacteristic`**: Relates matric potential to water content in mulch (Eq. 26)
6. **`MulchHeatWaterTransfer`**: Constitutive relations and single-node ODE for coupled heat and water transfer (Eq. 1-2, 6, 9)

### PDE Systems
7. **`MulchHeatWaterPDE`**: Coupled heat and water transfer through mulch (Eq. 1), discretized with MethodOfLines.jl
8. **`MulchSurfaceRunoffPDE`**: Saint-Venant surface runoff equations (Eq. 16-17), discretized with MethodOfLines.jl

### State Variables (Decomposition Model)

```@example mulch
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using Vegetation

sys = ResidueMulchDecomposition()
compiled = mtkcompile(sys)

vars = unknowns(compiled)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters (Decomposition Model)

```@example mulch
params = parameters(compiled)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape=false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example mulch
eqs = equations(sys)
```

## Analysis

### Decomposition Over 100 Days (cf. Figure 11)

The following example reproduces the qualitative behavior shown in Figure 11 of the paper,
demonstrating how CARB decomposes fastest, followed by CELL, while LIGN is most resistant.

```@example mulch
using OrdinaryDiffEqDefault
using Plots

tspan = (0.0, 100.0 * 86400.0)  # 100 days in seconds

prob = ODEProblem(compiled, [], tspan)
sol = solve(prob)

days = sol.t ./ 86400.0

# Convert from kg/m² to g/m²
RM_CARB = sol[compiled.RM_CARB] .* 1000
RM_CELL = sol[compiled.RM_CELL] .* 1000
RM_LIGN = sol[compiled.RM_LIGN] .* 1000

p1 = plot(days, RM_CARB, label="CARB", xlabel="Day", ylabel="Residue Mass (g/m²)",
    title="(a) Residue Mass Pools", linewidth=2)
plot!(p1, days, RM_CELL, label="CELL", linewidth=2)
plot!(p1, days, RM_LIGN, label="LIGN", linewidth=2, linestyle=:dash)

# Fraction remaining (cf. Figure 11c)
frac_CARB = sol[compiled.RM_CARB] ./ sol[compiled.RM_CARB][1]
frac_CELL = sol[compiled.RM_CELL] ./ sol[compiled.RM_CELL][1]
frac_LIGN = sol[compiled.RM_LIGN] ./ sol[compiled.RM_LIGN][1]

p2 = plot(days, frac_CARB, label="CARB", xlabel="Day", ylabel="Fraction Remaining",
    title="(c) Fraction Remaining", linewidth=2, ylim=(0, 1.2))
plot!(p2, days, frac_CELL, label="CELL", linewidth=2)
plot!(p2, days, frac_LIGN, label="LIGN", linewidth=2, linestyle=:dash)

plot(p1, p2, layout=(1, 2), size=(900, 400))
```

### N Mass Pools (cf. Figure 11d)

```@example mulch
RMN_CARB = sol[compiled.RMN_CARB] .* 1000
RMN_CELL = sol[compiled.RMN_CELL] .* 1000
RMN_LIGN = sol[compiled.RMN_LIGN] .* 1000

p3 = plot(days, RMN_CARB, label="CARB N", xlabel="Day", ylabel="N Mass (g/m²)",
    title="(d) N Mass Pools", linewidth=2)
plot!(p3, days, RMN_CELL, label="CELL N", linewidth=2)
plot!(p3, days, RMN_LIGN, label="LIGN N", linewidth=2, linestyle=:dash)

# Mineral N accumulation
N_inorg = sol[compiled.N_inorg] .* 1000
p4 = plot(days, N_inorg, label="Mineral N", xlabel="Day", ylabel="Mineral N (g/m²)",
    title="Mineral N Accumulation", linewidth=2, color=:red)

plot(p3, p4, layout=(1, 2), size=(900, 400))
```

### Temperature Sensitivity of Decomposition

The decomposition rate depends strongly on temperature through the MTRF factor (Eq. 22).
Below zero, decomposition stops entirely.

```@example mulch
temps = [278.15, 283.15, 288.15, 293.15, 298.15, 303.15]  # 5°C to 30°C
labels = ["5°C", "10°C", "15°C", "20°C", "25°C", "30°C"]

p5 = plot(xlabel="Day", ylabel="Total Residue Mass (g/m²)",
    title="Temperature Sensitivity of Decomposition", legend=:topright)

for (i, T) in enumerate(temps)
    prob_T = ODEProblem(compiled, [compiled.T_mulch => T], tspan)
    sol_T = solve(prob_T)
    RM_total = (sol_T[compiled.RM_CARB] .+ sol_T[compiled.RM_CELL] .+ sol_T[compiled.RM_LIGN]) .* 1000
    plot!(p5, sol_T.t ./ 86400, RM_total, label=labels[i], linewidth=2)
end

p5
```

### Shortwave Radiation Attenuation Through Mulch Layers

Demonstrates the geometric attenuation of shortwave radiation (Eq. 11) through
5 mulch layers with default properties.

```@example mulch
K = 5
ΔR = 0.3     # residue-area index
Ω_cl = 0.6   # clumping index
S_0 = 500.0  # incoming solar (W/m²)

# Eq. 11: S_d at each interface
S_d = zeros(K + 1)
S_d[K+1] = S_0
S_d[K] = S_0 * (1 - ΔR)
for k in K-1:-1:1
    S_d[k] = S_0 * (1 - ΔR) * (1 - Ω_cl * ΔR)^(K - k)
end

z_interfaces = range(0, 6, length=K+1)  # cm

p6 = plot(S_d, z_interfaces, xlabel="Shortwave Radiation (W/m²)",
    ylabel="Height in Mulch (cm)", title="Shortwave Attenuation (Eq. 11)",
    marker=:circle, linewidth=2, label="S_d", legend=:bottomright)

p6
```

### Wind Speed Profile in Mulch

Shows the exponential attenuation of wind speed from the mulch-air interface
to the mulch-soil interface (Eq. 3).

```@example mulch
Z_M = 0.06   # mulch thickness (m)
z_ref = 2.0  # reference height (m)
u_ref = 3.0  # wind speed at reference (m/s)
k_K = 0.4

d_disp = 0.87 * Z_M
z_r = 0.079 * Z_M

u_star = k_K * u_ref / log((z_ref - d_disp) / z_r)

z_centers = [(2k - 1) / (2K) * Z_M for k in 1:K]
u_layers = [0.21 * u_star * exp(2.2 / Z_M * z) for z in z_centers]

p7 = plot(u_layers, z_centers .* 100, xlabel="Wind Speed (m/s)",
    ylabel="Height in Mulch (cm)", title="Wind Speed Profile (Eq. 3)",
    marker=:circle, linewidth=2, label="u(z)", legend=:bottomright)

p7
```

### Heat and Water Transfer Constitutive Relations

The `MulchHeatWaterTransfer` component implements the constitutive relations for
coupled heat and water transfer through mulch (Eq. 1-2, 6, 9), including the water
characteristic function inverse, vapor transport coefficients, and effective thermal
conductivity.

```@example mulch
hw_sys = MulchHeatWaterTransfer()
hw_compiled = mtkcompile(hw_sys)

vars = unknowns(hw_compiled)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example mulch
params = parameters(hw_compiled)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape=false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### PDE Systems

The `MulchHeatWaterPDE` creates a `PDESystem` for coupled heat and water transfer
(Eq. 1), while `MulchSurfaceRunoffPDE` creates a `PDESystem` for the Saint-Venant
surface runoff equations (Eq. 16-17). Both systems are suitable for spatial
discretization with MethodOfLines.jl.

```@example mulch
using DomainSets

pde_hw = MulchHeatWaterPDE(0.06, 3600.0)
println("Heat-Water PDE: ", length(pde_hw.eqs), " equations, ",
    length(pde_hw.dvs), " dependent variables, ",
    length(pde_hw.ps), " parameters")

pde_sr = MulchSurfaceRunoffPDE(0.5, 60.0)
println("Surface Runoff PDE: ", length(pde_sr.eqs), " equations, ",
    length(pde_sr.dvs), " dependent variables, ",
    length(pde_sr.ps), " parameters")
```
