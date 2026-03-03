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

### Equations

**Equation 1** — Four favorability indices (Acock et al. 1985):

- ``f_1 = \frac{1}{2}(\psi_{trd} - 5.4|\psi|^{0.25} \exp(-10.58(1.7 - \rho_b))) - \frac{1}{4}(\psi_{trd} - \psi)``
- ``f_2 = (T/18)^{1.66}`` for ``0 < T < 18``°C; ``f_2 = 1`` for ``18 \leq T < 33``°C; ``f_2 = (T/33)^{-1.66}`` for ``T \geq 33``°C
- ``f_3 = ([O_2] - 0.02)^{7.14}``
- ``f_4 = 1 - \min(1, (M + Y)/0.03)``

**Equation 2** — Potential carbon assigned for root growth:

``\bar{R} = (M + Y) \times A \times \min\{f_1, f_2, f_3, f_4\}``

**Equations 3a, 3b** — Governing equations for diffusive root growth:

``\frac{\partial Y}{\partial t} = \nabla \cdot [D \nabla Y] + R - T_{Y \to M} Y``

``\frac{\partial M}{\partial t} = T_{Y \to M} Y``

**Equation 4** — Diffusion coefficient:

``D = D^0 \times \min\{\tilde{f}_1(\psi), \tilde{f}_2(T)\}``

where ``\tilde{f}_1(\psi) = \frac{1}{2}\sin\left(\frac{\pi(\psi - (\psi_s + \psi_r)/2)}{\psi_s - \psi_r}\right) + \frac{1}{2}`` and ``\tilde{f}_2(T) = \max\left\{\frac{(1 + e^{p - T/T_0}) e^{T/T_0 - p}}{1 + e^{q - u/T}}, 1.0\right\}``

### Soil Properties (Table 1 from Wang et al. 2021)

The following table reproduces the layered soil physical properties used in validation
simulations (Sections 3.2 and 4 of the paper, from Zhao et al. 2018):

| Property | Layer 1 (<20 cm) | Layer 2 (20–35 cm) | Layer 3 (35–55 cm) | Layer 4 (55–75 cm) | Layer 5 (>75 cm) | Fine Textured Soil Mulch (surface 3 cm on ridge) |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
| Residual Water Content θ_r (cm³/cm³) | 0.05 | 0.05 | 0.05 | 0.05 | 0.04 | 0.20 |
| Saturated Water Content θ_s (cm³/cm³) | 0.39 | 0.39 | 0.39 | 0.38 | 0.38 | 0.52 |
| van Genuchten α | 0.02 | 0.02 | 0.01 | 0.01 | 0.01 | 0.03 |
| van Genuchten n | 1.30 | 1.26 | 1.28 | 1.20 | 1.38 | 1.10 |
| Sat. Hydraulic Conductivity K_s (cm/day) | 140.30 | 144.70 | 155.80 | 170.50 | 122.60 | 1.00 |
| Soil Bulk Density ρ_b (g/cm³) | 1.38 | 1.37 | 1.39 | 1.38 | 1.44 | 1.20 |
| Soil Organic Matter (g/g) | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| Sand Fraction (g/g) | 0.40 | 0.40 | 0.42 | 0.38 | 0.49 | 0.25 |
| Silt Fraction (g/g) | 0.36 | 0.36 | 0.39 | 0.37 | 0.33 | 0.15 |

### Maize Yield and Water Balance (Table 2 from Wang et al. 2021)

The following table reproduces the simulated maize growth and evaporation-transpiration
results for the example in Section 3.2 of the paper:

| Treatment | Yield (kg/ha) | Total Dry Mass (kg/ha) | Shoot Dry Mass (kg/ha) | Root Dry Mass (kg/ha) | Root/Shoot Ratio | Evaporation (mm/cm²) | Transpiration (mm/cm²) |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Observed Precipitation | 9138 | 23,563 | 12,760 | 10,802 | 0.85 | 350 | 344 |
| Observed Precip. + Irrigation | 12,209 | 24,029 | 16,017 | 8011 | 0.50 | 419 | 456 |
| 50% Observed Precipitation | 6689 | 23,215 | 9555 | 13,660 | 1.43 | 238 | 167 |

### RFRH Management Results (Table 3 from Wang et al. 2021)

The following table reproduces Table 3 from the paper, showing the simulated maize growth
and evaporation-transpiration for Ridge-Furrow Rainfall Harvest (RFRH) management compared
to flat soil surface. These results were produced using the full 2D PDE diffusive root model
coupled with the MAIZSIM simulator and cannot be reproduced with the local ODE form implemented here.

| Soil Surface | Treatment | Yield (kg/ha) | Total Dry Mass (kg/ha) | Shoot Dry Mass (kg/ha) | Root Dry Mass (kg/ha) | Root/Shoot | Evaporation (mm) | Transpiration (mm) |
|:---|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Flat Soil | Observed Precip. | 9138 | 23,563 | 12,760 | 10,802 | 0.85 | 350 | 344 |
| RFRH, No Cover | No Cover | 9256 | 27,485 | 13,884 | 13,601 | 0.98 | 125 | 294 |
| RFRH, Plastic | Plastic Cover | 9334 | 26,437 | 14,817 | 11,620 | 0.78 | 76 | 320 |
| RFRH, Fine Soil | Fine-texture Soil Cover | 9426 | 25,821 | 14,471 | 11,350 | 0.78 | 109 | 321 |

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

### System Equations

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
prob = ODEProblem(compiled, Dict(), tspan)
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

### Temperature Favorability f₂ Response Curve (cf. Eq. 1)

Temperature favorability (f₂) follows a piecewise function from Eq. 1: growth increases
from 0 to 18°C as ``(T/18)^{1.66}``, is optimal (f₂ = 1) between 18–33°C, and decreases
above 33°C as ``(T/33)^{-1.66}``. The following figure shows the f₂ curve directly as a
function of temperature, matching the functional form in the paper.

```@example maize_root
T_celsius = 0.5:0.5:45.0
T_kelvin = T_celsius .+ 273.15
f2_vals = Float64[]
tspan_short = (0.0, 1.0)
for Tk in T_kelvin
    prob_f2 = ODEProblem(compiled, Dict(compiled.T_soil => Tk), tspan_short)
    sol_f2 = solve(prob_f2)
    push!(f2_vals, sol_f2[compiled.f2][1])
end

p = plot(T_celsius, f2_vals, xlabel="Temperature (°C)", ylabel="f₂",
    title="Temperature Favorability Index (Eq. 1, f₂)", linewidth=2, legend=false)
vline!(p, [18.0, 33.0], linestyle=:dash, color=:gray)
p
```

### Effect of Temperature on Root Growth (cf. Eq. 1, f₂)

The following figure shows the effect of different soil temperatures on total root
growth over 60 days, demonstrating how the f₂ index controls growth rates.

```@example maize_root
T_vals = [283.0, 288.0, 293.0, 298.0, 303.0, 308.0, 313.0]  # K
T_labels = round.(T_vals .- 273.15, digits=0)

tspan = (0.0, 60.0 * 86400.0)
p = plot(xlabel="Time (days)", ylabel="Total root density (g/m³)",
    title="Temperature Effect on Root Growth (Eq. 1, f₂)", legend=:topleft)

for (T, lab) in zip(T_vals, T_labels)
    prob_T = ODEProblem(compiled, Dict(compiled.T_soil => T), tspan)
    sol_T = solve(prob_T)
    days = sol_T.t ./ 86400.0
    total = (sol_T[compiled.Y] .+ sol_T[compiled.M]) .* 1000
    plot!(p, days, total, label="$(Int(lab))°C", linewidth=2)
end
p
```

### Effect of Soil Bulk Density on Root Growth (cf. Eq. 1, f₁)

Higher bulk density increases penetration resistance, reducing f₁ and slowing root growth.
The model captures the observed effect of compacted soil layers inhibiting root expansion
(Laboski et al. 1998). The paper validated this behavior against root measurements in
Minnesota soil with bulk densities up to 1.58 Mg/m³ (Section 3.1).

```@example maize_root
rho_vals = [1200.0, 1380.0, 1500.0, 1600.0]  # kg/m³

tspan = (0.0, 60.0 * 86400.0)
p = plot(xlabel="Time (days)", ylabel="Total root density (g/m³)",
    title="Bulk Density Effect on Root Growth (Eq. 1, f₁)", legend=:topleft)

for rho in rho_vals
    prob_rho = ODEProblem(compiled, Dict(compiled.ρ_b => rho), tspan)
    sol_rho = solve(prob_rho)
    days = sol_rho.t ./ 86400.0
    total = (sol_rho[compiled.Y] .+ sol_rho[compiled.M]) .* 1000
    plot!(p, days, total, label="ρ_b = $(Int(rho)) kg/m³", linewidth=2)
end
p
```

### Effect of Root Density on Growth Rate (cf. Eq. 1, f₄)

The root density favorability index f₄ = 1 - min(1, (M+Y)/0.03) decreases linearly
as root density approaches 0.03 kg/m³, providing self-limitation. This prevents
unrealistic root accumulation at any single location.

```@example maize_root
M_init_vals = [0.0, 0.01, 0.02, 0.025, 0.029]

tspan = (0.0, 30.0 * 86400.0)
p = plot(xlabel="Time (days)", ylabel="Total root density (g/m³)",
    title="Root Density Self-Limitation (Eq. 1, f₄)", legend=:topleft)

for M0 in M_init_vals
    prob_m = ODEProblem(compiled, Dict(compiled.Y => 0.001, compiled.M => M0), tspan)
    sol_m = solve(prob_m)
    days = sol_m.t ./ 86400.0
    total = (sol_m[compiled.Y] .+ sol_m[compiled.M]) .* 1000
    plot!(p, days, total, label="M₀ = $(M0*1000) g/m³", linewidth=2)
end
p
```

### Diffusion Coefficient Water Potential Response (cf. Eq. 4, f̃₁)

The water potential factor f̃₁(ψ) for the diffusion coefficient maps the soil
water potential from the dry limit (ψ_r = -500 cm head) to the wet limit
(ψ_s = -150 cm head) via a sinusoidal function, giving f̃₁ = 0 at the dry limit
and f̃₁ = 1 at the wet limit.

```@example maize_root
psi_cm_vals = -600:5:-50  # cm water head
psi_Pa_vals = psi_cm_vals .* 98.0665  # Convert to Pa

f_psi_vals = Float64[]
tspan = (0.0, 1.0)
for psi_Pa in psi_Pa_vals
    prob_p = ODEProblem(compiled, Dict(compiled.ψ_soil => psi_Pa), tspan)
    sol_p = solve(prob_p)
    push!(f_psi_vals, sol_p[compiled.f_tilde_psi][1])
end

p = plot(psi_cm_vals, f_psi_vals, xlabel="Soil Water Potential (cm head)",
    ylabel="f̃₁(ψ)", title="Diffusion Water Potential Factor (Eq. 4)",
    linewidth=2, legend=false)
vline!(p, [-150, -500], linestyle=:dash, color=:gray, label="ψ_s, ψ_r")
p
```
