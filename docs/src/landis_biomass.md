# LANDIS Biomass Module

## Overview

The LANDIS biomass module models forest growth, mortality, and dead woody biomass
decomposition for individual species-age cohorts. It integrates aboveground net
primary productivity (ANPP), biomass-related mortality, age-related mortality, and
dead wood decomposition into a cohort-level biomass dynamics model suitable for
landscape-scale forest simulations.

The model tracks two state variables per cohort: living aboveground biomass and
dead woody biomass. Growth follows a peaked function that increases at low biomass
and decreases as the cohort approaches its maximum potential biomass. Mortality
combines a logistic biomass-dependent term (increasing as the site fills up) and
an exponential age-dependent term (increasing sharply as the cohort approaches its
maximum lifespan). Dead biomass decays exponentially.

**Reference**: Scheller, R.M. and Mladenoff, D.J. (2004). A forest growth and
biomass module for a landscape simulation model, LANDIS: design, validation,
and application. *Ecological Modelling*, 180, 211-229.
doi:10.1016/j.ecolmodel.2004.01.022

```@docs
LANDISBiomass
```

## Implementation

### State Variables

```@example landis
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using Vegetation

sys = LANDISBiomass()
compiled = mtkcompile(sys)
vars = unknowns(compiled)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [string(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars],
)
```

### Parameters

```@example landis
params = parameters(compiled)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape=false)) for p in params],
    :Units => [string(ModelingToolkit.get_unit(p)) for p in params],
    :Default => [ModelingToolkit.hasdefault(p) ? ModelingToolkit.getdefault(p) : missing for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params],
)
```

### Equations

```@example landis
eqs = equations(sys)
```

## Analysis

### Fig. 3a: ANPP and Mortality vs. Actual/Potential Biomass Ratio

This figure reproduces Fig. 3a from the paper, showing the correlation between
ANPP (as a fraction of ANPP\_MAX) and biomass-related mortality (M\_BIO as a
fraction of ANPP\_MAX) as a function of B\_AP (the ratio of actual to potential
biomass). The ANPP curve peaks at B\_AP = 1 while mortality increases
logistically.

```@example landis
using Plots

B_AP_range = 0.0:0.01:1.0

# Eq. 4: ANPP fraction = e * B_AP * exp(-B_AP) (with B_PM = 1)
anpp_frac = [exp(1) * x * exp(-x) for x in B_AP_range]

# Eq. 5: Mortality fraction (r=0.08, y0=0.01)
r, y0 = 0.08, 0.01
mort_frac = [y0 / (y0 + (1 - y0) * exp(-r / y0 * x)) for x in B_AP_range]

p = plot(B_AP_range, anpp_frac, label="ANPP", linewidth=2,
    xlabel="B_AP (Actual biomass / Potential biomass)",
    ylabel="Fraction ANPP_MAX",
    title="Fig. 3a: Growth and Mortality vs. B_AP",
    legend=:right, ylim=(0, 1.05))
plot!(p, B_AP_range, mort_frac, label="Mortality", linewidth=2, linestyle=:dash)
p
```

### Fig. 3b: Age-Related Mortality

This figure reproduces Fig. 3b from the paper, showing the fraction of biomass
removed by age-related mortality as a function of the fraction of species
lifespan. With the default shape parameter d = 10, mortality begins to increase
noticeably around 50% of lifespan and reaches 100% at the maximum lifespan.

```@example landis
age_frac_range = 0.0:0.01:1.0
d = 10.0

# Eq. 6: M_AGE fraction = exp(age_frac * d) / exp(d)
age_mort_frac = [exp(x * d) / exp(d) for x in age_frac_range]

p = plot(age_frac_range, age_mort_frac, linewidth=2, label=nothing,
    xlabel="Fraction of species' lifespan",
    ylabel="Fraction of biomass removed",
    title="Fig. 3b: Age-Related Mortality",
    ylim=(0, 1.05))
p
```

### Fig. 6a: A. saccharum Single-Species Growth

This figure reproduces Fig. 6a from the paper, showing the growth trajectory of
a single *Acer saccharum* (sugar maple) cohort with default parameters
(ANPP\_MAX = 7.45 Mg/ha/yr, max\_age = 400 years). Living biomass reaches an
asymptote near 240 Mg/ha, while dead woody biomass accumulates over time. The
paper reports that the mortality rate does not exceed the decomposition rate
until around year 50, after which dead biomass remains roughly constant near
25 Mg/ha.

```@example landis
using OrdinaryDiffEqDefault

compiled = mtkcompile(LANDISBiomass())

yr_to_s = 3.15576e7
Mg_ha_to_kg_m2 = 0.1

tspan_s = 200.0 * yr_to_s
prob = ODEProblem(compiled, [], (0.0, tspan_s))
sol = solve(prob)

years = 1:200
B_vals = [sol(yr * yr_to_s; idxs=compiled.B) / Mg_ha_to_kg_m2 for yr in years]
D_vals = [sol(yr * yr_to_s; idxs=compiled.D_wood) / Mg_ha_to_kg_m2 for yr in years]

p = plot(years, B_vals, label="A. saccharum", linewidth=2,
    xlabel="Simulation year",
    ylabel="Biomass (Mg ha⁻¹)",
    title="Fig. 6a: Single Species Growth (A. saccharum)",
    legend=:right, ylim=(0, 300))
plot!(p, years, D_vals, label="Dead biomass (D)", linewidth=2, linestyle=:dash)
p
```

### Species Comparison: Short-Lived vs. Long-Lived Species

This analysis compares the growth trajectories of two contrasting species: a
short-lived, shade-intolerant species (*P. banksiana*, max\_age = 70 years,
ANPP\_MAX = 5.77 Mg/ha/yr) and a long-lived, shade-tolerant species
(*A. saccharum*, max\_age = 400 years, ANPP\_MAX = 7.45 Mg/ha/yr). The
short-lived species shows a rapid rise and decline in biomass due to age-related
mortality, while the long-lived species achieves higher peak biomass over a
longer time horizon.

```@example landis
# Short-lived species: P. banksiana
prob_short = ODEProblem(compiled,
    Dict(compiled.B => 5.0 * Mg_ha_to_kg_m2, compiled.D_wood => 0.0,
         compiled.max_age => 70.0 * yr_to_s,
         compiled.ANPP_MAX => 5.77 * Mg_ha_to_kg_m2 / yr_to_s),
    (0.0, tspan_s))
sol_short = solve(prob_short)

# Long-lived species: A. saccharum (defaults)
prob_long = ODEProblem(compiled, [], (0.0, tspan_s))
sol_long = solve(prob_long)

B_short = [sol_short(yr * yr_to_s; idxs=compiled.B) / Mg_ha_to_kg_m2 for yr in years]
B_long = [sol_long(yr * yr_to_s; idxs=compiled.B) / Mg_ha_to_kg_m2 for yr in years]

p = plot(years, B_long, label="A. saccharum (400 yr)", linewidth=2,
    xlabel="Simulation year",
    ylabel="Living Biomass (Mg ha⁻¹)",
    title="Species Longevity Comparison",
    legend=:topright, ylim=(0, 250))
plot!(p, years, B_short, label="P. banksiana (70 yr)", linewidth=2, linestyle=:dash)
p
```

### Competition Effect

This analysis demonstrates the effect of inter-cohort competition on growth. When
other cohorts occupy a significant fraction of the site's maximum biomass capacity
(B\_MAX\_site = 500 Mg/ha), the focal cohort's potential biomass (B\_POT) is
reduced, leading to lower ANPP and lower peak biomass. Competition only has an
effect when B\_other > B\_MAX\_site - B\_MAX (i.e., when the site is sufficiently
crowded that the growing space available to the cohort is less than its species
maximum).

```@example landis
tspan_comp = 150.0 * yr_to_s

# No competition
prob_no_comp = ODEProblem(compiled,
    Dict(compiled.B => 5.0 * Mg_ha_to_kg_m2, compiled.D_wood => 0.0,
         compiled.B_other => 0.0),
    (0.0, tspan_comp))
sol_no_comp = solve(prob_no_comp)

# High competition: 400 Mg/ha of other species
prob_comp = ODEProblem(compiled,
    Dict(compiled.B => 5.0 * Mg_ha_to_kg_m2, compiled.D_wood => 0.0,
         compiled.B_other => 400.0 * Mg_ha_to_kg_m2),
    (0.0, tspan_comp))
sol_comp = solve(prob_comp)

years_comp = 1:150
B_no_comp = [sol_no_comp(yr * yr_to_s; idxs=compiled.B) / Mg_ha_to_kg_m2 for yr in years_comp]
B_comp = [sol_comp(yr * yr_to_s; idxs=compiled.B) / Mg_ha_to_kg_m2 for yr in years_comp]

p = plot(years_comp, B_no_comp, label="No competition", linewidth=2,
    xlabel="Simulation year",
    ylabel="Living Biomass (Mg ha⁻¹)",
    title="Effect of Competition on Growth",
    legend=:topright)
plot!(p, years_comp, B_comp, label="B_other = 400 Mg/ha", linewidth=2, linestyle=:dash)
p
```

## Limitations

The following aspects of the paper are not implemented in this module:

- **Multi-cohort and multi-species interactions**: This implementation models a single cohort. The paper describes landscape-scale simulations with multiple species-age cohorts competing within sites.
- **Disturbance modules**: Fire, wind, and harvesting disturbance processes (Sections 4.2.1, 4.2.2) are not included.
- **Shade calculation**: The percent full sunlight calculation (Eq. 8) and shade tolerance classes are not implemented.
- **Leaf litter partitioning**: The fine/coarse dead biomass partitioning described in Section 3.2.4 is not modeled.
- **PnET-II coupling**: The ANPP\_MAX parameterization via PnET-II (Section 4.3) is external to this module.
