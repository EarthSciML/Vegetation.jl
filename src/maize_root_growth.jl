"""
    MaizeRootGrowth(; name=:MaizeRootGrowth)

A diffusive root growth model for maize (*Zea mays* L.) that treats root carbon
expansion as a diffusion process in a 2D soil domain.

The model consists of three coupled processes:
1. **Carbon Allocation**: Partitions photosynthetic carbon between shoot and root growth
   based on plant development stage and water stress
2. **Carbon Assignment**: Distributes allocated root carbon to spatial locations based on
   four favorability indices (penetration resistance, temperature, aeration, root density)
3. **Root Diffusion**: Solves the spatial distribution of young and mature root density
   using a diffusion PDE for young roots and a local ODE for mature roots

The diffusion coefficient depends on soil water potential and temperature. Young roots
diffuse spatially and mature into non-diffusing mature roots at rate T_YM.

Units: The paper uses mixed units (bar, cm, Mg/m³). This implementation uses SI base
units (Pa for pressure, m for length, kg/m³ for density, s for time). Empirical
equations from Acock et al. (1985) use non-dimensionalized inputs obtained by dividing
by SI reference values with appropriate units, following the non-dimensionalization
technique for fractional powers.

**Reference**: Wang, Z., Timlin, D., Li, S., Fleisher, D., Dathe, A., Luo, C., Dong, L.,
Reddy, V.R., and Tully, K. (2021). A diffusive model of maize root growth in MAIZSIM
and its applications in Ridge-Furrow Rainfall Harvesting. *Agricultural Water Management*,
254, 106966. doi:10.1016/j.agwat.2021.106966

$(SIGNATURES)
"""
@component function MaizeRootGrowth(; name = :MaizeRootGrowth)
    # ===== Reference constants for non-dimensionalization =====
    # Paper uses: bar for pressure, Mg/m³ for bulk density, mol/L for concentration.
    # We divide SI inputs by these reference values to obtain dimensionless ratios
    # before applying the empirical formulas.
    @constants begin
        P_ref = 1.0e5, [description = "Reference pressure (1 bar)", unit = u"Pa"]
        ρ_ref = 1000.0, [description = "Reference density (1 Mg/m³)", unit = u"kg/m^3"]
        C_ref = 1000.0, [description = "Reference concentration (1 mol/L)", unit = u"mol/m^3"]
        zero_diffusivity = 0.0, [description = "Zero diffusivity", unit = u"m^2/s"]
    end

    # ===== Physical constants =====
    @constants begin
        T_freeze = 273.15, [description = "Freezing point of water", unit = u"K"]
    end

    # ===== Empirical constants from Eq. (1) =====
    @constants begin
        # Penetration resistance constants (f₁, Acock et al. 1985)
        # f₁ = (1/2)(ψ_trd[bar] - 5.4|ψ[bar]|^0.25·exp(-10.58(1.7[Mg/m³] - ρ_b[Mg/m³]))) - (1/4)(ψ_trd[bar] - ψ[bar])
        # After non-dimensionalization by P_ref and ρ_ref, the coefficients are dimensionless.
        pen_coeff = 5.4, [description = "Penetration resistance coefficient in f₁ (dimensionless)"]
        pen_exp = 0.25, [description = "Exponent on |ψ/P_ref| in f₁ (dimensionless)"]
        pen_bd_coeff = 10.58, [description = "Bulk density coefficient in f₁ exponent (dimensionless)"]
        ρ_bd_crit = 1700.0, [description = "Critical bulk density in f₁ (1.7 Mg/m³)", unit = u"kg/m^3"]

        # Temperature favorability constants (f₂)
        T_opt_low = 291.15, [description = "Lower optimal temperature for f₂ (18°C)", unit = u"K"]
        T_opt_high = 306.15, [description = "Upper optimal temperature for f₂ (33°C)", unit = u"K"]
        T_exp = 1.66, [description = "Temperature exponent in f₂ (dimensionless)"]
        T_guard = 0.01, [description = "Guard temperature difference to prevent 0^exp in f₂", unit = u"K"]

        # Aeration constants (f₃)
        O2_thresh = 20.0, [description = "O₂ concentration threshold in f₃ (0.02 mol/L)", unit = u"mol/m^3"]
        O2_exp = 7.14, [description = "Exponent in f₃ (dimensionless)"]

        # Root density threshold (f₄)
        root_dens_thresh = 0.03, [description = "Root density threshold in f₄", unit = u"kg/m^3"]
    end

    # ===== Constants from Eq. (4) — Diffusion coefficient =====
    @constants begin
        ψ_s_diff = -14709.975, [description = "Wet limit water potential ψ_s for diffusion (-150 cm head)", unit = u"Pa"]
        ψ_r_diff = -49033.25, [description = "Dry limit water potential ψ_r for diffusion (-500 cm head)", unit = u"Pa"]
        T0_diff = 295.0, [description = "Reference temperature T₀ for diffusion", unit = u"K"]
        p_diff = 10000.0, [description = "Temperature parameter p in f̃₂ (dimensionless)"]
        q_diff = 1.0, [description = "Temperature parameter q in f̃₂ (dimensionless)"]
        u_diff = 18000.0, [description = "Temperature parameter u in f̃₂", unit = u"K"]
    end

    # ===== Parameters =====
    @parameters begin
        A_growth = 0.55 / 86400.0, [description = "Potential relative growth rate (Eq. 2, 0.55/day)", unit = u"s^-1"]
        T_YM = 0.1 / 86400.0, [description = "Young-to-mature root maturation rate", unit = u"s^-1"]
        D0_xx = 50.0 * 1.0e-4 / 86400.0, [description = "Potential horizontal diffusivity D⁰_x (50 cm²/day)", unit = u"m^2/s"]
        D0_zz = 3.0 * 1.0e-4 / 86400.0, [description = "Potential vertical diffusivity D⁰_z (3 cm²/day)", unit = u"m^2/s"]
        R_total = 0.001 / 86400.0, [description = "Total carbon input rate for root growth", unit = u"kg/m^3/s"]
        ψ_rtd = 5.0e5, [description = "Root turgor pressure at dawn", unit = u"Pa"]
        ρ_b = 1380.0, [description = "Soil bulk density", unit = u"kg/m^3"]
        ψ_soil = -3.0e4, [description = "Soil water potential", unit = u"Pa"]
        T_soil = 298.0, [description = "Soil temperature", unit = u"K"]
        O2_soil = 1020.0, [description = "Soil O₂ concentration (1.02 mol/L = well-aerated)", unit = u"mol/m^3"]
    end

    # ===== State Variables =====
    @variables begin
        Y(t) = 0.001, [description = "Young root density", unit = u"kg/m^3"]
        M(t) = 0.0, [description = "Mature root density", unit = u"kg/m^3"]
    end

    # ===== Auxiliary Variables =====
    @variables begin
        f1(t), [description = "Penetration resistance favorability index (dimensionless)"]
        f2(t), [description = "Temperature favorability index (dimensionless)"]
        f3(t), [description = "Aeration favorability index (dimensionless)"]
        f4(t), [description = "Root density favorability index (dimensionless)"]
        f_min(t), [description = "Minimum favorability index (dimensionless)"]
        R_bar(t), [description = "Potential C for root growth (Eq. 2)", unit = u"kg/m^3/s"]
        D_eff_xx(t), [description = "Effective horizontal diffusion coefficient", unit = u"m^2/s"]
        D_eff_zz(t), [description = "Effective vertical diffusion coefficient", unit = u"m^2/s"]
        f_tilde_psi(t), [description = "Water potential factor for diffusion (dimensionless)"]
        f_tilde_T(t), [description = "Temperature factor for diffusion (dimensionless)"]
    end

    eqs = [
        # --- Eq. (1): Favorability indices ---

        # f₁ - Penetration resistance (Eq. 1, Acock et al. 1985)
        # f₁ = (1/2)(ψ_trd[bar] - 5.4|ψ[bar]|^0.25·exp(-10.58(1.7[Mg/m³] - ρ_b[Mg/m³]))) - (1/4)(ψ_trd[bar] - ψ[bar])
        # Non-dimensionalize pressures by P_ref (1 bar) and densities by ρ_ref (1 Mg/m³)
        f1 ~ max(
            0.0, min(
                1.0,
                0.5 * (ψ_rtd / P_ref - pen_coeff * abs(ψ_soil / P_ref)^pen_exp * exp(-pen_bd_coeff * (ρ_bd_crit - ρ_b) / ρ_ref))
                    - 0.25 * (ψ_rtd / P_ref - ψ_soil / P_ref)
            )
        ),

        # f₂ - Temperature favorability (Eq. 1)
        # Piecewise: (T_C/18)^1.66 for 0<T_C<18°C, 1 for 18≤T_C<33°C, (T_C/33)^(-1.66) for T_C≥33°C
        # Non-dimensionalize: T_C/18 = (T_soil - T_freeze)/(T_opt_low - T_freeze), giving K/K = dimensionless
        f2 ~ max(
            0.0, min(
                1.0, ifelse(
                    T_soil < T_opt_low,
                    (max(T_guard, T_soil - T_freeze) / (T_opt_low - T_freeze))^T_exp,
                    ifelse(
                        T_soil < T_opt_high,
                        1.0,
                        ((T_soil - T_freeze) / (T_opt_high - T_freeze))^(-T_exp)
                    )
                )
            )
        ),

        # f₃ - Aeration / O₂ favorability (Eq. 1)
        # f₃ = ([O₂] - 0.02)^7.14 where [O₂] is in mol/L
        # Non-dimensionalize: ([O₂] - 0.02)[mol/L] = (O2_soil - O2_thresh)/C_ref, giving (mol/m³)/(mol/m³) = dimensionless
        f3 ~ max(
            0.0, min(
                1.0,
                (max(0.0, (O2_soil - O2_thresh) / C_ref))^O2_exp
            )
        ),

        # f₄ - Root density favorability (Eq. 1)
        # f₄ = 1 - min(1, (M+Y)/0.03) where M+Y in kg/m³
        f4 ~ max(0.0, 1.0 - min(1.0, (M + Y) / root_dens_thresh)),

        # Minimum favorability
        f_min ~ min(f1, min(f2, min(f3, f4))),

        # --- Eq. (2): Potential C assigned for root growth ---
        # R̄ = (M + Y) × A × min{f₁, f₂, f₃, f₄}
        R_bar ~ (M + Y) * A_growth * f_min,

        # --- Eq. (4): Diffusion coefficient ---
        # D = D⁰ × min{f̃₁(ψ), f̃₂(T)}

        # f̃₁(ψ) — water potential factor for diffusion
        # f̃₁(ψ) = (1/2)sin(π(ψ - (ψ_s+ψ_r)/2) / (ψ_s - ψ_r)) + 1/2
        # Work directly in Pa: (Pa - Pa)/(Pa - Pa) = dimensionless argument to sin
        f_tilde_psi ~ max(
            0.0, min(
                1.0,
                0.5 * sin(π * (ψ_soil - (ψ_s_diff + ψ_r_diff) / 2.0) / (ψ_s_diff - ψ_r_diff)) + 0.5
            )
        ),

        # f̃₂(T) — temperature factor for diffusion
        # f̃₂(T) = max{[(1 + e^(p-T/T₀))·e^(T/T₀-p)] / (1 + e^(q-u/T)), 1.0}
        # Algebraic simplification: (1+e^a)·e^(-a) = e^(-a) + 1 = 1 + e^(-a)
        # This avoids overflow from exp(p - T/T₀) when p ≫ T/T₀.
        # T/T₀ = K/K = dimensionless; u/T = u_diff/T_soil = K/K = dimensionless
        f_tilde_T ~ max(
            1.0,
            (1.0 + exp(T_soil / T0_diff - p_diff))
                / (1.0 + exp(q_diff - u_diff / T_soil))
        ),

        # Effective diffusion coefficients (Eq. 4: D_* = D⁰_* × min{f̃₁, f̃₂})
        D_eff_xx ~ max(zero_diffusivity, D0_xx * min(f_tilde_psi, f_tilde_T)),
        D_eff_zz ~ max(zero_diffusivity, D0_zz * min(f_tilde_psi, f_tilde_T)),

        # --- Eq. (3a): Young root dynamics ---
        # ∂Y/∂t = ∇·[D∇Y] + R - T_YM·Y
        # In the ODE (non-spatial) form, the diffusion term is zero.
        # The spatial PDE form should be implemented using MethodOfLines.jl
        # for 2D simulations. Here we provide the local (non-spatial) dynamics.
        D(Y) ~ min(R_bar, R_total) - T_YM * Y,  # Eq. 3a (local, without diffusion)

        # --- Eq. (3b): Mature root dynamics ---
        D(M) ~ T_YM * Y,  # Eq. 3b
    ]

    return System(eqs, t; name)
end
