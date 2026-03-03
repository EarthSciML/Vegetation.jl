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
units (Pa for pressure, m for length, kg/m³ for density, s for time).

**Reference**: Wang, Z., Timlin, D., Li, S., Fleisher, D., Dathe, A., Luo, C., Dong, L.,
Reddy, V.R., and Tully, K. (2021). A diffusive model of maize root growth in MAIZSIM
and its applications in Ridge-Furrow Rainfall Harvesting. *Agricultural Water Management*,
254, 106966. doi:10.1016/j.agwat.2021.106966

$(SIGNATURES)
"""
@component function MaizeRootGrowth(; name = :MaizeRootGrowth)
    # ===== Unit conversion constants =====
    # Paper uses: bar for pressure, cm for water potential head, Mg/m³ for bulk density,
    # cm²/day for diffusivity, °C for temperature (f₂), K for temperature (f̃₂).
    # We convert everything to SI base units.

    @constants begin
        one_bar = 1.0e5, [description = "1 bar in Pa", unit = u"Pa"]
        one_day = 86400.0, [description = "1 day in seconds", unit = u"s"]
        one_cm = 0.01, [description = "1 cm in meters", unit = u"m"]
        one_cm_sq = 1.0e-4, [description = "1 cm² in m²", unit = u"m^2"]
        one_Mg_per_m3 = 1000.0, [description = "1 Mg/m³ in kg/m³", unit = u"kg/m^3"]
        one_mol_per_L = 1000.0, [description = "1 mol/L in mol/m³", unit = u"mol/m^3"]
        one_K = 1.0, [description = "Unit temperature", unit = u"K"]
        one_Pa = 1.0, [description = "Unit pressure", unit = u"Pa"]
        one_kg_per_m3 = 1.0, [description = "Unit density", unit = u"kg/m^3"]
        zero_kg_per_m3 = 0.0, [description = "Zero density", unit = u"kg/m^3"]
        one_m2_per_s = 1.0, [description = "Unit diffusivity", unit = u"m^2/s"]
        zero_m2_per_s = 0.0, [description = "Zero diffusivity", unit = u"m^2/s"]
    end

    # ===== Constants from Eq. (1) =====
    @constants begin
        # Penetration resistance constants (f₁)
        pen_coeff = 5.4, [description = "Penetration resistance coefficient in f₁ (dimensionless)"]
        pen_bd_coeff = 10.58, [description = "Bulk density coefficient in f₁ exponent (dimensionless)"]
        pen_bd_ref = 1.7, [description = "Reference bulk density in f₁ (Mg/m³ equivalent, dimensionless)"]
        pen_quarter = 0.25, [description = "Quarter coefficient in f₁ (dimensionless)"]

        # Root density threshold (f₄)
        root_dens_thresh = 0.03, [description = "Root density threshold in f₄", unit = u"kg/m^3"]

        # Aeration threshold (f₃)
        O2_thresh = 0.02, [description = "O₂ concentration threshold in f₃", unit = u"mol/m^3"]
    end

    # ===== Constants from Eq. (4) — Diffusion coefficient =====
    @constants begin
        psi_s = -150.0, [description = "Water potential parameter ψ_s (cm equivalent, dimensionless)"]
        psi_r = -500.0, [description = "Water potential parameter ψ_r (cm equivalent, dimensionless)"]
        T0_diff = 295.0, [description = "Reference temperature T₀ for diffusion", unit = u"K"]
        p_diff = 10000.0, [description = "Temperature parameter p (dimensionless)"]
        q_diff = 1.0, [description = "Temperature parameter q (dimensionless)"]
        u_diff = 18000.0, [description = "Temperature parameter u (dimensionless)"]
    end

    # ===== Parameters =====
    @parameters begin
        A_growth = 0.55 / 86400.0, [description = "Potential relative growth rate (Eq. 2, 0.55/day)", unit = u"s^-1"]
        T_YM = 0.1 / 86400.0, [description = "Young-to-mature root maturation rate", unit = u"s^-1"]
        D0_xx = 50.0 * 1.0e-4 / 86400.0, [description = "Potential vertical diffusivity (50 cm²/day)", unit = u"m^2/s"]
        D0_zz = 3.0 * 1.0e-4 / 86400.0, [description = "Potential horizontal diffusivity (3 cm²/day)", unit = u"m^2/s"]
        R_total = 0.001 / 86400.0, [description = "Total carbon input rate for root growth", unit = u"kg/m^3/s"]
        ψ_rtd = 5.0e5, [description = "Root turgor pressure", unit = u"Pa"]
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
        D_eff(t), [description = "Effective diffusion coefficient", unit = u"m^2/s"]
        f_tilde_psi(t), [description = "Water potential factor for diffusion (dimensionless)"]
        f_tilde_T(t), [description = "Temperature factor for diffusion (dimensionless)"]
        T_celsius(t), [description = "Soil temperature in Celsius (dimensionless ratio)"]
        psi_cm(t), [description = "Soil water potential in cm head (dimensionless ratio)"]
        psi_rtd_bar(t), [description = "Root turgor pressure in bar (dimensionless ratio)"]
        psi_soil_bar(t), [description = "Soil water potential in bar (dimensionless ratio)"]
        rho_b_Mg(t), [description = "Bulk density in Mg/m³ (dimensionless ratio)"]
        O2_mol_per_L(t), [description = "O₂ in mol/L (dimensionless ratio)"]
    end

    eqs = [
        # --- Unit conversions for empirical equations ---
        T_celsius ~ T_soil / one_K - 273.15,                    # K to °C (dimensionless)
        psi_cm ~ ψ_soil / (one_Pa * 98.0665),                   # Pa to cm water head (dimensionless)
        psi_rtd_bar ~ ψ_rtd / one_bar,                          # Pa to bar (dimensionless)
        psi_soil_bar ~ ψ_soil / one_bar,                        # Pa to bar (dimensionless)
        rho_b_Mg ~ ρ_b / one_Mg_per_m3,                         # kg/m³ to Mg/m³ (dimensionless)
        O2_mol_per_L ~ O2_soil / one_mol_per_L,                 # mol/m³ to mol/L (dimensionless)

        # --- Eq. (1): Favorability indices ---

        # f₁ - Penetration resistance (Eq. 1, Acock et al. 1985)
        # Uses bar units internally; result clamped to [0,1]
        f1 ~ max(
            0.0, min(
                1.0,
                (psi_rtd_bar - pen_coeff * abs(psi_soil_bar)^0.25 * exp(-pen_bd_coeff * (pen_bd_ref - rho_b_Mg)))
                    - pen_quarter * (psi_rtd_bar - psi_soil_bar)
            )
        ),

        # f₂ - Temperature favorability (Eq. 1)
        # Piecewise: (T/18)^1.66 for T<18°C, 1 for 18-33°C, (T/33)^(-1.66) for T≥33°C
        f2 ~ max(
            0.0, min(
                1.0, ifelse(
                    T_celsius < 18.0,
                    (max(0.01, T_celsius) / 18.0)^1.66,
                    ifelse(
                        T_celsius < 33.0,
                        1.0,
                        (T_celsius / 33.0)^(-1.66)
                    )
                )
            )
        ),

        # f₃ - Aeration / O₂ favorability (Eq. 1)
        f3 ~ max(
            0.0, min(
                1.0,
                (max(0.0, O2_mol_per_L - 0.02))^7.14
            )
        ),

        # f₄ - Root density favorability (Eq. 1)
        f4 ~ max(0.0, 1.0 - min(1.0, (M + Y) / root_dens_thresh)),

        # Minimum favorability
        f_min ~ min(f1, min(f2, min(f3, f4))),

        # --- Eq. (2): Potential C assigned for root growth ---
        R_bar ~ (M + Y) * A_growth * f_min,

        # --- Eq. (4): Diffusion coefficient ---
        # f̃₁(ψ) — water potential factor for diffusion
        # sin argument uses cm water head values
        f_tilde_psi ~ max(
            0.0, min(
                1.0,
                0.5 * sin(π / ((psi_s - psi_r) * (psi_cm - (psi_s + psi_r) / 2.0))) + 0.5
            )
        ),

        # f̃₂(T) — temperature factor for diffusion
        f_tilde_T ~ max(
            1.0,
            ((1.0 + exp(p_diff - p_diff / (T_soil / T0_diff))) * exp(p_diff / (T_soil / T0_diff) - p_diff))
                / (1.0 + exp(q_diff - u_diff / (T_soil / one_K)))
        ),

        # Effective isotropic diffusion coefficient (using vertical D0 as representative)
        D_eff ~ max(zero_m2_per_s, D0_xx * min(f_tilde_psi, one_m2_per_s / one_m2_per_s * f_tilde_T) * one_m2_per_s / one_m2_per_s),

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
