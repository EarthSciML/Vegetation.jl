"""
    ResidueMulchDecomposition(; name=:ResidueMulchDecomposition)

A residue mulch decomposition and nitrogen mineralization model based on a modified
CERES-N approach.

This component implements the mulch decomposition model from Wang et al. (2021),
which tracks three residue mass pools (carbohydrate CARB, cellulose CELL, lignin LIGN)
and their associated nitrogen mass pools. Decomposition rates depend on temperature,
water potential, and C/N ratio through adjustment factors (MTRF and CNRF). The model
also computes nitrogen mineralization, immobilization, and humification.

The model represents the "contacting portion" of residue mulch near the soil surface,
where decomposition actively occurs. A feeding rate parameter controls transfer from
the non-contacting portion to the contacting portion.

Units: The paper uses g, m², days, °C, and meters of water. This implementation uses
SI base units (kg for mass, m² for area, seconds for time, K for temperature, m for
water potential). Conversion: 1 day = 86400 s.

**Reference**: Wang, Z., Thapa, R., Timlin, D., Li, S., Sun, W., Beegum, S., et al.
(2021). Simulations of water and thermal dynamics for soil surfaces with residue mulch
and surface runoff. *Water Resources Research*, 57, e2021WR030431.
https://doi.org/10.1029/2021WR030431

$(SIGNATURES)
"""
@component function ResidueMulchDecomposition(; name = :ResidueMulchDecomposition)

    @constants begin
        one_day = 86400.0, [description = "One day in seconds", unit = u"s"]
        one_kg_m2 = 1.0, [description = "Unit mass per area for nondimensionalization", unit = u"kg/m^2"]
        zero_rate = 0.0, [description = "Zero mass flux rate", unit = u"kg/m^2/s"]
        T_freeze = 273.15, [description = "Freezing point of water", unit = u"K"]
        C_frac = 0.41, [description = "Mass fraction of C in residue (dimensionless)"]
    end

    @parameters begin
        # Decay rate coefficients (Eq. 21, Table 2)
        # Paper values in d⁻¹, converted to s⁻¹
        D_CARB = 0.43 / 86400.0, [description = "CARB base decay coefficient (0.43 d⁻¹)", unit = u"s^-1"]
        D_CELL = 0.24 / 86400.0, [description = "CELL base decay coefficient (0.24 d⁻¹)", unit = u"s^-1"]
        D_LIGN = 0.0228 / 86400.0, [description = "LIGN decay coefficient (0.0228 d⁻¹)", unit = u"s^-1"]
        γ_lign = 12.0, [description = "Lignin inhibition parameter (dimensionless)"]

        # MTRF coefficients (Eq. 22, Table 2)
        # T in °C and h in meters of water in original; we convert internally
        a_MTRF = 0.384, [description = "MTRF temperature intercept (dimensionless)"]
        b_MTRF = 0.0187, [description = "MTRF temperature slope", unit = u"K^-1"]
        c_MTRF = 0.142, [description = "MTRF water potential intercept", unit = u"m^-1"]
        d_MTRF = 0.628, [description = "MTRF water potential-temperature interaction", unit = u"m^-1*K"]

        # CNRF coefficients (Eq. 22, Table 2)
        a_CNRF = 0.693, [description = "CNRF exponential coefficient (dimensionless)"]
        CNR_crit = 13.0, [description = "Critical C/N ratio (dimensionless)"]

        # N dynamics (Eq. 23)
        HUMF = 0.125, [description = "Humification factor (dimensionless)"]
        # 0.0213 = C_frac / target_CN for immobilization (dimensionless)
        immob_coeff = 0.0213, [description = "Immobilization coefficient (dimensionless)"]

        # Feeding rate (Eq. 24)
        α_F = 0.1, [description = "Feeding rate from non-contacting to contacting portion (dimensionless)"]

        # Initial mass fractions (Table 2)
        F_CARB_ini = 0.2, [description = "Initial CARB mass fraction (dimensionless)"]
        F_CELL_ini = 0.7, [description = "Initial CELL mass fraction (dimensionless)"]
        F_LIGN_ini = 0.1, [description = "Initial LIGN mass fraction (dimensionless)"]

        # Initial N mass fractions (Table 2)
        FN_CARB_ini = 0.08, [description = "Initial N fraction in CARB (dimensionless)"]
        FN_CELL_ini = 0.01, [description = "Initial N fraction in CELL (dimensionless)"]
        FN_LIGN_ini = 0.01, [description = "Initial N fraction in LIGN (dimensionless)"]

        # Initial total residue mass per unit area (Table 2: 1200 g/m² = 0.12 kg/m²)
        RM_T_init = 0.12, [description = "Initial total residue mass (1200 g/m²)", unit = u"kg/m^2"]

        # Forcing inputs (external conditions)
        T_mulch = 293.15, [description = "Mulch temperature (forcing)", unit = u"K"]
        h_mulch = -1.0, [description = "Mulch matric potential (forcing)", unit = u"m"]
    end

    @variables begin
        # Contacting portion residue mass pools (Eq. 18, 25)
        RM_CARB(t) = 0.12 * 0.2, [description = "CARB residue mass in contacting portion", unit = u"kg/m^2"]
        RM_CELL(t) = 0.12 * 0.7, [description = "CELL residue mass in contacting portion", unit = u"kg/m^2"]
        RM_LIGN(t) = 0.12 * 0.1, [description = "LIGN residue mass in contacting portion", unit = u"kg/m^2"]

        # Non-contacting portion residue mass pools (Eq. 24)
        RM_CARB_nc(t) = 0.0, [description = "CARB residue mass in non-contacting portion", unit = u"kg/m^2"]
        RM_CELL_nc(t) = 0.0, [description = "CELL residue mass in non-contacting portion", unit = u"kg/m^2"]
        RM_LIGN_nc(t) = 0.0, [description = "LIGN residue mass in non-contacting portion", unit = u"kg/m^2"]

        # Contacting portion N mass pools (Eq. 19, 25)
        RMN_CARB(t) = 0.12 * 0.2 * 0.08, [description = "N mass in CARB contacting portion", unit = u"kg/m^2"]
        RMN_CELL(t) = 0.12 * 0.7 * 0.01, [description = "N mass in CELL contacting portion", unit = u"kg/m^2"]
        RMN_LIGN(t) = 0.12 * 0.1 * 0.01, [description = "N mass in LIGN contacting portion", unit = u"kg/m^2"]

        # Non-contacting portion N mass pools (Eq. 24)
        RMN_CARB_nc(t) = 0.0, [description = "N mass in CARB non-contacting portion", unit = u"kg/m^2"]
        RMN_CELL_nc(t) = 0.0, [description = "N mass in CELL non-contacting portion", unit = u"kg/m^2"]
        RMN_LIGN_nc(t) = 0.0, [description = "N mass in LIGN non-contacting portion", unit = u"kg/m^2"]

        # Mineral N available (forcing/state)
        N_inorg(t) = 0.001, [description = "Mineral (microbe-available) N", unit = u"kg/m^2"]

        # Diagnostic variables
        RM_T(t), [description = "Total residue mass (contacting)", unit = u"kg/m^2"]
        RMN_T(t), [description = "Total N mass (contacting)", unit = u"kg/m^2"]
        f_LIGN(t), [description = "LIGN mass fraction (dimensionless)"]
        CNR(t), [description = "C/N ratio (dimensionless)"]

        # Decay coefficients (Eq. 21)
        d_CARB(t), [description = "CARB effective decay rate", unit = u"s^-1"]
        d_CELL(t), [description = "CELL effective decay rate", unit = u"s^-1"]
        d_LIGN(t), [description = "LIGN effective decay rate", unit = u"s^-1"]

        # Adjustment factors (Eq. 22)
        MTRF(t), [description = "Moisture-temperature reduction factor (dimensionless)"]
        CNRF(t), [description = "C/N ratio reduction factor (dimensionless)"]

        # Decomposition rates (Eq. 20)
        R_dcmp_CARB(t), [description = "CARB decomposition rate", unit = u"kg/m^2/s"]
        R_dcmp_CELL(t), [description = "CELL decomposition rate", unit = u"kg/m^2/s"]
        R_dcmp_LIGN(t), [description = "LIGN decomposition rate", unit = u"kg/m^2/s"]
        R_dcmp_total(t), [description = "Total decomposition rate", unit = u"kg/m^2/s"]

        # N mineralization rates (Eq. 20)
        RN_mine_CARB(t), [description = "CARB gross N mineralization rate", unit = u"kg/m^2/s"]
        RN_mine_CELL(t), [description = "CELL gross N mineralization rate", unit = u"kg/m^2/s"]
        RN_mine_LIGN(t), [description = "LIGN gross N mineralization rate", unit = u"kg/m^2/s"]
        RN_mine_total(t), [description = "Total gross N mineralization rate", unit = u"kg/m^2/s"]

        # N dynamics (Eq. 23)
        N_Im(t), [description = "N immobilization rate", unit = u"kg/m^2/s"]
        RN_Net(t), [description = "Net N mineralization rate (mulch to soil)", unit = u"kg/m^2/s"]
        N_Humi(t), [description = "N humification rate", unit = u"kg/m^2/s"]

        # Feeding rates from non-contacting to contacting (Eq. 24)
        feed_RM_CARB(t), [description = "CARB mass feeding rate", unit = u"kg/m^2/s"]
        feed_RM_CELL(t), [description = "CELL mass feeding rate", unit = u"kg/m^2/s"]
        feed_RM_LIGN(t), [description = "LIGN mass feeding rate", unit = u"kg/m^2/s"]
        feed_RMN_CARB(t), [description = "CARB N feeding rate", unit = u"kg/m^2/s"]
        feed_RMN_CELL(t), [description = "CELL N feeding rate", unit = u"kg/m^2/s"]
        feed_RMN_LIGN(t), [description = "LIGN N feeding rate", unit = u"kg/m^2/s"]
    end

    # Use T_C for temperature in Celsius (for MTRF calculation)
    # Use h_m for water potential in meters (already in SI)
    eqs = [
        # Eq. 18 - Total residue mass in contacting portion
        RM_T ~ RM_CARB + RM_CELL + RM_LIGN,

        # Eq. 19 - Total N mass in contacting portion
        RMN_T ~ RMN_CARB + RMN_CELL + RMN_LIGN,

        # LIGN mass fraction (used in Eq. 21, 26)
        f_LIGN ~ RM_LIGN / max(RM_T, 1.0e-10 * one_kg_m2),

        # C/N ratio (Eq. 22): CNR = 0.41 * RM_T / (RMN_T + N_inorg)
        CNR ~ C_frac * RM_T / max(RMN_T + N_inorg, 1.0e-10 * one_kg_m2),

        # Eq. 21 - Decay coefficients
        d_CARB ~ D_CARB * exp(-γ_lign * f_LIGN),
        d_CELL ~ D_CELL * exp(-γ_lign * f_LIGN),
        d_LIGN ~ D_LIGN,

        # Eq. 22 - MTRF: moisture-temperature reduction factor
        # T in °C for the empirical formula, h in meters
        MTRF ~ ifelse(
            T_mulch > T_freeze,
            (a_MTRF + b_MTRF * (T_mulch - T_freeze)) *
                exp((c_MTRF + d_MTRF / (T_mulch - T_freeze)) * h_mulch),
            0.0
        ),

        # Eq. 22 - CNRF: C/N ratio reduction factor
        CNRF ~ ifelse(
            CNR > CNR_crit,
            exp(-a_CNRF * (CNR - CNR_crit) / CNR_crit),
            1.0
        ),

        # Eq. 20 - Decomposition rates
        R_dcmp_CARB ~ d_CARB * RM_CARB * max(MTRF * CNRF, 0.0),
        R_dcmp_CELL ~ d_CELL * RM_CELL * max(MTRF * CNRF, 0.0),
        R_dcmp_LIGN ~ d_LIGN * RM_LIGN * max(MTRF * CNRF, 0.0),
        R_dcmp_total ~ R_dcmp_CARB + R_dcmp_CELL + R_dcmp_LIGN,

        # Eq. 20 - Gross N mineralization rates
        RN_mine_CARB ~ d_CARB * RMN_CARB * max(MTRF * CNRF, 0.0),
        RN_mine_CELL ~ d_CELL * RMN_CELL * max(MTRF * CNRF, 0.0),
        RN_mine_LIGN ~ d_LIGN * RMN_LIGN * max(MTRF * CNRF, 0.0),
        RN_mine_total ~ RN_mine_CARB + RN_mine_CELL + RN_mine_LIGN,

        # Eq. 23 - N immobilization
        # N_Im = max(min(0.0213 * Σ R_dcmp - Σ RN_mine, N_inorg/Δt), 0)
        # In continuous form: N_Im is a rate (kg/m²/s)
        N_Im ~ max(min(immob_coeff * R_dcmp_total - RN_mine_total, N_inorg / one_day), zero_rate),

        # Eq. 23 - Net N mineralization (mulch to soil)
        RN_Net ~ (1 - HUMF) * RN_mine_total - N_Im,

        # Eq. 23 - N humification
        N_Humi ~ HUMF * RN_mine_total,

        # Feeding rates from non-contacting to contacting (Eq. 24)
        # Rate = α_F * total decomposition rate, distributed by initial fractions
        feed_RM_CARB ~ α_F * F_CARB_ini * R_dcmp_total,
        feed_RM_CELL ~ α_F * F_CELL_ini * R_dcmp_total,
        feed_RM_LIGN ~ α_F * F_LIGN_ini * R_dcmp_total,
        feed_RMN_CARB ~ α_F * FN_CARB_ini * RN_mine_total,
        feed_RMN_CELL ~ α_F * FN_CELL_ini * RN_mine_total,
        feed_RMN_LIGN ~ α_F * FN_LIGN_ini * RN_mine_total,

        # Eq. 25 - Contacting portion mass dynamics
        # dRM/dt = -R_dcmp + feeding_in + (for CARB: N_Im contribution via mass)
        D(RM_CARB) ~ -R_dcmp_CARB + feed_RM_CARB,
        D(RM_CELL) ~ -R_dcmp_CELL + feed_RM_CELL,
        D(RM_LIGN) ~ -R_dcmp_LIGN + feed_RM_LIGN,

        # Eq. 24 - Non-contacting portion loses mass via feeding
        D(RM_CARB_nc) ~ -feed_RM_CARB,
        D(RM_CELL_nc) ~ -feed_RM_CELL,
        D(RM_LIGN_nc) ~ -feed_RM_LIGN,

        # Eq. 25 - Contacting portion N dynamics
        # CARB receives immobilized N
        D(RMN_CARB) ~ -RN_mine_CARB + feed_RMN_CARB + N_Im,
        D(RMN_CELL) ~ -RN_mine_CELL + feed_RMN_CELL,
        D(RMN_LIGN) ~ -RN_mine_LIGN + feed_RMN_LIGN,

        # Eq. 24 - Non-contacting portion N dynamics
        D(RMN_CARB_nc) ~ -feed_RMN_CARB,
        D(RMN_CELL_nc) ~ -feed_RMN_CELL,
        D(RMN_LIGN_nc) ~ -feed_RMN_LIGN,

        # Mineral N dynamics: gains from net mineralization, consumed by immobilization
        D(N_inorg) ~ RN_Net,
    ]

    return System(eqs, t; name)
end

"""
    MulchRadiationAttenuation(; name=:MulchRadiationAttenuation, K=5)

Radiation attenuation through residue mulch layers.

This component implements the shortwave and longwave radiation attenuation model
from Wang et al. (2021), Equations 11-15. Shortwave radiation follows a geometric
attenuation pattern (Ross, 1976), while longwave radiation follows the
Stefan-Boltzmann law with layer-by-layer accumulation.

The model computes net radiation received by each mulch elemental layer, given
the incoming solar radiation, atmospheric conditions, and layer temperatures.

**Reference**: Wang, Z., Thapa, R., Timlin, D., et al. (2021). *Water Resources Research*,
57, e2021WR030431.

$(SIGNATURES)
"""
@component function MulchRadiationAttenuation(; name = :MulchRadiationAttenuation, K = 5)

    @constants begin
        σ_SB = 5.67e-8, [description = "Stefan-Boltzmann constant", unit = u"W/m^2/K^4"]
        one_W_m2 = 1.0, [description = "Unit radiation flux", unit = u"W/m^2"]
    end

    @parameters begin
        ΔR = 0.3, [description = "Residue-area index (dimensionless)"]
        Ω_cl = 0.6, [description = "Clumping index (dimensionless)"]
        α_m = 0.3, [description = "Mulch shortwave reflectivity (dimensionless)"]
        α_s = 0.15, [description = "Soil shortwave reflectivity (dimensionless)"]
        ε_m = 1.0, [description = "Mulch longwave emissivity (dimensionless)"]
        ε_s = 1.0, [description = "Soil longwave emissivity (dimensionless)"]
        ε_a = 0.8, [description = "Atmospheric longwave emissivity (dimensionless)"]
    end

    @variables begin
        # Forcing inputs
        S_0(t), [description = "Solar irradiance above mulch", unit = u"W/m^2"]
        T_a(t), [description = "Air temperature above mulch", unit = u"K"]
        T_s(t), [description = "Soil surface temperature", unit = u"K"]
        T_layer(t)[1:K], [description = "Temperature of each mulch layer", unit = u"K"]

        # Shortwave radiation at interfaces (Eq. 11, 12)
        S_d(t)[1:(K + 1)], [description = "Downward shortwave at each interface", unit = u"W/m^2"]
        S_u(t)[1:(K + 1)], [description = "Upward shortwave at each interface", unit = u"W/m^2"]

        # Longwave radiation at interfaces (Eq. 13, 14)
        L_d(t)[1:(K + 1)], [description = "Downward longwave at each interface", unit = u"W/m^2"]
        L_u(t)[1:(K + 1)], [description = "Upward longwave at each interface", unit = u"W/m^2"]

        # Net radiation (Eq. 15)
        R_net_iface(t)[1:(K + 1)], [description = "Net radiation at each interface", unit = u"W/m^2"]
        R_elem(t)[1:K], [description = "Net radiation received by each element", unit = u"W/m^2"]
    end

    # Transmissivity through n consecutive layers: τ_n = (1-ΔR)(1-Ω_cl*ΔR)^(n-1)
    # We precompute these symbolically
    τ = [(1 - ΔR) * (1 - Ω_cl * ΔR)^(n - 1) for n in 1:K]

    eqs = Equation[]

    # --- Eq. 11: Downward shortwave radiation ---
    # S_d[K+1] = S_0 (above mulch-air interface)
    push!(eqs, S_d[K + 1] ~ S_0)
    # S_d[K] = S_0 * (1 - ΔR)
    push!(eqs, S_d[K] ~ S_0 * (1 - ΔR))
    # S_d[k] = S_0 * (1-ΔR) * (1-Ω_cl*ΔR)^(K-k) for k = 1,...,K-1
    for k in 1:(K - 1)
        push!(eqs, S_d[k] ~ S_0 * τ[K - k])
    end

    # --- Eq. 12: Upward shortwave reflection ---
    # S_u[1] = S_d[1] * α_s  (reflection from soil)
    push!(eqs, S_u[1] ~ S_d[1] * α_s)
    if K >= 2
        # S_u[2] = S_d[1]*α_s*(1-ΔR) + α_m*(S_d[2] - S_d[1])
        push!(eqs, S_u[2] ~ S_d[1] * α_s * (1 - ΔR) + α_m * (S_d[2] - S_d[1]))
    end
    for k in 3:(K + 1)
        # S_u[k] = S_d[1]*α_s*τ_{k-1} + α_m*Σ_{j=2}^{k-1}(S_d[j]-S_d[j-1])*τ_{k-j} + α_m*(S_d[k]-S_d[k-1])
        reflected_soil = S_d[1] * α_s * τ[k - 1]
        reflected_mulch_layers = sum(α_m * (S_d[j] - S_d[j - 1]) * τ[k - j] for j in 2:(k - 1); init = 0.0 * one_W_m2)
        reflected_top = α_m * (S_d[k] - S_d[k - 1])
        push!(eqs, S_u[k] ~ reflected_soil + reflected_mulch_layers + reflected_top)
    end

    # --- Eq. 13: Downward longwave radiation ---
    # L_d[K+1] = ε_a * σ * T_a^4  (from atmosphere)
    push!(eqs, L_d[K + 1] ~ ε_a * σ_SB * T_a^4)
    if K >= 1
        # L_d[K] = ε_a*σ*T_a^4*τ_1 + ε_m*σ*T_layer[K]^4*(1-τ_1)
        push!(eqs, L_d[K] ~ ε_a * σ_SB * T_a^4 * τ[1] + ε_m * σ_SB * T_layer[K]^4 * (1 - τ[1]))
    end
    for k in 1:(K - 1)
        # L_d[k] = ε_a*σ*T_a^4*τ_{K-k+1} + ε_m*Σ_{j=k+1}^{K}σ*T_layer[j]^4*(τ_{j-k}-τ_{j-k+1}) + ε_m*σ*T_layer[k]^4*(1-τ_1)
        atmos_term = ε_a * σ_SB * T_a^4 * τ[K - k + 1]
        # Note: the sum goes from j=k+1 to K, and the paper says τ_{j-k} - τ_{j-k+1}
        mulch_sum = sum(ε_m * σ_SB * T_layer[j]^4 * (τ[j - k] - τ[j - k + 1]) for j in (k + 1):K; init = 0.0 * one_W_m2)
        adjacent_term = ε_m * σ_SB * T_layer[k]^4 * (1 - τ[1])
        push!(eqs, L_d[k] ~ atmos_term + mulch_sum + adjacent_term)
    end

    # --- Eq. 14: Upward longwave radiation ---
    # L_u[1] = ε_s * σ * T_s^4  (from soil)
    push!(eqs, L_u[1] ~ ε_s * σ_SB * T_s^4)
    if K >= 1
        # L_u[2] = ε_s*σ*T_s^4*τ_1 + ε_m*σ*T_layer[1]^4*(1-τ_1)
        push!(eqs, L_u[2] ~ ε_s * σ_SB * T_s^4 * τ[1] + ε_m * σ_SB * T_layer[1]^4 * (1 - τ[1]))
    end
    for k in 3:(K + 1)
        # L_u[k] = ε_s*σ*T_s^4*τ_k + ε_m*Σ_{j=1}^{k-2}σ*T_layer[j]^4*(τ_{k-j-1}-τ_{k-j}) + ε_m*σ*T_layer[k-1]^4*(1-τ_1)
        soil_term = ε_s * σ_SB * T_s^4 * τ[k - 1]
        mulch_sum = sum(ε_m * σ_SB * T_layer[j]^4 * (τ[k - j - 1] - τ[k - j]) for j in 1:(k - 2); init = 0.0 * one_W_m2)
        adjacent_term = ε_m * σ_SB * T_layer[k - 1]^4 * (1 - τ[1])
        push!(eqs, L_u[k] ~ soil_term + mulch_sum + adjacent_term)
    end

    # --- Eq. 15: Net radiation ---
    for k in 1:(K + 1)
        push!(eqs, R_net_iface[k] ~ S_d[k] - S_u[k] + L_d[k] - L_u[k])
    end
    for k in 1:K
        push!(eqs, R_elem[k] ~ R_net_iface[k] - R_net_iface[k + 1])
    end

    return System(eqs, t; name)
end

"""
    MulchWindProfile(; name=:MulchWindProfile, K=5)

Wind speed distribution within residue mulch.

Implements the vertical wind speed profile within residue mulch from Wang et al. (2021),
Equation 3. Mean wind speed attenuates exponentially from the mulch-air interface to
the mulch-soil interface, based on Novak et al. (2000a).

**Reference**: Wang, Z., Thapa, R., Timlin, D., et al. (2021). *Water Resources Research*,
57, e2021WR030431.

$(SIGNATURES)
"""
@component function MulchWindProfile(; name = :MulchWindProfile, K = 5)

    @constants begin
        k_K = 0.4, [description = "Von Karman's constant (dimensionless)"]
    end

    @parameters begin
        Z_M = 0.06, [description = "Mulch thickness", unit = u"m"]
        z_ref = 2.0, [description = "Reference height for wind measurement", unit = u"m"]
    end

    @variables begin
        u_ref(t), [description = "Wind speed at reference height", unit = u"m/s"]
        u_star(t), [description = "Friction velocity", unit = u"m/s"]
        u_layer(t)[1:K], [description = "Wind speed in each mulch layer", unit = u"m/s"]
    end

    # Surface displacement and roughness length
    d_disp = 0.87 * Z_M
    z_r = 0.079 * Z_M

    eqs = Equation[]

    # Eq. 3a - Friction velocity
    push!(eqs, u_star ~ k_K * u_ref / log((z_ref - d_disp) / z_r))

    # Eq. 3b - Wind speed in each layer
    # Layer k has center at z = (z_k + z_{k+1})/2
    # With uniform layers: z_k = (k-1)/K * Z_M, z_{k+1} = k/K * Z_M
    # Center = (2k-1)/(2K) * Z_M
    for k in 1:K
        z_center = (2k - 1) / (2K) * Z_M
        push!(eqs, u_layer[k] ~ 0.21 * u_star * exp(2.2 / Z_M * z_center))
    end

    return System(eqs, t; name)
end

"""
    MulchHeatVaporFluxes(; name=:MulchHeatVaporFluxes)

Heat and vapor flux parameterizations for residue mulch.

Implements the diffusive (Eq. 6), free convective (Eq. 7), and forced convective (Eq. 8)
flux formulations from Wang et al. (2021), along with the Rayleigh and Richardson number
diagnostics (Eq. 4) and flux regime selection (Eq. 5).

This component computes fluxes for the whole mulch layer based on the temperature
and wind speed differences between the mulch-air and mulch-soil interfaces.

**Reference**: Wang, Z., Thapa, R., Timlin, D., et al. (2021). *Water Resources Research*,
57, e2021WR030431.

$(SIGNATURES)
"""
@component function MulchHeatVaporFluxes(; name = :MulchHeatVaporFluxes)

    @constants begin
        ν_air = 1.5e-5, [description = "Kinematic viscosity of air", unit = u"m^2/s"]
        D_hm = 2.2e-5, [description = "Molecular thermal diffusivity of air", unit = u"m^2/s"]
        g_acc = 9.81, [description = "Gravitational acceleration", unit = u"m/s^2"]
        # α_c_fr has units m K^(-1/2) s^(-1) in the paper (5.6e-3).
        # Since fractional power units aren't supported, we nondimensionalize:
        # α_c_fr * sqrt(|ΔT|) = α_c_fr_nd * sqrt(|ΔT/one_K|), where α_c_fr_nd is in m/s
        α_c_fr_nd = 5.6e-3, [description = "Free convection coefficient (nondimensionalized, m/s)", unit = u"m/s"]
        M_H2O = 0.01802, [description = "Molecular weight of water", unit = u"kg/mol"]
        R_gas = 8.314, [description = "J/mol/K gas constant", unit = u"J/mol/K"]
        C_as = 718.0, [description = "Specific heat of dry air", unit = u"J/kg/K"]
        ρ_a = 1.204, [description = "Dry air density at 293.75K", unit = u"kg/m^3"]
        k_K = 0.4, [description = "Von Karman's constant (dimensionless)"]
        Ra_crit = 1706.0, [description = "Critical Rayleigh number (dimensionless)"]
        Ri_crit = 1.0, [description = "Critical Richardson number (dimensionless)"]
        one_K = 1.0, [description = "Unit temperature", unit = u"K"]
        one_ms = 1.0, [description = "Unit velocity", unit = u"m/s"]
    end

    @parameters begin
        Z_M = 0.06, [description = "Mulch thickness", unit = u"m"]
    end

    @variables begin
        # Interface temperatures and wind speeds
        T_top(t), [description = "Temperature at mulch-air interface", unit = u"K"]
        T_bot(t), [description = "Temperature at mulch-soil interface", unit = u"K"]
        u_top(t), [description = "Wind speed at top layer", unit = u"m/s"]
        u_bot(t), [description = "Wind speed at bottom layer", unit = u"m/s"]

        # Diagnostic numbers (Eq. 4)
        Ra(t), [description = "Rayleigh number (dimensionless)"]
        Ri(t), [description = "Richardson number (dimensionless)"]

        # Convective conductances (Eq. 7, 8)
        C_h_conv(t), [description = "Heat convective conductance", unit = u"W/m^2/K"]
        C_w_conv(t), [description = "Water vapor convective conductance", unit = u"kg/m^2/s/Pa"]

        # Average temperature
        T_avg(t), [description = "Average mulch temperature", unit = u"K"]
        ΔT(t), [description = "Temperature difference across mulch", unit = u"K"]
        Δu(t), [description = "Wind speed difference across mulch", unit = u"m/s"]
    end

    eqs = [
        # Average and differences
        T_avg ~ (T_top + T_bot) / 2,
        ΔT ~ T_top - T_bot,
        Δu ~ u_top - u_bot,

        # Eq. 4 - Rayleigh number: Ra = 2g|ΔT|Z_M³ / (T̄ · ν · D_hm)
        # ΔT(K) * Z_M³(m³) / (T_avg(K) * ν(m²/s) * D_hm(m²/s)) → m³/(m⁴/s²) * K/K → s²/m → need g(m/s²)
        # Full: m/s² * K * m³ / (K * m²/s * m²/s) = m⁴/(s² * m⁴/s²) = dimensionless ✓
        Ra ~ 2 * g_acc * abs(ΔT) * (Z_M^3) / (T_avg * ν_air * D_hm),

        # Eq. 4 - Richardson number: Ri = 2g|ΔT|Z_M / (T̄ · Δu²)
        # m/s² * K * m / (K * m²/s²) = m²/s² / (m²/s²) = dimensionless ✓
        Ri ~ 2 * g_acc * abs(ΔT) * Z_M / (T_avg * max(abs(Δu), 1.0e-6 * one_ms)^2),

        # Flux regime selection (Eq. 5) and convective conductances
        # Combined: if Ra < 1706 → no convection; Ra ≥ 1706 and Ri ≥ 1 → free; Ri < 1 → forced
        # Note: α_c_fr_nd * sqrt(|ΔT/one_K|) = original α_c_fr * sqrt(|ΔT|) with nondimensionalized units
        # Flux regime selection (Eq. 5) and convective conductances
        # For the zero branches, multiply by a small value instead of literal 0 to maintain units
        C_h_conv ~ ifelse(
            Ra < Ra_crit,
            1.0e-30 * α_c_fr_nd * C_as * ρ_a,
            ifelse(
                Ri >= Ri_crit,
                # Eq. 7 - Free convection: C_h^fr = α_c^fr * sqrt(|ΔT|) * C_as * ρ_a
                α_c_fr_nd * sqrt(max(abs(ΔT / one_K), 1.0e-10)) * C_as * ρ_a,
                # Eq. 8 - Forced convection: C_h^fo = 0.155 * k_K² * u * C_as * ρ_a
                0.155 * k_K^2 * abs(Δu) * C_as * ρ_a
            )
        ),

        C_w_conv ~ ifelse(
            Ra < Ra_crit,
            1.0e-30 * α_c_fr_nd * M_H2O / (R_gas * one_K),
            ifelse(
                Ri >= Ri_crit,
                # Eq. 7 - Free convection: C_w^fr = α_c^fr * sqrt(|ΔT|) * M_H2O / (R_g * T̄)
                α_c_fr_nd * sqrt(max(abs(ΔT / one_K), 1.0e-10)) * M_H2O / (R_gas * T_avg),
                # Eq. 8 - Forced convection: C_w^fo = 0.155 * k_K² * u * M_H2O / (R_g * T̄)
                0.155 * k_K^2 * abs(Δu) * M_H2O / (R_gas * T_avg)
            )
        ),
    ]

    return System(eqs, t; name)
end

"""
    MulchWaterCharacteristic(; name=:MulchWaterCharacteristic)

Water characteristic function for residue mulch.

Implements Equation 26 from Wang et al. (2021), which relates matric potential to
gravimetric water content in residue mulch. The water characteristic function depends
on the lignin mass fraction, reflecting changes in mulch structure during decomposition.

**Reference**: Wang, Z., Thapa, R., Timlin, D., et al. (2021). *Water Resources Research*,
57, e2021WR030431.

$(SIGNATURES)
"""
@component function MulchWaterCharacteristic(; name = :MulchWaterCharacteristic)

    @constants begin
        ρ_n = 1.0e6, [description = "Liquid water density", unit = u"g/m^3"]
        one_Pa = 1.0, [description = "Unit pressure for nondimensionalization", unit = u"Pa"]
        g_per_m2 = 1.0, [description = "Unit mass per area", unit = u"g/m^2"]
    end

    @parameters begin
        ρ_m = 20000.0, [description = "Mulch residue density", unit = u"g/m^3"]
        # Water characteristic function fitting parameters (Table 2, Eq. 26)
        # Original a_wrc1 = 20.1 MPa = 20.1e6 Pa
        a_wrc1 = 20.1e6, [description = "WRC parameter a1 (20.1 MPa)", unit = u"Pa"]
        a_wrc2 = 0.249, [description = "WRC parameter a2 (dimensionless)"]
        b_wrc1 = 0.324, [description = "WRC parameter b1 (dimensionless)"]
        b_wrc2 = 0.124, [description = "WRC parameter b2 (dimensionless)"]
        # Saturated water content fitting parameters
        a_sat1 = 7.1, [description = "Saturated gravimetric WC parameter a1 (dimensionless)"]
        a_sat2 = 0.079, [description = "Saturated gravimetric WC parameter a2 (dimensionless)"]
    end

    @variables begin
        θ_grav(t), [description = "Gravimetric water content (dimensionless)"]
        θ_grav_sat(t), [description = "Saturated gravimetric water content (dimensionless)"]
        f_LIGN(t), [description = "LIGN mass fraction (dimensionless)"]
        a_m_wrc(t), [description = "WRC parameter a_m", unit = u"Pa"]
        b_m_wrc(t), [description = "WRC exponent b_m (dimensionless)"]
        h_Pa(t), [description = "Matric potential", unit = u"Pa"]
    end

    eqs = [
        # Eq. 26 - WRC parameters depend on lignin fraction
        a_m_wrc ~ -a_wrc1 * exp(-a_wrc2 * f_LIGN),

        b_m_wrc ~ b_wrc1 + b_wrc2 * f_LIGN,

        # Eq. 26 - Saturated gravimetric water content
        θ_grav_sat ~ a_sat1 * exp(-a_sat2 * f_LIGN),

        # Eq. 26 - Water characteristic: h = a_m * (θ*ρ_n/ρ_m)^(-b_m)
        # θ here is volumetric; gravimetric = θ * ρ_n / ρ_m
        h_Pa ~ a_m_wrc * (max(θ_grav, 0.01))^(-b_m_wrc),
    ]

    return System(eqs, t; name)
end
