using SurfaceFluxes
using Thermodynamics
using StaticArrays

# Note that this assumes that earth_param_ste and LSMP have been defined
# already!
thermo_params = LSMP.thermodynamic_parameters(earth_param_set)
surface_flux_params = LSMP.surface_fluxes_parameters(earth_param_set)

h_canopy = FT(18.5)
z_0m = FT(0.1)*h_canopy # Bonan's book
z_0b = FT(0.1)*z_0m # Bonan's book
d_sfc = FT(0.67)*h_canopy # Bonan's book
_σ = LSMP.Stefan(earth_param_set)

dleaf = 0.04 # CLM Table 2.5.1
Cv = 0.01 # CLM right after 2.5.122
Cs_canopy = 0.004 # CLM 2.5.120

# Define the functions that compute the surface fluxes.
function atmos_airspace_fluxes(ts_sfc, ts_in, atmos_h, d_sfc, u_air,z_0m, z_0b, surface_flux_params)
    state_sfc =
        SurfaceFluxes.SurfaceValues(FT(0), SVector{2, FT}(0, 0), ts_sfc)
    state_in =
        SurfaceFluxes.InteriorValues(atmos_h - d_sfc, SVector{2, FT}(u_air, 0), ts_in)
        sc = SurfaceFluxes.ValuesOnly{FT}(;
                                          state_in,
                                          state_sfc,
                                          z0m = z_0m,
                                          z0b = z_0b,
                                          )
        
        # Canopy airspace at z0+d and atmos at h
        conditions = SurfaceFluxes.surface_conditions(
            surface_flux_params,
            sc;
            tol_neutral = SFP.cp_d(surface_flux_params) / 10000,
        )
    cp_d = Thermodynamics.Parameters.cp_d(thermo_params)
    grav = SurfaceFluxes.Parameters.grav(surface_flux_params)
    R_d = Thermodynamics.Parameters.R_d(thermo_params)
    T_0 = LSMP.T_0(earth_param_set)
    T_sfc = Thermodynamics.air_temperature(thermo_params, ts_sfc)
    hd_sfc = cp_d * (T_sfc - T_0) + R_d * T_0
    E0 = SurfaceFluxes.evaporation(surface_flux_params, sc, conditions.Ch)
    r_a = 1 / (conditions.Ch * SurfaceFluxes.windspeed(sc)) # definition
    return (r_a = r_a, ustar = conditions.ustar, shf = conditions.shf + hd_sfc*E0, lhf = conditions.lhf- hd_sfc*E0, evaporation = conditions.evaporation, L_MO = conditions.L_MO)
end

function airspace_canopy_fluxes(ustar, T_airspace, q_airspace, T_canopy, ρ_air, cp_m, LAI, AI, r_canopy)
    r_b = 1/Cv/sqrt(ustar/dleaf)^0.5 # CLM Equation 2.5.122
    q_canopy =  
        Thermodynamics.q_vap_saturation_generic(
            thermo_params,
            T_canopy,
            ρ_air,
            Thermodynamics.Liquid(),
        )
    ΔT = (T_airspace - T_canopy)
    canopy_shf  = -ρ_air * cp_m * ΔT / (r_b/AI)
    canopy_evap = -ρ_air * (q_airspace - q_canopy) / (r_canopy + r_b/LAI)
    return (;shf = canopy_shf, transpiration = canopy_evap)
end
function land_lw_fluxes(LW_d, T_canopy, T_soil, ϵ_soil, ϵ_canopy, _σ)
    LW_d_canopy = (1 - ϵ_canopy) * LW_d + ϵ_canopy * _σ * T_canopy^4
    LW_u_soil = ϵ_soil * _σ * T_soil^4 + (1 - ϵ_soil) * LW_d_canopy
    LW_canopy = ϵ_canopy * LW_d - 2 * ϵ_canopy * _σ * T_canopy^4 + ϵ_canopy * LW_u_soil
    LW_out = (1 - ϵ_canopy) * LW_u_soil + ϵ_canopy * _σ * T_canopy^4
    LW_soil = ϵ_soil * LW_d_canopy - ϵ_soil * _σ * T_soil^4
    return (;LW_canopy = LW_canopy, LW_out = LW_out, LW_soil = LW_soil)
end

function airspace_soil_fluxes(ustar, T_airspace, q_airspace, W, T_soil, q_soil, r_soil, ρ_air, cp_m, ΔΦ)
    Cs_bare = 0.4/0.13*(0.01*ustar/(1.5e-5))^(-0.45) # CLM Equation 2.5.121
    Cs = @. W*Cs_bare + (1-W)*Cs_canopy # CLM Equation 2.5.118
    r_ah_cs = @. 1/Cs/ustar # CLm Eq. 2.5.116 and 2.5.117
    ΔT = (T_airspace - T_soil)
    soil_shf = -ρ_air * (cp_m * ΔT + ΔΦ) / r_ah_cs
    soil_evap = -ρ_air * (q_airspace - q_soil) / (r_ah_cs + r_soil)
    return (;shf = soil_shf, evaporation = soil_evap)
end


function atmos_sfc_fluxes(ts_sfc, ts_in, atmos_h, d_sfc, u_air,z_0m, z_0b, surface_flux_params, thermo_params, r_sfc)
    state_sfc =
        SurfaceFluxes.SurfaceValues(FT(0), SVector{2, FT}(0, 0), ts_sfc)
    state_in =
        SurfaceFluxes.InteriorValues(atmos_h - d_sfc, SVector{2, FT}(u_air, 0), ts_in)
    sc = SurfaceFluxes.ValuesOnly{FT}(;
                                      state_in,
                                      state_sfc,
                                      z0m = z_0m,
                                      z0b = z_0b,
                                      )
    
    # Surface at z0+d and atmos at h
    conditions = SurfaceFluxes.surface_conditions(
        surface_flux_params,
        sc;
        tol_neutral = SFP.cp_d(surface_flux_params) / 10000,
    )
    cp_d = Thermodynamics.Parameters.cp_d(thermo_params)
    R_d = Thermodynamics.Parameters.R_d(thermo_params)
    T_0 = LSMP.T_0(earth_param_set)
    T_sfc = Thermodynamics.air_temperature(thermo_params, ts_sfc)
    hd_sfc = cp_d * (T_sfc - T_0) + R_d * T_0
    E0 = SurfaceFluxes.evaporation(surface_flux_params, sc, conditions.Ch)
    r_a = 1 / (conditions.Ch * SurfaceFluxes.windspeed(sc)) # definition
    E = E0 * r_a / (r_sfc + r_a)
    H = conditions.shf + hd_sfc * E0
    LH = conditions.lhf * r_a / (r_sfc + r_a) - hd_sfc*E

    return (r_a = r_a, ustar = conditions.ustar, shf = H, lhf = LH, evaporation = E, L_MO = conditions.L_MO)
end
