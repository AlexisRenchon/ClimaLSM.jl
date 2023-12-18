"""Constant parameters used in the integrated model to run on Fluxnet sites. The 
locations of these parameters should eventually be in CliMA parameters but for 
now are here just to remove them from the site-specific parameter files."""

# Constant thermal properties of the soil

κ_quartz = FT(7.7) # W/m/K
κ_minerals = FT(2.5) # W/m/K
κ_om = FT(0.25) # W/m/K
κ_liq = FT(0.57) # W/m/K
κ_ice = FT(2.29) # W/m/K
κ_air = FT(0.025); #W/m/K
ρp = FT(2700); # kg/m^3
κ_solid = Soil.κ_solid(ν_ss_om, ν_ss_quartz, κ_om, κ_quartz, κ_minerals)
κ_dry = Soil.κ_dry(ρp, soil_ν, κ_solid, κ_air)
κ_sat_frozen = Soil.κ_sat_frozen(κ_solid, soil_ν, κ_ice)
κ_sat_unfrozen = Soil.κ_sat_unfrozen(κ_solid, soil_ν, κ_liq);
ρc_ds = FT((1 - soil_ν) * 4e6); # J/m^3/K
z_0m_soil = FT(0.01)
z_0b_soil = FT(0.001)

# Constants Attributes for Plant Hydraulics
conductivity_model =
    PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
