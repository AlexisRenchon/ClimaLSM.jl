"""Site-specific model parameters for running CliMA LSM on the Harvard Forest
fluxtower site."""

# Data File download link
data_link = "https://caltech.box.com/shared/static/xixaod6511cutz51ag81k1mtvy05hbol.csv"

# Timezone (offset from UTC in hrs)
time_offset = 5

# Site latitude and longitude
lat = FT(42.5378) # degree
long = FT(-72.1715) # degree

# Height of sensor
atmos_h = FT(30)
# https://atmos.seas.harvard.edu/research-harvard_forest-instrumentation

## LAI data from MODIS
# Heterotrophic respiration parameters
# DAMM
α_sx = FT(194e3)
Ea_sx = FT(61e3)
kM_sx = FT(5e-3)
kM_o2 = FT(0.004)
O2_a = FT(0.209)
p_sx = FT(0.024)

# Autotrophic respiration parameters
ηsl = FT(0.01)
μr = FT(1.0)
μs = FT(0.1)

# Soil parameters
soil_ν = FT(0.5) # m3/m3
soil_K_sat = FT(4e-7) # m/s, matches Natan
soil_S_s = FT(1e-3) # 1/m, guess
soil_vg_n = FT(2.05) # unitless
soil_vg_α = FT(0.04) # inverse meters
θ_r = FT(0.067) # m3/m3, from Wang et al. 2021 https://doi.org/10.5194/gmd-14-6741-2021

# Soil heat transfer parameters; not needed for hydrology only test
ν_ss_quartz = FT(0.1)
ν_ss_om = FT(0.1)
ν_ss_gravel = FT(0.0);
soil_ϵ = FT(0.98)
soil_α_PAR = FT(0.2)
soil_α_NIR = FT(0.2)

# TwoStreamModel parameters
Ω = FT(0.69)
α_PAR_leaf = FT(0.1)
τ_PAR_leaf = FT(0.05)
α_NIR_leaf = FT(0.45)
τ_NIR_leaf = FT(0.25)
ϵ_canopy = FT(0.97)

# Energy Balance model
ac_canopy = FT(2.5e3)

# Conductance Model
g1 = FT(141) # Wang et al: 141 sqrt(Pa) for Medlyn model; Natan used 300.
g0 = FT(1e-4)

#Photosynthesis model
oi = FT(0.209)
ϕ = FT(0.6)
θj = FT(0.9)
sc = FT(2e-6) # Bonan's book: range of 2-5e-6
pc = FT(-2e6) # Bonan's book: -2e6
Vcmax25 = FT(9e-5) # from Yujie's paper 4.5e-5 , Natan used 9e-5
Γstar25 = FT(4.275e-5)
Kc25 = FT(4.049e-4)
Ko25 = FT(0.2874)

# Plant Hydraulics and general plant parameters
SAI = FT(1.0) # m2/m2 or: estimated from Wang et al, FT(0.00242) ?
f_root_to_shoot = FT(3.5)
K_sat_plant = 5e-9 # m/s # seems much too small?
ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value is -4 MPa
Weibull_param = FT(4) # unitless, Holtzman's original c param value
a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
conductivity_model =
    PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
capacity = FT(10) # kg/m^2
plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
rooting_depth = FT(0.5) # from Natan
z0_m = FT(0.13) * h_canopy
z0_b = FT(0.1) * z0_m
