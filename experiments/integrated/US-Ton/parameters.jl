
global parsed_values
global LAI_data
##PARAMETERS FOR CALIBRATION
g1 = FT(100) # Wang et al: 141 sqrt(Pa) for Medlyn model; Natan used 300.
g0 = FT(FT(0.001))
soil_K_sat = FT(4e-7) # m/s, matches Natan
soil_vg_n = FT(2.05) # unitless
soil_vg_α = FT(0.04) # inverse meters
K_sat_plant = 5e-9 # m/s # seems much too small?
rooting_depth = FT(0.5) # from Natan
pc = FT(-2e6) # Bonan's book: -2e6
global Vcmax25=quantile(skipmissing(parsed_values),0.95)*10^-6;
ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value is -4 MPa
Weibull_param = FT(4) # unitless, Holtzman's original c param value
# domain and time step parameters

global end_date
global start_date
global spinup
global nspinup
if spinup==1
global t0 = FT((start_date +350*nspinup)* 3600 * 24)# start mid year
else
global t0 = FT(1* 3600 * 24)# start mid year
end
global N_days = end_date-start_date
global tf = t0 + FT(3600 * 24 * N_days)
global dt = FT(60)
global n = FT(1) #

# soil domain
# For soil column
global nelements = 60#number of grids in the domain
global zmin = FT(-10)# max depth of the domain
global zmax = FT(0)# depth at the surface,
global dz_bottom = FT(0.5)
global dz_top = FT(0.025)
# Number of stem and leaf compartments. Leaf compartments are stacked on top of stem compartments
global n_stem = Int64(1);
global n_leaf = Int64(1);
global h_stem = FT(1) # m, from Wang et al.
global h_leaf = FT(6.5) # m from Wang et al.



# Autotrophic respiration parameters
ne = FT(8 * 1e-4)
ηsl = FT(0.01)
σl = FT(0.05)
μr = FT(1.0)
μs = FT(0.1)
f1 = FT(0.012)
f2 = FT(0.25)

# Soil parameters
soil_ν = FT(0.55) # m3/m3
soil_S_s = FT(1e-3) # 1/m, guess
θ_r = FT(0.067) # m3/m3, from Wang et al. 2021 https://doi.org/10.5194/gmd-14-6741-2021

# Soil heat transfer parameters; not needed for hydrology only test
ν_ss_quartz = FT(0.1)
ν_ss_om = FT(0.1)
ν_ss_gravel = FT(0.0);
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
z_0m_soil = FT(0.1)
z_0b_soil = FT(0.1)
soil_ϵ = FT(0.98)
soil_α_PAR = FT(0.2)
soil_α_NIR = FT(0.2)

# TwoStreamModel parameters
Ω = FT(0.69)
ld = FT(0.5)#Leaf angle distribution no unit, this is Xl in the text book, but it actually does not matter much
α_PAR_leaf = FT(0.1)
λ_γ_PAR = FT(5e-7)
λ_γ_NIR = FT(1.65e-6)
τ_PAR_leaf = FT(0.05)
α_NIR_leaf = FT(0.45)
τ_NIR_leaf = FT(0.25)
n_layers = UInt64(20)
ϵ_canopy = FT(0.97)

# Conductance Model

Drel = FT(1.6)


#Photosynthesis model
oi = FT(0.209)
ϕ = FT(0.6)
θj = FT(0.9)
f = FT(0.015)
sc = FT(2e-6) # Bonan's book: range of 2-5e-6

Γstar25 = FT(4.275e-5)
Kc25 = FT(4.049e-4)
Ko25 = FT(0.2874)
To = FT(298.15)#Standard temperature
ΔHkc = FT(79430)
ΔHko = FT(36380)
ΔHVcmax = FT(58520)
ΔHΓstar = FT(37830)
ΔHJmax = FT(43540)
ΔHRd = FT(46390)

# Plant Hydraulics and general plant parameters
maxLAI = FT(quantile(skipmissing(LAI_data.LAI),0.95)) # m2/m2, from Wang et al.
SAI = FT(maxLAI/4) 
f_root_to_shoot = FT(3.5)
RAI = (SAI + maxLAI) * f_root_to_shoot # CLM
a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
conductivity_model =
    PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
capacity = FT(10) # kg/m^2
plant_ν = capacity / (maxLAI / 2 * h_leaf + SAI * h_stem) / FT(1000)
plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
z0_m = FT(2)
z0_b = FT(0.2)


#=
global parsed_values
global LAI_data
##PARAMETERS FOR CALIBRATION
g1 = FT(141) # Wang et al: 141 sqrt(Pa) for Medlyn model; Natan used 300.
g0 = FT(1e-4)
soil_K_sat = FT(4e-7) # m/s, matches Natan
soil_vg_n = FT(2.6257) # unitless
soil_vg_α = FT(1.368) # inverse meters
K_sat_plant = 1.8e-9 # m/s # seems much too small?
rooting_depth = FT(1) # from Natan
pc = FT(-2e5)# Bonan's book: -2e6
global Vcmax25=FT(5e-5)#global Vcmax25=quantile(skipmissing(parsed_values),0.95)*10^-6;
ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value is -4 MPa
Weibull_param = FT(4) # unitless, Holtzman's original c param value
# domain and time step parameters

global end_date
global start_date
global spinup
global nspinup
if spinup==1
global t0 = FT((start_date +350*nspinup)* 3600 * 24)# start mid year
else
global t0 = FT(7* 3600 * 24)# start mid year
end
global N_days = end_date-start_date
global tf = t0 + FT(3600 * 24 * N_days)
global dt = FT(70)
global n = FT(1) 

# soil domain
# For soil column
global nelements = 10#number of grids in the domain
global zmin = FT(-2)# max depth of the domain
global zmax = FT(0)# depth at the surface,
global dz_bottom = FT(0.5)
global dz_top = FT(0.025)
# Number of stem and leaf compartments. Leaf compartments are stacked on top of stem compartments
global n_stem = Int64(1);
global n_leaf = Int64(1);
global h_stem = FT(9)  # m, from Wang et al.
global h_leaf = FT(9) # m from Wang et al.



# Autotrophic respiration parameters
ne = FT(8 * 1e-4)
ηsl = FT(0.01)
σl = FT(0.05)
μr = FT(1.0)
μs = FT(0.1)
f1 = FT(0.012)
f2 = FT(0.25)

# Soil parameters
soil_ν = FT(maximum(SWC)) # m3/m3
soil_S_s = FT(1e-3) # 1/m, guess
θ_r = FT(0.067) # m3/m3, from Wang et al. 2021 https://doi.org/10.5194/gmd-14-6741-2021

# Soil heat transfer parameters; not needed for hydrology only test
ν_ss_quartz = FT(0.1)
ν_ss_om = FT(0.1)
ν_ss_gravel = FT(0.0);
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
z_0m_soil = FT(0.1)
z_0b_soil = FT(0.1)
soil_ϵ = FT(0.98)
soil_α_PAR = FT(0.2)
soil_α_NIR = FT(0.2)

# TwoStreamModel parameters
Ω = FT(0.69)
ld = FT(0.5)#Leaf angle distribution no unit, this is Xl in the text book, but it actually does not matter much
α_PAR_leaf = FT(0.11)
λ_γ_PAR = FT(4e-7)#FT(5e-7)#
λ_γ_NIR = FT(1.65e-6)
τ_PAR_leaf = FT(0.05)
α_NIR_leaf = FT(0.35)
τ_NIR_leaf = FT(0.34)
#n_layers = UInt64(20)
n_layers = UInt64(5)
ϵ_canopy = FT(0.97)



# Conductance Model

Drel = FT(1.6)


#Photosynthesis model
oi = FT(0.209)
ϕ = FT(0.6)
θj = FT(0.9)
f = FT(0.015)
sc = FT(2e-6) # Bonan's book: range of 2-5e-6

Γstar25 = FT(4.275e-5)
Kc25 = FT(4.049e-4)
Ko25 = FT(0.2874)
To = FT(298.15)#Standard temperature
ΔHkc = FT(79430)
ΔHko = FT(36380)
ΔHVcmax = FT(58520)
ΔHΓstar = FT(37830)
ΔHJmax = FT(43540)
ΔHRd = FT(46390)

# Plant Hydraulics and general plant parameters
#maxLAI = FT(quantile(skipmissing(LAI_data.LAI),0.95)) # m2/m2, from Wang et al.
maxLAI = FT(maximum(LAI_data.LAI))
SAI = 1#FT(maxLAI/4) 
f_root_to_shoot = FT(3.5)
RAI = (SAI + maxLAI) * f_root_to_shoot # CLM
a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
conductivity_model =
    PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)
capacity = FT(7) # kg/m^2
plant_ν = capacity / (maxLAI / 2 * h_leaf + SAI * h_stem) / FT(1000)
plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
z0_m = FT(2)
z0_b = FT(0.2)
=#

