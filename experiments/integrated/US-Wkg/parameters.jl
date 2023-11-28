global parsed_values
global LAI_data
##PARAMETERS FOR CALIBRATION
g1 = FT(80) # Wang et al: 141 sqrt(Pa) for Medlyn model; Natan used 300.
g0 = FT(FT(0.0001))
soil_K_sat = FT(2e-7) # m/s, matches Natan
soil_vg_n = FT(2.05) # unitless
soil_vg_α = FT(0.02) # inverse meters
K_sat_plant = 5e-9 # m/s # seems much too small?
rooting_depth = FT(1) # from Natan
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
global dt = FT(50)
global n = FT(1) #

# soil domain
# For soil column
global nelements = 60#number of grids in the domain
global zmin = FT(-3)# max depth of the domain
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
global SWC
# Soil parameters
soil_ν = FT(maximum(SWC).+0.05) # m3/m3
soil_S_s = FT(1e-3) # 1/m, guess
θ_r = FT(0) # m3/m3, from Wang et al. 2021 https://doi.org/10.5194/gmd-14-6741-2021

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

capacity = FT(22) # kg/m^2
ϕ = FT(0.6) #The quantum yied of photosystem II : this parameter changes in simulation and as a function of GPP
Ω = FT(0.7)#clumping index can also be a function of effective LAI and number of trees a coefficient that light extiction coefficient is devided by
θj = FT(0.9)#Curvature parameter, a fitting constant to compute  𝐽, this has been reported 0.7 in other works
f = FT(0.015)#Constant factor appearing the dark respiration term
Drel = FT(1.6)# relative diffusivity
Γstar25 = FT(4.275e-5)#mol/mol, CO2 compensation point, it is to compensate for photorespiration, this parameter is constant and similar to other literature
Kc25 = FT(4.049e-4)#mol/mol This is a constant value: half rate constant of carboxylase in leaf at 25 dergree C
#Kc the lower it is the higher the GPP rate will be, as it shows with lower CO2 concentration leaf reachees its maximum potentia
Ko25 = FT(0.2874)#mol/mol This is a constant value:half rate constant of oxylase in leaf at 25 dergree C
#Rubisco can either enter oxygenation process or carboxylation, if oxygenation, O2 is added to Rubisco, the process is named as photo respiration and consumes energy
# a higher value of ko indicates that it takes a higher internal concentration of O2 to reach the photorespiration process
oi = FT(0.209)#Intercellular 𝑂2 concentration unit mol/mol
#sc 5e-6 causes a faster transition
#sc very small causes values of close to 1 for all points
#beta =(1 ) ./( exp.(-sc.*(pc .- pl)) .+ 1);
#beta = (1 + exp(sc*pc)) ./( exp.(-sc.*(pc .- pl)) .+ 1)
To = FT(298.15) 





##########
# model parameters varying for different sites
#########
############################
#parameters from dataset
############################
# Soil parameters
global SWC
soil_ν = FT(round(maximum(SWC),digits=2)) # m3/m3
println("saturated water content is $soil_ν")
soil_K_sat = FT(4.1e-9) # m/s, matches Natan



############################
# parameters with high impact that are tuned
############################
#beta the moisture stress factor
sc = FT(2e-6)#FT(2e-6)#FT(2e6) # Bonan's book: range of 2-5e-6 Pa^{-1}
SAI = FT(maxLAI/4) # m2/m2 or: estimated from Wang et al, FT(0.00242) ?
# Plant Hydraulics and general plant parameters
f_root_to_shoot = FT(1)
RAI = (SAI + maxLAI) * f_root_to_shoot # CLM

plant_ν = capacity / (maxLAI / 2 * h_leaf + SAI * h_stem) / FT(1000)
plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
#to be checked
soil_S_s = FT(1e-3) # 1/m, guess


rooting_depth = FT(7) # from Natan in m
ψ63 = FT(-3 / 0.0098) # / MPa to m, Holtzman's original parameter value is -4 MPa
pc = FT(-1e6) # Bonan's book: -2e6 (Pa): this is a threshols pressure that if leaf potential drops below it, its photosynthesis will be affected by droughts
K_sat_plant = 5e-8 # m/s # seems much too small?
a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
Weibull_param = FT(4) # unitless, Holtzman's original c param value
soil_vg_n = FT(1.28) # unitless
soil_vg_α = FT(1.9) # inverse meters


=#