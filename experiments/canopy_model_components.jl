using ClimaCore
using Insolation
using Dates
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM: PrescribedAtmosphere, PrescribedRadiativeFluxes
using ClimaLSM.Canopy
using ClimaLSM.Canopy.PlantHydraulics
using ClimaLSM.Domains: Point
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

FT = Float64 
domain = Point(; z_sfc = FT(0.0))

RTparams = BeerLambertParameters{FT}()
photosynthesis_params = FarquharParameters{FT}(C3();)
stomatal_g_params = MedlynConductanceParameters{FT}()

stomatal_model = MedlynConductanceModel{FT}(stomatal_g_params)
photosynthesis_model = FarquharModel{FT}(photosynthesis_params)
rt_model = BeerLambertModel{FT}(RTparams)

earth_param_set = create_lsm_parameters(FT)
LAI = FT(8.0) # m2 [leaf] m-2 [ground]
z_0m = FT(2.0) # m, Roughness length for momentum - value from tall forest ChatGPT
z_0b = FT(0.1) # m, Roughness length for scalars - value from tall forest ChatGPT
h_c = FT(20.0) # m, canopy height
h_sfc = FT(20.0) # m, canopy height 
h_int = FT(30.0) # m, "where measurements would be taken at a typical flux tower of a 20m canopy"
shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(LAI, h_c, z_0m, z_0b, earth_param_set)
lat = FT(0.0) # degree
long = FT(-180) # degree

function zenith_angle(t::FT; latitude = lat, longitude = long, insol_params = earth_param_set.insol_params) where {FT}    
    return FT(instantaneous_zenith_angle(DateTime(t), longitude, latitude, insol_params)[1])
end

function shortwave_radiation(t::FT; latitude = lat, longitude = long, insol_params = earth_param_set.insol_params) where {FT}
    #θs = FT(instantaneous_zenith_angle(DateTime(t), longitude, latitude, insol_params)[1])
    return FT(1000) # W/m^2
end

function longwave_radiation(t::FT) where {FT}
    return FT(200) # W/m^2
end

u_atmos = t -> eltype(t)(3)

liquid_precip = (t) -> eltype(t)(0) # m
snow_precip = (t) -> eltype(t)(0) # m
T_atmos = t -> eltype(t)(290) # Kelvin
q_atmos = t -> eltype(t)(0.001) # kg/kg
P_atmos = t -> eltype(t)(1e5) # Pa
h_atmos = h_int # m
c_atmos = (t) -> eltype(t)(4.11e-4) # mol/mol
atmos = PrescribedAtmosphere(
            liquid_precip, 
            snow_precip,
            T_atmos,
            u_atmos,
            q_atmos,
            P_atmos,
            c_atmos,
            h_atmos,
        )
radiation = PrescribedRadiativeFluxes(FT, shortwave_radiation, longwave_radiation, zenith_angle)

# Plant Hydraulics
RAI = FT(1) # root area index, m2/m2
SAI = FT(1) # stem area index, m2/m2 TODO: this needs to be consistent with the radiative transfer model
# we already defined LAI, above
area_index = (root = RAI, stem = SAI, leaf = LAI)
K_sat_plant = 1.8e-8 # m/s. Typical conductivity range is [1e-8, 1e-5] m/s. See Kumar, 2008 and
# Pierre Gentine's database for total global plant conductance (1/resistance) 
# (https://github.com/yalingliu-cu/plant-strategies/blob/master/Product%20details.pdf)
K_sat_root = FT(K_sat_plant) # m/s
K_sat_stem = FT(K_sat_plant)
K_sat_leaf = FT(K_sat_plant)
K_sat = (root = K_sat_root, stem = K_sat_stem, leaf = K_sat_leaf)
plant_vg_α = FT(0.002) # 1/m
plant_vg_n = FT(4.2) # unitless
plant_vg_m = FT(1) - FT(1) / plant_vg_n
plant_ν = FT(0.7) # m3/m3
plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
root_depths = -Array(10:-1:1.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0 # 1st element is the deepest root depth 
function root_distribution(z::T) where {T}
    return T(1.0 / 0.5) * exp(z / T(0.5)) # (1/m)
end
Δz = FT(1.0) # height of compartments
n_stem = Int64(0) # number of stem elements
n_leaf = Int64(1) # number of leaf elements
compartment_centers = Vector(
    range(
        start = Δz / 2,
        step = Δz,
        stop = Δz * (n_stem + n_leaf) - (Δz / 2),
    ),
)
compartment_faces =
    Vector(range(start = 0.0, step = Δz, stop = Δz * (n_stem + n_leaf)))
earth_param_set = create_lsm_parameters(FT)

plant_hydraulics_domain = domain
param_set = PlantHydraulics.PlantHydraulicsParameters{
    FT,
    typeof(earth_param_set),
}(
    area_index,
    K_sat,
    plant_vg_α,
    plant_vg_n,
    plant_vg_m,
    plant_ν,
    plant_S_s,
    root_distribution,
    earth_param_set,
)
function leaf_transpiration(t::FT) where {FT}
    T = FT(0)
end

ψ_soil0 = FT(0.0)
transpiration =DiagnosticTranspiration{FT}()
    #PrescribedTranspiration{FT}((t::FT) -> leaf_transpiration(t))
root_extraction = PrescribedSoilPressure{FT}((t::FT) -> ψ_soil0)

plant_hydraulics = PlantHydraulics.PlantHydraulicsModel{FT}(;
                                                            domain = plant_hydraulics_domain,
                                                            parameters = param_set,
                                                            root_extraction = root_extraction,
                                                            transpiration = transpiration,
                                                            root_depths = root_depths,
                                                            n_stem = n_stem,
                                                            n_leaf = n_leaf,
                                                            compartment_surfaces = compartment_faces,
                                                            compartment_midpoints = compartment_centers,
                                                            )
canopy = ClimaLSM.Canopy.CanopyModel{FT}(; parameters = shared_params,
                                         domain = domain,
                                         radiative_transfer = rt_model,
                                         photosynthesis = photosynthesis_model,
                                         conductance = stomatal_model,
                                         hydraulics = plant_hydraulics,
                                         atmos = atmos,
                                         radiation = radiation)
Y,p,coords = ClimaLSM.initialize(canopy)
Y.canopy.hydraulics[1] = plant_ν
update_aux! = make_update_aux(canopy)
t0 = FT(0.0)
update_aux!(p,Y,t0)


# Check values:
# If we write a unit test for this, should we give a range?
parent(p.canopy.photosynthesis.An) * 1e6 # should be around 10 (umol CO2 m-2 s-1)
parent(p.canopy.photosynthesis.GPP) * 1e6





# 1. How to update the auxiliary variables for photosynthesis, RT, stomatal conductance (src) [X]
# 2. Ingest the drivers (atmospheric and radiation) (src) [X]
# 3. DiagnosticTranspiration -> sets the boundary condition for the hydraulics model with the right value (src) [x]
@show propertynames(Y.canopy)
@show propertynames(Y.canopy.hydraulics),Y.canopy.hydraulics

@show propertynames(p.canopy)
@show propertynames(p.canopy.hydraulics), p.canopy.hydraulics
@show propertynames(p.canopy.conductance),p.canopy.conductance
@show propertynames(p.canopy.photosynthesis),p.canopy.photosynthesis
@show propertynames(p.canopy.radiative_transfer),p.canopy.radiative_transfer


# The function needed for timestepping
ode_function! = make_ode_function(canopy)
dY = similar(Y)
ode_function!(dY,Y,p,t0)

#TODO:
 # carefully test all of the pipes and make sure everything is updated when we think it is   
