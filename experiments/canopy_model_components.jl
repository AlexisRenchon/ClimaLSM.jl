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
LAI = FT(8.0)
z_0m = FT(3.0) # 10% of a 30m canopy
z_0b = FT(3.0)
shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(LAI, z_0m, z_0b, earth_param_set)
lat = FT(0.0)
long = FT(-180)

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

liquid_precip = (t) -> eltype(t)(0)
snow_precip = (t) -> eltype(t)(0)
T_atmos = t -> eltype(t)(290)
q_atmos = t -> eltype(t)(0.01) # kg/kg
P_atmos = t -> eltype(t)(1e5) # Pa
h_atmos = FT(2)
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
RAI = FT(1) # m2/m2
SAI = FT(1) # m2/m2
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
compartment_midpoints = Vector(
    range(
        start = Δz / 2,
        step = Δz,
        stop = Δz * (n_stem + n_leaf) - (Δz / 2),
    ),
)
compartment_surfaces =
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
transpiration =
    PrescribedTranspiration{FT}((t::FT) -> leaf_transpiration(t))
root_extraction = PrescribedSoilPressure{FT}((t::FT) -> ψ_soil0)

plant_hydraulics = PlantHydraulics.PlantHydraulicsModel{FT}(;
                                                            domain = plant_hydraulics_domain,
                                                            parameters = param_set,
                                                            root_extraction = root_extraction,
                                                            transpiration = transpiration,
                                                            root_depths = root_depths,
                                                            n_stem = n_stem,
                                                            n_leaf = n_leaf,
                                                            compartment_surfaces = compartment_surfaces,
                                                            compartment_midpoints = compartment_midpoints,
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

# 1. How to update the auxiliary variables for photosynthesis, RT, stomatal conductance (src) [X]
# 2. Ingest the drivers (atmospheric and radiation) (src) [X]
# 3. DiagnosticTranspiration -> sets the boundary condition for the hydraulics model with the right value (src)
@show propertynames(Y.canopy)
@show Y.canopy.hydraulics
@show propertynames(p.canopy)
@show p.canopy.hydraulics
@show p.canopy.conductance
@show p.canopy.photosynthesis
@show p.canopy.radiative_transfer

#C_d * |u| = g_ae -> g_eff
ClimaLSM.Canopy.canopy_surface_fluxes(canopy.atmos, canopy, Y, p, FT(0.0))