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
using ClimaLSM.Domains: Point, Plane
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

FT = Float64
domains = [
    Point(; z_sfc = FT(0.0)),
    Plane(;
        xlim = (0.0, 1.0),
        ylim = (0.0, 1.0),
        nelements = (2, 2),
        periodic = (true, true),
        npolynomial = 1,
    ),
]

RTparams = BeerLambertParameters{FT}()
photosynthesis_params = FarquharParameters{FT}(C3();)
stomatal_g_params = MedlynConductanceParameters{FT}()

stomatal_model = MedlynConductanceModel{FT}(stomatal_g_params)
photosynthesis_model = FarquharModel{FT}(photosynthesis_params)
rt_model = BeerLambertModel{FT}(RTparams)

earth_param_set = create_lsm_parameters(FT)
LAI = FT(0.5)
shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(LAI, earth_param_set)
lat = FT(0.0)
long = FT(0.0)
function θs(t::FT; latitude = lat, longitude = long, insol_params = earth_param_set.insol_params) where {FT}
    return FT(instantaneous_zenith_angle(DateTime(t), longitude, latitude, insol_params)[1])
end

function SW_d(t::FT; latitude = lat, longitude = long, insol_params = earth_param_set.insol_params) where {FT}
    θs = FT(instantaneous_zenith_angle(DateTime(t), longitude, latitude, insol_params)[1])
    return cos(θs) * FT(500) # W/m^2
end

function LW_d(t::FT; latitude = lat, longitude = long, insol_params = earth_param_set.insol_params) where {FT}
    θs = FT(instantaneous_zenith_angle(DateTime(t), longitude, latitude, insol_params)[1])
    return cos(θs) * FT(500) # W/m^2
end

u(t) = t -> eltype(t)(3)
function u(t)
    return 3
end

T(t) = t -> eltype(t)(290)
q(t) = t -> eltype(t)()
atmos = PrescribedAtmos{FT}()
radiation = PrescribedRadiativeFluxes{FT}()

# Plant Hydraulics
RAI = FT(1) # m2/m2
SAI = FT(1) # m2/m2
LAI = FT(1) # m2/m2
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

plant_hydraulics_domain = domains[1]
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
args = (shared_params, domains[1], rt_model, photosynthesis_model, stomatal_model, plant_hydraulics)
typeargs = (rt_model, photosynthesis_model, stomatal_model, plant_hydraulics)
canopy = ClimaLSM.Canopy.CanopyModel{FT}(; parameters = shared_params,
                                         domain =domains[1],
                                         radiative_transfer = rt_model,
                                         photosynthesis = photosynthesis_model,
                                         conductance = stomatal_model,
                                         hydraulics = plant_hydraulics)
Y,p,coords = ClimaLSM.initialize(canopy)


# 1. How to update the auxiliary variables for photosynthesis, RT, stomatal conductance (src)
# 2. Ingest the drivers (atmospheric and radiation) (src)



# 3. DiagnosticTranspiration -> sets the boundary condition for the hydraulics model with the right value (src)
@show propertynames(Y.canopy)
@show Y.canopy.hydraulics
@show propertynames(p.canopy)
@show p.canopy.hydraulics
@show p.canopy.conductance
@show p.canopy.photosynthesis
@show p.canopy.radiative_transfer

#C_d * |u| = g_ae -> g_eff
