using DiffEqCallbacks
import OrdinaryDiffEq as ODE
import ClimaTimeSteppers as CTS
using ClimaCore
import CLIMAParameters as CP
using Statistics
using Dates
using Insolation

using ClimaLSM
using ClimaLSM.Domains: LSMSingleColumnDomain
using ClimaLSM.Soil
using ClimaLSM.Canopy
using ClimaLSM.Canopy.PlantHydraulics
import ClimaLSM
import ClimaLSM.Parameters as LSMP

const FT = Float64

include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))
earth_param_set = create_lsm_parameters(FT)
climalsm_dir = pkgdir(ClimaLSM)

function ClimaLSM_ozark(g1)

#=
# This reads in the data from the flux tower site and creates
# the atmospheric and radiative driver structs for the model
include(
    joinpath(
        climalsm_dir,
        "experiments/web_dashboards/ozark_met_drivers_FLUXNET.jl",
    ),
)
include(joinpath(climalsm_dir, "experiments/web_dashboards/ozark_domain.jl"))
include(
    joinpath(climalsm_dir, "experiments/web_dashboards/ozark_parameters.jl"),
)
include(
    joinpath(climalsm_dir, "experiments/web_dashboards/ozark_simulation.jl"),
)
=#

# Now we set up the model. For the soil model, we pick
# a model type and model args:
soil_domain = land_domain.subsurface
soil_ps = Soil.EnergyHydrologyParameters{FT}(;
    κ_dry = FT(κ_dry_val),
    κ_sat_frozen = FT(κ_sat_frozen_val),
    κ_sat_unfrozen = FT(κ_sat_unfrozen_val),
    ρc_ds = ρc_ds,
    ν = soil_ν,
    ν_ss_om = ν_ss_om,
    ν_ss_quartz = ν_ss_quartz,
    ν_ss_gravel = ν_ss_gravel,
    hydrology_cm = vanGenuchten(; α = soil_vg_α, n = soil_vg_n),
    K_sat = soil_K_sat,
    S_s = soil_S_s,
    θ_r = θ_r,
    earth_param_set = earth_param_set,
    z_0m = z_0m_soil,
    z_0b = z_0b_soil,
    emissivity = soil_ϵ,
    PAR_albedo = soil_α_PAR,
    NIR_albedo = soil_α_NIR,
);

soil_args = (domain = soil_domain, parameters = soil_ps)
soil_model_type = Soil.EnergyHydrology{FT}

# Now we set up the canopy model, which we set up by component:
# Component Types
canopy_component_types = (;
    radiative_transfer = Canopy.TwoStreamModel{FT},
    photosynthesis = Canopy.FarquharModel{FT},
    conductance = Canopy.MedlynConductanceModel{FT},
    hydraulics = Canopy.PlantHydraulicsModel{FT},
)
# Individual Component arguments
# Set up radiative transfer
radiative_transfer_args = (;
    parameters = TwoStreamParameters{FT}(;
        Ω = Ω,
        ld = ld,
        α_PAR_leaf = α_PAR_leaf,
        λ_γ_PAR = λ_γ_PAR,
        λ_γ_NIR = λ_γ_NIR,
        τ_PAR_leaf = τ_PAR_leaf,
        α_NIR_leaf = α_NIR_leaf,
        τ_NIR_leaf = τ_NIR_leaf,
        n_layers = n_layers,
        diff_perc = diff_perc,
    )
)
# Set up conductance
conductance_args = (;
    parameters = MedlynConductanceParameters{FT}(;
        g1 = g1,
        Drel = Drel,
        g0 = g0,
    )
)
# Set up photosynthesis
photosynthesis_args = (;
    parameters = FarquharParameters{FT}(
        Canopy.C3();
        oi = oi,
        ϕ = ϕ,
        θj = θj,
        f = f,
        sc = sc,
        pc = pc,
        Vcmax25 = Vcmax25,
        Γstar25 = Γstar25,
        Kc25 = Kc25,
        Ko25 = Ko25,
        To = To,
        ΔHkc = ΔHkc,
        ΔHko = ΔHko,
        ΔHVcmax = ΔHVcmax,
        ΔHΓstar = ΔHΓstar,
        ΔHJmax = ΔHJmax,
        ΔHRd = ΔHRd,
    )
)
# Set up plant hydraulics
ai_parameterization = PrescribedSiteAreaIndex{FT}(LAIfunction, SAI, RAI)

function root_distribution(z::T; rooting_depth = rooting_depth) where {T}
    return T(1.0 / rooting_depth) * exp(z / T(rooting_depth)) # 1/m
end

plant_hydraulics_ps = PlantHydraulics.PlantHydraulicsParameters(;
    ai_parameterization = ai_parameterization,
    ν = plant_ν,
    S_s = plant_S_s,
    root_distribution = root_distribution,
    conductivity_model = conductivity_model,
    retention_model = retention_model,
)
plant_hydraulics_args = (
    parameters = plant_hydraulics_ps,
    n_stem = n_stem,
    n_leaf = n_leaf,
    compartment_midpoints = compartment_midpoints,
    compartment_surfaces = compartment_surfaces,
)

# Canopy component args
canopy_component_args = (;
    radiative_transfer = radiative_transfer_args,
    photosynthesis = photosynthesis_args,
    conductance = conductance_args,
    hydraulics = plant_hydraulics_args,
)
# Other info needed
shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
    z0_m,
    z0_b,
    earth_param_set,
)

canopy_model_args = (; parameters = shared_params, domain = land_domain.surface)

# Integrated plant hydraulics and soil model
land_input = (atmos = atmos, radiation = radiation)
land = SoilCanopyModel{FT}(;
    land_args = land_input,
    soil_model_type = soil_model_type,
    soil_args = soil_args,
    canopy_component_types = canopy_component_types,
    canopy_component_args = canopy_component_args,
    canopy_model_args = canopy_model_args,
)
Y, p, cds = initialize(land)
exp_tendency! = make_exp_tendency(land)

#Initial conditions
Y.soil.ϑ_l = SWC[1 + Int(round(t0 / 1800))] # Get soil water content at t0
# recalling that the data is in intervals of 1800 seconds. Both the data
# and simulation are reference to 2005-01-01-00 (LOCAL)
# or 2005-01-01-06 (UTC)
Y.soil.θ_i = FT(0.0)
T_0 = TS[1 + Int(round(t0 / 1800))] # Get soil temperature at t0
ρc_s =
    volumetric_heat_capacity.(Y.soil.ϑ_l, Y.soil.θ_i, Ref(land.soil.parameters))
Y.soil.ρe_int =
    volumetric_internal_energy.(
        Y.soil.θ_i,
        ρc_s,
        T_0,
        Ref(land.soil.parameters),
    )
ψ_stem_0 = FT(-1e5 / 9800)
ψ_leaf_0 = FT(-2e5 / 9800)

S_l_ini =
    inverse_water_retention_curve.(
        retention_model,
        [ψ_stem_0, ψ_leaf_0],
        plant_ν,
        plant_S_s,
    )

for i in 1:2
    Y.canopy.hydraulics.ϑ_l[i] .=
        augmented_liquid_fraction.(plant_ν, S_l_ini[i])
end

set_initial_aux_state! = make_set_initial_aux_state(land)
set_initial_aux_state!(p, Y, t0);

# Simulation
sv = (;
    t = Array{FT}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
)
cb = ClimaLSM.NonInterpSavingCallback(sv, saveat)

prob =
    ODE.ODEProblem(CTS.ClimaODEFunction(T_exp! = exp_tendency!), Y, (t0, tf), p);
sol = ODE.solve(
    prob,
    ode_algo;
    dt = dt,
    callback = cb,
    adaptive = false,
    saveat = saveat,
)

#daily = sol.t ./ 3600 ./ 24
modeled_GPP = [parent(sv.saveval[k].canopy.photosynthesis.GPP)[1] for k in 1:length(sv.saveval)]

return modeled_GPP
end