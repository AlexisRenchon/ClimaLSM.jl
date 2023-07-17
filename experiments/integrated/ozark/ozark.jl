using DiffEqCallbacks
import OrdinaryDiffEq as ODE
import ClimaTimeSteppers as CTS
using ClimaCore
import CLIMAParameters as CP
using Plots
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
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))
const FT = Float64
earth_param_set = create_lsm_parameters(FT)
climalsm_dir = pkgdir(ClimaLSM)
# This reads in the data from the flux tower site and creates
# the atmospheric and radiative driver structs for the model
include(
    joinpath(
        climalsm_dir,
        "experiments/integrated/ozark/ozark_met_drivers_FLUXNET.jl",
    ),
)
include(joinpath(climalsm_dir, "experiments/integrated/ozark/ozark_domain.jl"))
include(
    joinpath(climalsm_dir, "experiments/integrated/ozark/ozark_parameters.jl"),
)
include(
    joinpath(climalsm_dir, "experiments/integrated/ozark/ozark_simulation.jl"),
)

# Now we set up the model. For the soil model, we pick
# a model type and model args:
soil_domain = land_domain.subsurface
soil_ps = Soil.EnergyHydrologyParameters{FT}(;
    κ_dry = κ_dry,
    κ_sat_frozen = κ_sat_frozen,
    κ_sat_unfrozen = κ_sat_unfrozen,
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
    albedo = soil_α,
);

soil_args = (domain = soil_domain, parameters = soil_ps)
soil_model_type = Soil.EnergyHydrology{FT}

# Now we set up the canopy model, which we set up by component:
# Component Types
canopy_component_types = (;
    radiative_transfer = Canopy.BeerLambertModel{FT},
    photosynthesis = Canopy.FarquharModel{FT},
    conductance = Canopy.MedlynConductanceModel{FT},
    hydraulics = Canopy.PlantHydraulicsModel{FT},
)
# Individual Component arguments
# Set up radiative transfer
radiative_transfer_args = (;
    parameters = BeerLambertParameters{FT}(;
        Ω = Ω,
        ld = ld,
        ρ_leaf = ρ_leaf,
        λ_γ = λ_γ,
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
area_index = (root = RAI, stem = SAI, leaf = LAI)

function root_distribution(z::T; rooting_depth = rooting_depth) where {T}
    return T(1.0 / rooting_depth) * exp(z / T(rooting_depth)) # 1/m
end

plant_hydraulics_ps = PlantHydraulics.PlantHydraulicsParameters(;
    area_index = area_index,
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
    LAI,
    h_stem + h_leaf,
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
Y.soil.ϑ_l = SWC[Int(round(t0 / 1800))] # Get soil water content at t0
# recalling that the data is in intervals of 1800 seconds. Both the data
# and simulation are reference to 2005-01-01-00 (LOCAL)
# or 2005-01-01-06 (UTC)
Y.soil.θ_i = FT(0.0)
T_0 = TS[Int(round(t0 / 1800))] # Get soil temperature at t0
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

update_aux! = make_update_aux(land)
update_aux!(p, Y, t0)

# Simulation
sv = (;
    t = Array{FT}(undef, length(saveat)),
    saveval = Array{ClimaCore.Fields.FieldVector}(undef, length(saveat)),
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

# Plotting
daily = sol.t ./ 3600 ./ 24
savedir = joinpath(climalsm_dir, "experiments/integrated/ozark/")

# GPP
model_GPP = [
    parent(sv.saveval[k].canopy.photosynthesis.GPP)[1] for
    k in 1:length(sv.saveval)
]

plt1 = Plots.plot(size = (500, 400))
Plots.plot!(
    plt1,
    daily,
    model_GPP .* 1e6,
    label = "Model",
    xlim = [minimum(daily), maximum(daily)],
)
Plots.plot!(
    plt1,
    seconds ./ 3600 ./ 24,
    GPP .* 1e6,
    label = "Data",
    lalpha = 0.3,
)

plt2 = Plots.plot(size = (500, 400))
Plots.plot!(
    plt2,
    daily,
    model_GPP .* 1e6,
    label = "",
    xlim = [maximum(daily) - 7, maximum(daily)],
    xlabel = "Day of year",
    ylabel = "GPP [μmol/mol]",
)
Plots.plot!(plt2, seconds ./ 3600 ./ 24, GPP .* 1e6, label = "", lalpha = 0.3)
Plots.plot(plt1, plt2, layout = (2, 1))
Plots.savefig(joinpath(savedir, "GPP.png"))

# SW_IN

plt1 = Plots.plot(size = (500, 400))
Plots.plot!(
    plt1,
    seconds ./ 3600 ./ 24,
    FT.(SW_IN),
    label = "Data",
    xlim = [maximum(daily) - 7, maximum(daily)],
)

plt2 = Plots.plot(size = (500, 400))
Plots.plot!(
    plt2,
    seconds ./ 3600 ./ 24,
    FT.(SW_IN),
    label = "",
    xlim = [minimum(daily), maximum(daily)],
    ylabel = " SW radiation (W/m^2)",
    xlabel = "Day of year",
)
Plots.plot(plt1, plt2, layout = (2, 1))

Plots.savefig(joinpath(savedir, "SW_IN.png"))

# VPD

plt1 = Plots.plot(size = (500, 400))
Plots.plot!(
    plt1,
    seconds ./ 3600 ./ 24,
    FT.(VPD),
    label = "Data",
    xlim = [maximum(daily) - 7, maximum(daily)],
)

plt2 = Plots.plot(size = (500, 400))
Plots.plot!(
    plt2,
    seconds ./ 3600 ./ 24,
    FT.(VPD),
    label = "",
    xlim = [minimum(daily), maximum(daily)],
    ylabel = "VPD (Pa)",
    xlabel = "Day of year",
)
Plots.plot(plt1, plt2, layout = (2, 1))
Plots.savefig(joinpath(savedir, "VPD.png"))



T =
    [
        parent(sv.saveval[k].canopy.conductance.transpiration)[1] for
        k in 1:length(sol.t)
    ] .* (1e3 * 24 * 3600)
E =
    [parent(sv.saveval[k].soil_evap)[1] for k in 1:length(sol.t)] .* (1e3 * 24 * 3600)
measured_T = LE ./ (LSMP.LH_v0(earth_param_set) * 1000) .* (1e3 * 24 * 3600)

plt1 = Plots.plot(size = (500, 400))
Plots.plot!(
    plt1,
    seconds ./ 3600 ./ 24,
    measured_T,
    label = "Data",
    margins = 10Plots.mm,
)
Plots.plot!(
    plt1,
    daily,
    T,
    label = "Model T",
    xlim = [minimum(daily), maximum(daily)],
    ylim = [0, 30],
)

Plots.plot!(
    plt1,
    daily,
    E,
    label = "Model E",
    ylim = [0, 30],
    legend = :topright,
)

plt2 = Plots.plot(size = (500, 400))
Plots.plot!(plt2, seconds ./ 3600 ./ 24, measured_T, label = "")
Plots.plot!(
    plt2,
    daily,
    T,
    xlim = [maximum(daily) - 7, maximum(daily)],
    label = "",
    xlabel = "Day of year",
    ylabel = "Vapor Flux [mm/day]",
    ylim = [0, 30],
)

Plots.plot!(plt2, daily, E, label = "", ylim = [0, 30])
Plots.plot(plt1, plt2, layout = (2, 1))
Plots.savefig(joinpath(savedir, "ET.png"))

β = [parent(sv.saveval[k].canopy.hydraulics.β)[1] for k in 1:length(sol.t)]
plt1 = Plots.plot(size = (500, 400))
Plots.plot!(
    plt1,
    daily,
    β,
    label = "Model",
    xlim = [minimum(daily), maximum(daily)],
)
#i_week = Int(round((7 * 24 * 3600) / (sol.t[2] - sol.t[1])))
plt2 = Plots.plot(size = (500, 400))
Plots.plot!(
    plt2,
    daily,
    β,
    label = "",
    xlim = [maximum(daily) - 7, maximum(daily)],
    xlabel = "Day of year",
    #ylim = [minimum(β[1:i_week]), maximum(β[1:i_week])],
    ylabel = "Moisture stress factor",
)

Plots.plot(plt1, plt2, layout = (2, 1))
Plots.savefig(joinpath(savedir, "moisture_stress.png"))


g_stomata =
    [parent(sv.saveval[k].canopy.conductance.gs)[1] for k in 1:length(sol.t)]
plt1 = Plots.plot(size = (500, 400))
Plots.plot!(
    plt1,
    daily,
    g_stomata,
    label = "Model",
    xlim = [minimum(daily), maximum(daily)],
)
i_week = Int(round((7 * 24 * 3600) / (sol.t[2] - sol.t[1])))
plt2 = Plots.plot(size = (500, 400))
Plots.plot!(
    plt2,
    daily,
    g_stomata,
    label = "",
    xlim = [maximum(daily) - 7, maximum(daily)],
    xlabel = "Day of year",
    ylim = [minimum(g_stomata[1:i_week]), maximum(g_stomata[1:i_week])],
    ylabel = "Stomatal conductance",
)

Plots.plot(plt1, plt2, layout = (2, 1))
Plots.savefig(joinpath(savedir, "stomatal_conductance.png"))

plt1 = Plots.plot(size = (500, 400))
Plots.plot!(
    plt1,
    daily,
    [parent(sol.u[k].soil.ϑ_l)[end] for k in 1:1:length(sol.t)],
    label = "10cm",
    xtickfontsize = 5,
    ytickfontsize = 5,
    xlim = [minimum(daily), maximum(daily)],
    ylim = [0.05, 0.55],
    xlabel = "Days",
    ylabel = "SWC [m/m]",
    color = "blue",
)

#plot!(
#    plt1,
#    daily,
#    [parent(sol.u[k].soil.θ_i)[end] for k in 1:1:length(sol.t)],
#    color = "cyan",
#    label = "Ice, 10cm",
#)

plot!(
    plt1,
    daily,
    [parent(sol.u[k].soil.ϑ_l)[end - 2] for k in 1:1:length(sol.t)],
    label = "30cm",
    color = "purple",
)
#plot!(
#    plt1,
#    daily,
#    [parent(sol.u[k].soil.θ_i)[end - 2] for k in 1:1:length(sol.t)],
#    color = "magenta",
#    label = "Ice, 30cm",
#)

Plots.plot!(plt1, seconds ./ 3600 ./ 24, SWC, label = "Data")
plt2 = Plots.plot(
    seconds ./ 3600 ./ 24,
    P .* (-1e3 * 24 * 3600),
    label = "Data",
    ylabel = "Precipitation [mm/day]",
    xlim = [minimum(daily), maximum(daily)],
    ylim = [-200, 0],
    size = (500, 400),
)
Plots.plot(plt2, plt1, layout = (2, 1))
Plots.savefig(joinpath(savedir, "soil_water_content.png"))


SHF = [parent(sv.saveval[k].soil_shf)[1] for k in 1:length(sol.t)]
plt1 = Plots.plot(size = (500, 700))
Plots.plot!(
    plt1,
    seconds ./ 3600 ./ 24,
    FT.(H),
    label = "Data",
    margins = 10Plots.mm,
)
Plots.plot!(
    plt1,
    seconds ./ 3600 ./ 24,
    FT.(H_CORR),
    label = "Data Corr",
    margins = 10Plots.mm,
)
Plots.plot!(
    plt1,
    daily,
    SHF,
    label = "Model",
    xlim = [maximum(daily) - 7, maximum(daily)],
    xlabel = "days",
    ylabel = "Soil Sensible Heat Flux [W/m^2]",
)
Plots.savefig(joinpath(savedir, "shf_week.png"))


dt_model = sol.t[2] - sol.t[1]
dt_data = seconds[2] - seconds[1]
# Find which index in the data our simulation starts at:
idx = argmin(abs.(seconds .- sol.t[1]))
Plots.plot(
    seconds ./ 24 ./ 3600,
    cumsum(measured_T[:]) * dt_data,
    label = "Data ET",
)
Plots.plot!(
    seconds ./ 24 ./ 3600,
    cumsum(P[:]) * dt_data * (1e3 * 24 * 3600),
    label = "Data P",
)
Plots.plot!(
    daily,
    cumsum(T .+ E) * dt_model .+ cumsum(measured_T[:])[idx] * dt_data,
    label = "Model ET",
)

Plots.plot!(ylabel = "∫ Water fluxes dt", xlabel = "Days", margins = 10Plots.mm)
Plots.savefig(joinpath(savedir, "cumul_p_et.png"))


# Leaf water potential data from Pallardy et al (2018)
# Predawn Leaf Water Potential of Oak-Hickory Forest at Missouri Ozark (MOFLUX) Site: 2004-2020
# https://doi.org/10.3334/CDIAC/ORNLSFA.004
lwp_filename = "MOFLUX_PredawnLeafWaterPotential_2020_20210125.csv"
lwp_artifact = ArtifactFile(
    url = "https://caltech.box.com/shared/static/d2nbhezw1q99vslnh5qfwnqrnp3p4edo.csv",
    filename = lwp_filename,
)
lwp_dataset = ArtifactWrapper(
    @__DIR__,
    "lwp_pallardy_etal2018",
    ArtifactFile[lwp_artifact],
);

lwp_path = joinpath(get_data_folder(lwp_dataset), lwp_filename)
lwp_data = readdlm(lwp_path, ',', skipstart = 1)
# We are using 2005 data in this test, so restrict to this year
YEAR = lwp_data[:, 1]
DOY = lwp_data[YEAR .== 2005, 2]
# Our t0 = Dec 31, midnight, 2005. Predawn = guess of 0600 hours
seconds_since_t0 = FT.(DOY) * 24 .* 3600 .+ (6 * 3600)
lwp_measured = lwp_data[YEAR .== 2005, 7] .* 1e6 # MPa to Pa


root_stem_flux = [
    sum(sv.saveval[k].root_extraction) .* (1e3 * 3600 * 24) for
    k in 1:length(sol.t)
]

stem_leaf_flux = [
    parent(sv.saveval[k].canopy.hydraulics.fa)[1] .* (1e3 * 3600 * 24) for
    k in 1:length(sol.t)
]
leaf_air_flux = [
    parent(sv.saveval[k].canopy.hydraulics.fa)[2] .* (1e3 * 3600 * 24) for
    k in 1:length(sol.t)
]
lwp = [
    parent(sv.saveval[k].canopy.hydraulics.ψ)[2] * 9800 for k in 1:length(sol.t)
]
swp = [
    parent(sv.saveval[k].canopy.hydraulics.ψ)[1] * 9800 for k in 1:length(sol.t)
]
ψ_soil = [
    sum(
        parent(sv.saveval[k].soil.ψ) .*
        root_distribution.(parent(cds.subsurface.z)),
    ) / sum(root_distribution.(parent(cds.subsurface.z))) * 9800 for
    k in 1:length(sol.t)
]

plt1 = Plots.plot(size = (500, 400))
Plots.plot!(
    plt1,
    daily,
    lwp,
    label = "Model, Leaf",
    xlim = [minimum(daily), maximum(daily)],
)
Plots.plot!(plt1, daily, swp, label = "Model, Stem")
Plots.plot!(plt1, daily, ψ_soil, label = "Model, Mean soil")
Plots.plot!(
    plt1,
    seconds_since_t0 ./ 24 ./ 3600,
    lwp_measured,
    label = "Data; all species",
    legend = :bottomleft,
)

plt2 = Plots.plot(size = (500, 400))
Plots.plot!(
    plt2,
    daily,
    lwp,
    label = "",
    xlim = [maximum(daily) - 7, maximum(daily)],
    xlabel = "Day of year",
    ylim = [-5e6, 0],
    ylabel = "Water potentials",
)

Plots.plot!(plt2, daily, swp, label = "")
Plots.plot!(plt2, daily, ψ_soil, label = "")
Plots.plot!(plt2, seconds_since_t0 ./ 24 ./ 3600, lwp_measured, label = "")

Plots.plot(plt1, plt2, layout = (2, 1))
Plots.savefig(joinpath(savedir, "plant_hydraulics.png"))

plt2 = Plots.plot(
    daily,
    leaf_air_flux,
    label = "Leaf-air flux",
    xlim = [minimum(daily), maximum(daily)],
    size = (500, 400),
)
Plots.plot!(plt2, daily, stem_leaf_flux, label = "Stem-leaf flux")
Plots.plot!(plt2, daily, root_stem_flux, label = "Soil-root-stem flux")

plt1 = Plots.plot(
    daily,
    leaf_air_flux,
    label = "",
    xlabel = "Day of year",
    ylabel = "Within plant fluxes[mm/day]",
    xlim = [maximum(daily) - 7, maximum(daily)],
    size = (500, 400),
)
Plots.plot!(plt1, daily, stem_leaf_flux, label = "")
Plots.plot!(plt1, daily, root_stem_flux, label = "")
Plots.plot(plt2, plt1, layout = (2, 1))
Plots.savefig(joinpath(savedir, "water_fluxes.png"))

soil_T_sfc = [parent(sv.saveval[k].soil.T)[end] for k in 1:length(sol.t)]
soil_T_20 = [parent(sv.saveval[k].soil.T)[end - 1] for k in 1:length(sol.t)]
soil_T_30 = [parent(sv.saveval[k].soil.T)[end - 2] for k in 1:length(sol.t)]

plt2 = Plots.plot()
Plots.plot!(
    plt2,
    daily,
    soil_T_sfc,
    color = "blue",
    label = "",
    xtickfontsize = 5,
    ytickfontsize = 5,
    xlim = [maximum(daily) - 7, maximum(daily)],
    ylim = [260, 320],
    xlabel = "Day of year",
    ylabel = "Soil Temperature [K]",
)

plot!(plt2, daily, soil_T_30, color = "purple", label = "")
Plots.plot!(plt2, seconds ./ 3600 ./ 24, TA, label = "", color = "green")
LW_OUT[LW_OUT .< 0] .= 0
ϵ_canopy = 0.95
Plots.plot!(
    plt2,
    seconds ./ 3600 ./ 24,
    (LW_OUT ./ (ϵ_canopy * 5.67e-8)) .^ 0.25,
    label = "",
    color = "orange",
)

Plots.plot!(plt2, seconds ./ 3600 ./ 24, TS, label = "", color = "red")



plt3 = Plots.plot()
Plots.plot!(
    plt3,
    daily,
    soil_T_sfc,
    color = "blue",
    label = "T_soil, 10cm",
    xtickfontsize = 5,
    ytickfontsize = 5,
    xlim = [minimum(daily), maximum(daily)],
    ylim = [260, 320],
)

plot!(plt3, daily, soil_T_30, color = "purple", label = "T_soil, 30cm")
Plots.plot!(
    plt3,
    seconds ./ 3600 ./ 24,
    TA,
    label = "T_air (data)",
    color = "green",
)
LW_OUT[LW_OUT .< 0] .= 0
ϵ_canopy = 0.95
Plots.plot!(
    plt3,
    seconds ./ 3600 ./ 24,
    (LW_OUT ./ (ϵ_canopy * 5.67e-8)) .^ 0.25,
    label = "T_sfc (data)",
    color = "orange",
)

Plots.plot!(
    plt3,
    seconds ./ 3600 ./ 24,
    TS,
    label = "T_soil (data)",
    color = "red",
    legend = :bottomright,
)
Plots.plot!(plt3, legend = :bottomleft)
Plots.plot(plt3, plt2, layout = (2, 1))
Plots.savefig(joinpath(savedir, "soil_temperature.png"))

rm(joinpath(savedir, "Artifacts.toml"))