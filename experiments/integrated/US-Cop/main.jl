# THIS SITE IS A GRASSLAND 
cur_dir = @__DIR__
#cur_dir="/Users/mitraasadollahi/Projects/CliMA/ClimaLSM.jl/experiments/integrated/US-Cop"
include(
    joinpath(
        cur_dir,
        "load_packages.jl",
    ),
)
function increment_datetime_by_days(base_datetime::DateTime, days::Float64)::DateTime
    whole_days = floor(Int, days)
    remaining = days - whole_days
    hours = floor(Int, remaining * 24)
    remaining = (remaining * 24) - hours
    minutes = floor(Int, remaining * 60)
    remaining = (remaining * 60) - minutes
    seconds = round(Int, remaining * 60) # rounding to nearest second

    return base_datetime + Day(whole_days) + Hour(hours) + Minute(minutes) + Second(seconds)
end

function diurnal_avg(series)
    num_days = Int64(N_days - N_spinup_days)
    daily_points = Int64(length(series) ÷ num_days) # Num datapoints per day
    daily_data = [
        series[i:1:(i+daily_points-1)] for
        i in 1:daily_points:(daily_points*num_days)
    ]
    daily_avgs =
        [mean([daily_data[i][j] for i in 1:num_days]) for j in 1:daily_points]
    return daily_avgs
end

global site_name = "US-Cop";
substring_to_remove = "/experiments/integrated/$site_name"
climalsm_dir = replace(cur_dir, substring_to_remove => "")
metadata_dir = joinpath(replace(cur_dir, "$site_name" => ""), "input")
global main_dir = replace(climalsm_dir, "ClimaLSM.jl" => "")
include(joinpath(climalsm_dir, "parameters", "create_parameters.jl"))
global start_year = 1#minimum(daily)
global end_year = 1#maximum(daily)
global start_date = 1 + (start_year - 1) * 365#maximum(daily)
global end_date = 350 + (end_year - 1) * 365#maximum(daily)
savedir = joinpath(climalsm_dir, "experiments/integrated/" * site_name * "/results/")
global savedir_input = joinpath(climalsm_dir, "experiments/integrated/" * site_name * "/input_drivers/")
global plot_input_var=false
if !isdir(savedir_input)
    mkdir(savedir_input)
    global plot_input_var=true
end
curdir = joinpath(climalsm_dir, "experiments/integrated/" * site_name * "/")
const FT = Float64
earth_param_set = create_lsm_parameters(FT)
df = CSV.File(joinpath(metadata_dir, "fluxnet_site_location.csv")) |> DataFrame

# Filter the dataframe for the specific site_name (assuming you have a variable `site_name` with the desired name)
site_row = df[lowercase.(df.Site).==lowercase(site_name), :]
# Extract Longitude and Latitude
global long = FT(site_row[1, :Longitude])
global lat = FT(site_row[1, :Latitude])
global spinup = 0 #allows 1 year of spinup
global tf
global t_spinup
# This reads in the data from the flux tower site and creates
# the atmospheric and radiative driver structs for the model
include(
    joinpath(
        cur_dir,
        "met_drivers_FLUXNET.jl",
    ),
)
include(
    joinpath(cur_dir, "parameters.jl"),
)
include(joinpath(cur_dir, "domain.jl"))
include(
    joinpath(cur_dir, "simulation.jl"),
)

# Now we set up the model. For the soil model, we pick
# a model type and model args:
soil_domain = land_domain
soil_ps = Soil.EnergyHydrologyParameters{FT}(;
    κ_dry=κ_dry,
    κ_sat_frozen=κ_sat_frozen,
    κ_sat_unfrozen=κ_sat_unfrozen,
    ρc_ds=ρc_ds,
    ν=soil_ν,
    ν_ss_om=ν_ss_om,
    ν_ss_quartz=ν_ss_quartz,
    ν_ss_gravel=ν_ss_gravel,
    hydrology_cm=vanGenuchten(; α=soil_vg_α, n=soil_vg_n),
    K_sat=soil_K_sat,
    S_s=soil_S_s,
    θ_r=θ_r,
    earth_param_set=earth_param_set,
    z_0m=z_0m_soil,
    z_0b=z_0b_soil,
    emissivity=soil_ϵ,
    PAR_albedo=soil_α_PAR,
    NIR_albedo=soil_α_NIR
);

soil_args = (domain=soil_domain, parameters=soil_ps)
soil_model_type = Soil.EnergyHydrology{FT}

# Now we set up the canopy model, which we set up by component:
# Component Types
canopy_component_types = (;
    autotrophic_respiration=Canopy.AutotrophicRespirationModel{FT},
    radiative_transfer=Canopy.TwoStreamModel{FT},
    photosynthesis=Canopy.FarquharModel{FT},
    conductance=Canopy.MedlynConductanceModel{FT},
    hydraulics=Canopy.PlantHydraulicsModel{FT}
)
# Individual Component arguments
# Set up autotrophic respiration
autotrophic_respiration_args = (;
    parameters=AutotrophicRespirationParameters{FT}(;
        ne=ne,
        ηsl=ηsl,
        σl=σl,
        μr=μr,
        μs=μs,
        f1=f1,
        f2=f2
    )
)
# Set up radiative transfer
radiative_transfer_args = (;
    parameters=TwoStreamParameters{FT}(;
        Ω=Ω,
        ld=ld,
        α_PAR_leaf=α_PAR_leaf,
        λ_γ_PAR=λ_γ_PAR,
        λ_γ_NIR=λ_γ_NIR,
        τ_PAR_leaf=τ_PAR_leaf,
        α_NIR_leaf=α_NIR_leaf,
        τ_NIR_leaf=τ_NIR_leaf,
        n_layers=n_layers,
        ϵ_canopy=ϵ_canopy
    )
)
# Set up conductance
conductance_args = (;
    parameters=MedlynConductanceParameters{FT}(;
        g1=g1,
        Drel=Drel,
        g0=g0
    )
)
# Set up photosynthesis
photosynthesis_args = (;
    parameters=FarquharParameters{FT}(
        Canopy.C3();
        oi=oi,
        ϕ=ϕ,
        θj=θj,
        f=f,
        sc=sc,
        pc=pc,
        Vcmax25=Vcmax25,
        Γstar25=Γstar25,
        Kc25=Kc25,
        Ko25=Ko25,
        To=To,
        ΔHkc=ΔHkc,
        ΔHko=ΔHko,
        ΔHVcmax=ΔHVcmax,
        ΔHΓstar=ΔHΓstar,
        ΔHJmax=ΔHJmax,
        ΔHRd=ΔHRd
    )
)
# Set up plant hydraulics
ai_parameterization = PrescribedSiteAreaIndex{FT}(LAIfunction, SAI, RAI)

function root_distribution(z::T; rooting_depth=rooting_depth) where {T}
    return T(1.0 / rooting_depth) * exp(z / T(rooting_depth)) # 1/m
end

plant_hydraulics_ps = PlantHydraulics.PlantHydraulicsParameters(;
    ai_parameterization=ai_parameterization,
    ν=plant_ν,
    S_s=plant_S_s,
    root_distribution=root_distribution,
    conductivity_model=conductivity_model,
    retention_model=retention_model
)
plant_hydraulics_args = (
    parameters=plant_hydraulics_ps,
    n_stem=n_stem,
    n_leaf=n_leaf,
    compartment_midpoints=compartment_midpoints,
    compartment_surfaces=compartment_surfaces,
)

# Canopy component args
canopy_component_args = (;
    autotrophic_respiration=autotrophic_respiration_args,
    radiative_transfer=radiative_transfer_args,
    photosynthesis=photosynthesis_args,
    conductance=conductance_args,
    hydraulics=plant_hydraulics_args
)
# Other info needed
shared_params = SharedCanopyParameters{FT,typeof(earth_param_set)}(
    z0_m,
    z0_b,
    earth_param_set,
)

canopy_model_args = (; parameters=shared_params, domain=canopy_domain)

# Integrated plant hydraulics and soil model
land_input = (atmos=atmos, radiation=radiation)
land = SoilCanopyModel{FT}(;
    land_args=land_input,
    soil_model_type=soil_model_type,
    soil_args=soil_args,
    canopy_component_types=canopy_component_types,
    canopy_component_args=canopy_component_args,
    canopy_model_args=canopy_model_args
)
Y, p, cds = initialize(land)
exp_tendency! = make_exp_tendency(land)

#Initial conditions
Y.soil.ϑ_l = SWC[1+Int(round(t0 / DATA_DT))] # Get soil water content at t0
println("initial saturation is $(Y.soil.ϑ_l)")
# Both data and simulation are reference to 2005-01-01-00 (LOCAL)
# or 2005-01-01-06 (UTC)
Y.soil.θ_i = FT(0.0)
T_0 = TS[1+Int(round(t0 / DATA_DT))] # Get soil temperature at t0
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
    Y.canopy.hydraulics.ϑ_l.:($i) .=
        augmented_liquid_fraction.(plant_ν, S_l_ini[i])
end

set_initial_aux_state! = make_set_initial_aux_state(land)
set_initial_aux_state!(p, Y, t0);

# Simulation
sv = (;
    t=Array{FT}(undef, length(saveat)),
    saveval=Array{NamedTuple}(undef, length(saveat))
)
cb = ClimaLSM.NonInterpSavingCallback(sv, saveat)

prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction((T_exp!)=exp_tendency!),
    Y,
    (t0, tf),
    p,
);
sol = SciMLBase.solve(
    prob,
    ode_algo;
    dt=dt,
    callback=cb,
    adaptive=false,
    saveat=saveat
)

# Plotting
global n
daily = sol.t ./ 3600 ./ 24
# Number of datapoints per day
data_daily_points = Int64(86400 / DATA_DT)
# Number of model points per day
model_daily_points = Int64(86400 / n / dt)
# Scales data indices 0 to 24 in a day
data_daily_indices = range(0, step=DATA_DT / 3600, length=data_daily_points)
model_daily_indices =
    range(0, step=dt * n / 3600, length=model_daily_points)
# Number of data points per each model point (Ratio of data dt to model dt)
data_per_model = Int64(dt * n ÷ DATA_DT)

# This function will be used to compute averages over diurnal cycles. Input a
# data series to average over and the total number of days in the data series, 
# and it will return a vector of the average value at each timestamp in the 
# day averaged over every day in the series.
global N_spinup_days
global N_days
#measured data
measured_ET = LE ./ (LSMP.LH_v0(earth_param_set) * 1000) .* (1e3 * 24 * 3600)
tmpt = sol.t ./ 3600 ./ 24
#create measurement dataframe
global LOCAL_DATETIME
model_data = DataFrame()
model_data[!, "time"] = sol.t ./ 3600 ./ 24
model_data[!, "date"] = increment_datetime_by_days.(Ref(LOCAL_DATETIME[1]), tmpt)
num_rows = length(vec(tmpt))
#Meas = DataFrame(SWC = Symbol[], GPP = Symbol[], ET = Symbol[])
Obs = ["GPP", "SWC", "ET", "T", "E", "SHF"]
for col_name in Obs
    model_data[!, col_name] = zeros(Float64, num_rows)  # Replace `num_rows` with the desired length
end
model_data[!, "GPP"] = [
    parent(sv.saveval[k].canopy.photosynthesis.GPP)[1] for
    k in 1:length(sv.saveval)]
model_data[!, "SWC"] = [parent(sol.u[k].soil.ϑ_l)[end-1] for k in 1:1:length(sol.t)];
model_data[!, "E"] = [parent(sv.saveval[k].soil_evap)[1] for k in 1:length(sol.t)] .* (1e3 * 24 * 3600)
model_data[!, "T"] = [
    parent(sv.saveval[k].canopy.conductance.transpiration)[1] for
    k in 1:length(sol.t)
] .* (1e3 * 24 * 3600)
model_data[!, "ET"] = model_data[!, "E"] .+ model_data[!, "T"]
model_data[!, "SHF"] = [parent(sv.saveval[k].soil_shf)[1] for k in 1:length(sol.t)]

# Write all data to a csv file
CSV.write(joinpath(savedir, "model_output.csv"), model_data)
all_columns = names(model_data)

# Exclude date and time columns
columns_to_plot = setdiff(all_columns, ["date", "time"])
measured_ET = LE ./ (LSMP.LH_v0(earth_param_set) * 1000) .* (1e3 * 24 * 3600)
#=
println(" length of T is $(length(measured_ET))")
println(" length of LOCAL_DATETIME is $(length(LOCAL_DATETIME))")
println(" tf is $tf")
println(" t_spinup is $t_spinup")
println(" DATA_DT is $DATA_DT")
println(" length of ET is $(length(measured_ET))")
println(" length of GPP is $(length(GPP))")
println(" length of LOCAL_DATETIME is $(length((LOCAL_DATETIME[1:24:end])))")
=#
ls=Int64(round((LOCAL_DATETIME[end].-LOCAL_DATETIME[1]) ./ Millisecond(60*60*24*1000)))
sp=Int64(round.(length(LOCAL_DATETIME)/ls))
global id=1:sp:length(LOCAL_DATETIME)
# Plot each column and save the plot
for col in columns_to_plot
    println(col)
    global LOCAL_DATETIME
    global GPP
    global SWC
    global ET
    global id
    prnd = plot(model_data.date, model_data[!, col], xticks=:out, yticks=:out, xlabel="Time", ylabel=col,label="CliMa",dpi=400)
    if col.=="GPP"
        plot!(prnd,LOCAL_DATETIME[1:24:end],GPP[1:24:end],markershape=:star,color="red",label="daily FLUXNET data",dpi=400)
    elseif col.=="ET"
        plot!(prnd,LOCAL_DATETIME[1:24:end],measured_ET[1:24:end],markershape=:star,color="red",label="daily FLUXNET data",dpi=400,ylim=(0,maximum(measured_ET)))
    elseif col.=="SWC"
        plot!(prnd,LOCAL_DATETIME[1:24:end],SWC[1:24:end],markershape=:star,color="red",label="daily FLUXNET data",dpi=400)
    end
    savefig(prnd, joinpath(savedir, "$(col).png"))
end
@info "Saved model output to $(savedir)model_output.csv"
#end

file_path = joinpath(cur_dir, "Artifacts.toml")

if isfile(file_path)
    rm(file_path)
    println("File removed: ", file_path)
else
    println("File does not exist: ", file_path)
end
