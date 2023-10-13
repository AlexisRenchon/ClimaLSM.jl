# Uses data and computes fluxes using a CLM-style model
using ClimaCore
import CLIMAParameters as CP
using Plots
using Statistics
using Dates
using Insolation
using StatsBase
using StaticArrays
using Thermodynamics
using SurfaceFluxes
using NLsolve

using ClimaLSM
using ClimaLSM.Domains: Column
using ClimaLSM.Soil
using ClimaLSM.Canopy
using ClimaLSM.Canopy.PlantHydraulics
import ClimaLSM
import ClimaLSM.Parameters as LSMP
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))
const FT = Float64
earth_param_set = create_lsm_parameters(FT)
climalsm_dir = pkgdir(ClimaLSM)
include(
    joinpath(
        climalsm_dir,
        "experiments/integrated/ozark/ozark_met_drivers_FLUXNET.jl",
    ),
)
include(
    joinpath(climalsm_dir, "experiments/integrated/ozark/ozark_domain.jl"),
)
include(
    joinpath(climalsm_dir, "experiments/integrated/ozark/ozark_parameters.jl"),
)
include(joinpath(climalsm_dir, "experiments/integrated/flux_helper_functions.jl"))

# Met/Soil/LAI data - get timeseries for dates of interest
P_air = atmos.P.(seconds)
q_air = atmos.q.(seconds)
T_air = atmos.T.(seconds)
u_air = atmos.u.(seconds)
h = atmos_h
LW_d = radiation.LW_d.(seconds)
SW_d = radiation.SW_d.(seconds)
T_soil = TS
ϑ_soil = SWC
LAI = LAIfunction.(seconds)
ref_time = atmos.ref_time
θs = radiation.θs.(seconds, Ref(radiation.orbital_data), ref_time)
c_co2_air = atmos.c_co2.(seconds)


# Constants and model parameters:
_ρ_liq = LSMP.ρ_cloud_liq(earth_param_set)
thermo_params = LSMP.thermodynamic_parameters(earth_param_set)
_LH_v0 = LSMP.LH_v0(earth_param_set)
surface_flux_params = LSMP.surface_fluxes_parameters(earth_param_set)
c = FT(LSMP.light_speed(earth_param_set))
planck_h = FT(LSMP.planck_constant(earth_param_set))
N_a = FT(LSMP.avogadro_constant(earth_param_set))
grav = FT(LSMP.grav(earth_param_set))
R = FT(LSMP.gas_constant(earth_param_set))
M_w = LSMP.molar_mass_water(earth_param_set)
hcm = ClimaLSM.Soil.vanGenuchten(; α = soil_vg_α, n = soil_vg_n)
# soil params
soil_params = Soil.EnergyHydrologyParameters{FT}(;
    κ_dry = κ_dry,
    κ_sat_frozen = κ_sat_frozen,
    κ_sat_unfrozen = κ_sat_unfrozen,
    ρc_ds = ρc_ds,
    ν = soil_ν,
    ν_ss_om = ν_ss_om,
    ν_ss_quartz = ν_ss_quartz,
    ν_ss_gravel = ν_ss_gravel,
    hydrology_cm = hcm,
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


# Canopy Parameterizations - Only radiative_transfer, Photosynthesis, and stomatal conductance needed here.
radiative_transfer = Canopy.TwoStreamModel{FT}(TwoStreamParameters{FT}(;Ω = Ω,
                                                                       ld = ld,
                                                                       α_PAR_leaf = α_PAR_leaf,
                                                                       λ_γ_PAR = λ_γ_PAR,
                                                                       λ_γ_NIR = λ_γ_NIR,
                                                                       τ_PAR_leaf = τ_PAR_leaf,
                                                                       α_NIR_leaf = α_NIR_leaf,
                                                                       τ_NIR_leaf = τ_NIR_leaf,
                                                                       n_layers = n_layers,
                                                                       ϵ_canopy = ϵ_canopy,
                                                                       )
                                               )
conductance = Canopy.MedlynConductanceModel{FT}(MedlynConductanceParameters{FT}(;
                                                                                g1 = g1,
                                                                                Drel = Drel,
                                                                                g0 = g0,
                                                                                )
                                                )

photosynthesis = Canopy.FarquharModel{FT}(FarquharParameters{FT}(Canopy.C3();
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

energy_per_photon_PAR = planck_h * c / λ_γ_PAR
energy_per_photon_NIR = planck_h * c / λ_γ_NIR

# Compute other diagnostic variables as timeseries:
K = @. extinction_coeff(ld, θs)
PAR = ClimaLSM.Canopy.compute_PAR.(Ref(radiative_transfer), Ref(radiation), seconds)
NIR = ClimaLSM.Canopy.compute_NIR.(Ref(radiative_transfer), Ref(radiation), seconds)
e_sat =
    Thermodynamics.saturation_vapor_pressure.(
        Ref(thermo_params),
        T_air,
        Ref(Thermodynamics.Liquid()),
    )
e =
    Thermodynamics.partial_pressure_vapor.(
        Ref(thermo_params),
        P_air,
        PhasePartition.(q_air),
    )
rel_hum = e ./ e_sat
DOY = @. Dates.dayofyear(ref_time + Dates.Second(floor(Int64, seconds)))
frac_diff = @. ClimaLSM.Canopy.diffuse_fraction(DOY, T_air, PAR + NIR, rel_hum, θs)
# Compute the absorbed, reflected, and transmitted PAR and NIR. All positive by defn, in moles of photons/m^2/s
absorbance_tuple = ClimaLSM.Canopy.compute_absorbances.(
    Ref(radiative_transfer),
    PAR / (energy_per_photon_PAR * N_a),
    NIR / (energy_per_photon_NIR * N_a),
    LAI,
    K,
    soil_α_PAR,
    soil_α_NIR,
    θs,
    frac_diff,
)
# It's returned as a Vector of tuples, so unpack that
APAR = [absorbance_tuple[k].par[1].abs for k in 1:length(absorbance_tuple)]
ANIR = [absorbance_tuple[k].nir[1].abs  for k in 1:length(absorbance_tuple)]
TPAR = [absorbance_tuple[k].par[1].trans  for k in 1:length(absorbance_tuple)]
TNIR = [absorbance_tuple[k].nir[1].trans  for k in 1:length(absorbance_tuple)]
RPAR = [absorbance_tuple[k].par[1].refl for k in 1:length(absorbance_tuple)]
RNIR = [absorbance_tuple[k].nir[1].refl for k in 1:length(absorbance_tuple)]

# Compute outgoing SW; net SW of the canopy, and net SW_soil in W/m^2. These are all positive by definition.
SW_out_model = energy_per_photon_NIR * N_a .* RNIR .+ energy_per_photon_PAR * N_a .* RPAR
SW_canopy = energy_per_photon_NIR * N_a .* ANIR .+ energy_per_photon_PAR * N_a .* APAR
SW_soil = energy_per_photon_NIR * N_a .* TNIR .* (1-soil_α_NIR) .+ energy_per_photon_PAR * N_a .* TPAR .* (1-soil_α_PAR)

# Here we make some approximations. This shouldnt be a big deal, though
# because it affects GPP and SW absorbed, and those are two metrics where
# we do fairly well.
# Stomatal conductance. This is inconsistent with T_canopy_airspace, below.
medlyn_factor =
    medlyn_term.(g1, T_air, P_air, q_air, Ref(thermo_params))
# Treat as fixed.
β = 0.96
# Evaluate at T_air. This is inconsistent with T_canopy, below, but 
Rd = @. ClimaLSM.Canopy.dark_respiration(Vcmax25, β, f, ΔHRd, T_air, To, R)
An =
    ClimaLSM.Canopy.compute_photosynthesis.(
        Ref(photosynthesis),
        T_air,
        medlyn_factor,
        APAR,
        c_co2_air,
        β,
        R,
        Rd,
    )

GPP = @. ClimaLSM.Canopy.compute_GPP(An, K, LAI, Ω)
gs = @. ClimaLSM.Canopy.medlyn_conductance(g0, Drel, medlyn_factor, An, c_co2_air)
g_canopy = @. ClimaLSM.Canopy.upscale_leaf_conductance(gs, LAI, T_air, R, P_air)
r_canopy = 1 ./ g_canopy
ψ_soil = @. ClimaLSM.Soil.matric_potential(hcm, (SWC - θ_r)/(0.52-θ_r))
r_soil = ClimaLSM.Soil.soil_resistance.(SWC, SWC, 0.0, soil_params)
ϵ_soil = 0.96
# Needed to compute resistances, below.
AI = SAI .+ LAI
W = @. exp(-AI) # CLM Equation 2.5.119


ts_in = @. Thermodynamics.PhaseEquil_pTq(thermo_params, P_air, T_air, q_air)
ρ_air = @. Thermodynamics.air_density(thermo_params, ts_in)
cp_m = @. Thermodynamics.cp_m(thermo_params, ts_in)
Rm_int = @. Thermodynamics.gas_constant_air(thermo_params, ts_in)
# Extrapolate to the surface to get the density of the air right above the ground.
ρ_soil = @. ρ_air *
        (T_soil / T_air)^(Thermodynamics.cv_m(thermo_params, ts_in) / Rm_int)
q_sat_soil = @. Thermodynamics.q_vap_saturation_generic(
    thermo_params,
    T_soil,
    ρ_soil,
    Thermodynamics.Liquid(),
)
q_soil= @. (q_sat_soil * exp(grav * ψ_soil * M_w / (R * T_soil))) # Bonan's book



# Variables we will save
T_airspace_model = []
q_airspace_model = []
T_canopy_model = []
shf = []
lhf = []
soil_shf = []
soil_lhf = []
canopy_lhf = []
canopy_shf = []
ustar_model = []
LW_out_model = []
LW_canopy = []
LW_soil = []
L_MO = []


# step through the data and compute surface fluxes according to the model
steps = 120*48:1:180*48
# Errors on near step 7314
for step in 120*48:1:7313
    @info(step)
    # Tairspace, qairspace, Tcanopy, ustar, shf, lhf, LMO. Only the first three are optimized, so we set F[4:7] = 0.
    initial_guess = [0.5 .*(T_air[step] .+ T_soil[step]), 0.5 .*(q_air[step] .+ q_soil[step]), 0.5 .*(T_air[step] .+ T_soil[step]), u_air[step], 0.0, 0.0, 0.0]
    function flux_equality(F, x)
        # Compute the surface fluxes between the airspace and the atmosphere. CLM 2.5.100 and CLM2.5.87. See discussion right after 2.5.92.
        # The surface corresponds to the canopy airspace, with values at height z_0+d
        T_airspace = x[1]
        q_airspace = x[2]
        T_canopy = x[3]
        
        ts_sfc = Thermodynamics.PhaseEquil_ρTq(thermo_params, ρ_air[step], T_airspace, q_airspace)
        airspace_atmos_fluxes =  atmos_airspace_fluxes(ts_sfc, ts_in[step], h, d_sfc, u_air[step],z_0m, z_0b, surface_flux_params)
        ustar = airspace_atmos_fluxes.ustar
        soil_fluxes = airspace_soil_fluxes(ustar, T_airspace, q_airspace, W[step], T_soil[step], q_soil[step], r_soil[step], ρ_air[step], cp_m[step], h-d_sfc)
        canopy_fluxes = airspace_canopy_fluxes(ustar, T_airspace, q_airspace, T_canopy, ρ_air[step], cp_m[step], LAI[step], AI[step], r_canopy[step])
        lw_fluxes = land_lw_fluxes(LW_d[step], T_canopy, T_soil[step], ϵ_soil, ϵ_canopy, _σ)
        # Functions to find the roots of:
        F[1] = canopy_fluxes.shf + soil_fluxes.shf - airspace_atmos_fluxes.shf
        F[2] = (canopy_fluxes.transpiration + soil_fluxes.evaporation - airspace_atmos_fluxes.evaporation) * _LH_v0
        F[3] = (-lw_fluxes.LW_canopy - SW_canopy[step] + canopy_fluxes.transpiration * _LH_v0 + canopy_fluxes.shf)
        F[4] = 0.0
        F[5] = 0.0
        F[6] = 0.0
        F[7] = 0.0

        # Store latest values from this iteration.
        x[4] = ustar
        x[5] = airspace_atmos_fluxes.shf
        x[6] = airspace_atmos_fluxes.lhf
        x[7] = airspace_atmos_fluxes.L_MO
    end
    soln = nlsolve(
        flux_equality,
        initial_guess,
        ftol = 1e-8,
        iterations = 100
    )
    T_airspace = soln.zero[1]
    q_airspace = soln.zero[2]
    T_canopy = soln.zero[3]
    ustar= soln.zero[4]

    soil_fluxes = airspace_soil_fluxes(ustar, T_airspace, q_airspace, W[step], T_soil[step], q_soil[step], r_soil[step], ρ_air[step], cp_m[step], h-d_sfc)
    canopy_fluxes = airspace_canopy_fluxes(ustar, T_airspace, q_airspace, T_canopy, ρ_air[step], cp_m[step], LAI[step], AI[step], r_canopy[step])
    lw_fluxes = land_lw_fluxes(LW_d[step], T_canopy, T_soil[step], ϵ_soil, ϵ_canopy, _σ)
                       

    # Save output
    push!(T_airspace_model, T_airspace)
    push!(q_airspace_model, q_airspace)
    push!(ustar_model, ustar)
    push!(shf, soln.zero[5])
    push!(lhf, soln.zero[6])
    push!(L_MO, soln.zero[7])

    push!(soil_shf, soil_fluxes.shf)
    push!(soil_lhf, soil_fluxes.evaporation * _LH_v0)
    push!(canopy_lhf, canopy_fluxes.transpiration * _LH_v0)
    push!(canopy_shf,  canopy_fluxes.shf)
    push!(LW_canopy, lw_fluxes.LW_canopy)
    push!(T_canopy_model, T_canopy)
    push!(LW_out_model,lw_fluxes.LW_out)
    push!(LW_soil, lw_fluxes.LW_soil)
end


# Plotting
idx_end = length(T_airspace_model)
days = seconds[steps[1:idx_end]] ./ 24 ./ 3600

# Weird spikes but not wrong immediately. due to MO solve?
plt_temp = Plots.plot(days, T_airspace_model, label = "Airspace", title = "Temperature")
Plots.plot!(plt_temp, days, T_air[steps[1:idx_end]], label = "Atmos")
Plots.plot!(plt_temp, days, T_soil[steps[1:idx_end]], label = "Soil")
Plots.plot!(plt_temp, days, T_canopy_model, label = "Canopy")

# Weird spikes but not wrong immediately
plt_lw = Plots.plot(days, LW_out_model, label = "Model", title = "LW out")
Plots.plot!(plt_lw, days, LW_OUT[steps[1:idx_end]], label = "Data")

plt_sw = Plots.plot(days, SW_out_model[steps[1:idx_end]], label = "Model", title = "SW out")
Plots.plot!(plt_sw, days, SW_OUT[steps[1:idx_end]], label = "Data")

# same spikes
plot_q = Plots.plot(days, q_airspace_model, label = "Airspace", title = "q")
Plots.plot!(plot_q, days, q_air[steps[1:idx_end]], label = "Atmos")
Plots.plot!(plot_q, days, q_soil[steps[1:idx_end]], label = "Soil")

# G still looks bad
# For us, G is positive = Loss of heat from the ground -> add a sign to data
plot_G = Plots.plot(days, soil_shf, label = "Soil SHF", title = "Ground heat flux")
Plots.plot!(plot_G, days, -1 .* LW_soil .- SW_soil[steps[1:idx_end]] .+ soil_shf .+ soil_lhf, label = "Net soil heat flux")
Plots.plot!(plot_G, days, -1 .* G[steps[1:idx_end]], label = "Data")

# Looks wrong all around
plot_shf = Plots.plot(days, soil_shf, label  = "Soil-Airspace", title = "SHF")
Plots.plot!(plot_shf, days, canopy_shf, label  = "Canopy-Airspace")
Plots.plot!(plot_shf, days, shf, label  = "Airspace-Atmos")
Plots.plot!(plot_shf, days, H[steps[1:idx_end]], label = "data",)

# Looks pretty good
plot_lhf = Plots.plot(days, soil_lhf, label  = "Soil-Airspace", title = "LHF")
Plots.plot!(plot_lhf, days, canopy_lhf, label  = "Canopy-Airspace")
Plots.plot!(plot_lhf, days, lhf, label  = "Airspace-Atmos")
Plots.plot!(plot_lhf, days, LE[steps[1:idx_end]], label = "data")


# Spikes when ustar -> 0?
plot_ustar = Plots.plot(days, ustar_model, title = "Ustar")


# Energy balance at the flux tower site:

RminusG = @.(LW_d-LW_OUT) + (SW_d - SW_OUT)
HplusL = LE .+ H
Plots.scatter(RminusG, HplusL)
