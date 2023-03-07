module Canopy
using Thermodynamics
using ClimaLSM
using ClimaCore
using ClimaLSM: AbstractModel, AbstractRadiativeDrivers, AbstractAtmosphericDrivers
import ..Parameters as LSMP

import ClimaLSM:
    name,
    domain,
    prognostic_vars,
    prognostic_types,
    auxiliary_vars,
    auxiliary_types,
    initialize_prognostic,
    initialize_auxiliary,
    make_update_aux,
    make_ode_function
    
using ClimaLSM.Domains: Point, Plane, SphericalSurface
using DocStringExtensions
export SharedCanopyParameters, CanopyModel, DiagnosticTranspiration
include("./PlantHydraulics.jl")
using .PlantHydraulics
import .PlantHydraulics: transpiration, AbstractTranspiration
include("./stomatalconductance.jl")
include("./photosynthesis.jl")
include("./radiation.jl")
include("./canopy_parameterizations.jl")

"""
    SharedCanopyParameters{FT <: AbstractFloat, PSE}

A place to store shared parameters that are required by all canopy components.
$(DocStringExtensions.FIELDS)
"""
struct SharedCanopyParameters{FT <: AbstractFloat, PSE}
    "Leaf Area Index"
    LAI::FT
    "Canopy height"
    h_c::FT
    #"Rate of change of roughness length for momentum with canopy height"
    # ω_m::FT
    "Roughness length for momentum"
    z_0m::FT
    "Roughness length for scalars"
    z_0b::FT
    "Earth param set"
    earth_param_set::PSE
end

struct CanopyModel{FT,RM, PM, SM, PHM, A, R, PS, D} <: AbstractModel{FT}
    radiative_transfer::RM
    photosynthesis::PM
    conductance::SM
    hydraulics::PHM
    atmos::A
    radiation::R
    parameters::PS
    domain::D
end
    
function CanopyModel{FT}(;
                         radiative_transfer::AbstractRadiationModel{FT},
                         photosynthesis::AbstractPhotosynthesisModel{FT},
                         conductance::AbstractStomatalConductanceModel{FT},
                         hydraulics::AbstractPlantHydraulicsModel{FT},
                         atmos::AbstractAtmosphericDrivers{FT},
                         radiation::AbstractRadiativeDrivers{FT},
                         parameters::SharedCanopyParameters{FT, PSE},
                         domain::Union{ClimaLSM.Domains.Point, ClimaLSM.Domains.Plane, ClimaLSM.Domains.SphericalSurface}) where {FT, PSE}
    if hydraulics.domain != domain
        throw(
            AssertionError(
                "The provided canopy model domain must be the same as the hydraulics model domain.",
            ),
        )
    end
    args = (radiative_transfer, photosynthesis, conductance, hydraulics,atmos, radiation, parameters, domain)
    return CanopyModel{FT, typeof.(args)...}(args...)
end

ClimaLSM.name(::CanopyModel) = :canopy
ClimaLSM.domain(::CanopyModel) = :surface
canopy_components(::CanopyModel) = (:hydraulics, :conductance, :photosynthesis, :radiative_transfer)

function prognostic_vars(canopy::CanopyModel)
    components = canopy_components(canopy)
    prognostic_list = map(components) do model
        prognostic_vars(getproperty(canopy, model))
    end
    return NamedTuple{components}(prognostic_list)
end

function prognostic_types(canopy::CanopyModel)
    components = canopy_components(canopy)
    prognostic_list = map(components) do model
        prognostic_types(getproperty(canopy, model))
    end
    return NamedTuple{components}(prognostic_list)
end

function auxiliary_vars(canopy::CanopyModel)
    components = canopy_components(canopy)
    auxiliary_list = map(components) do model
        auxiliary_vars(getproperty(canopy, model))
    end
    return NamedTuple{components}(
        auxiliary_list,
    )
end

function auxiliary_types(canopy::CanopyModel)
    components = canopy_components(canopy)
    auxiliary_list = map(components) do model
        auxiliary_types(getproperty(canopy, model))
    end
    return NamedTuple{components}(
        auxiliary_list,
    )
end


function initialize_prognostic(model::CanopyModel{FT}, coords::ClimaCore.Fields.Field) where {FT}
    components = canopy_components(model)
    Y_state_list = map(components) do (component)
        submodel = getproperty(model, component)
        zero_state =
            map(_ -> zero(FT), coords)
        getproperty(initialize_prognostic(submodel, zero_state), component)
    end
    Y = ClimaCore.Fields.FieldVector(; name(model) => NamedTuple{components}(Y_state_list))
    return Y
end

function initialize_auxiliary(model::CanopyModel{FT}, coords::ClimaCore.Fields.Field) where {FT}
    components = canopy_components(model)
    p_state_list = map(components) do (component)
        submodel = getproperty(model, component)
        zero_state =
            map(_ -> zero(FT), coords)
        getproperty(initialize_auxiliary(submodel, zero_state), component)
    end
    p = ClimaCore.Fields.FieldVector(; name(model) => NamedTuple{components}(p_state_list))
    return p
end

function ClimaLSM.make_update_aux(canopy::CanopyModel{FT,
                                                      <: BeerLambertModel,
                                                      <: FarquharModel,
                                                      <: MedlynConductanceModel,
                                                      <: PlantHydraulicsModel 
                                                      }) where {FT}
    function update_aux!(p, Y, t)
        # unpack params
        earth_param_set = canopy.parameters.earth_param_set
        R = FT(LSMP.gas_constant(earth_param_set))
        c = FT(LSMP.light_speed(earth_param_set))
        h = FT(LSMP.planck_constant(earth_param_set))
        N_a = FT(LSMP.avogadro_constant(earth_param_set))
        thermo_params = canopy.parameters.earth_param_set.thermo_params

        (; Vcmax25, Γstar25, ΔHJmax, ΔHVcmax, ΔHΓstar, f, ΔHRd, To, θj, ϕ, mechanism, sc, ψc, oi, Kc25, Ko25, ΔHkc, ΔHko) = canopy.photosynthesis.parameters 
        (; g1, g0, Drel) = canopy.conductance.parameters
        (; ld, Ω, ρ_leaf) = canopy.radiative_transfer.parameters
        (; LAI) = canopy.parameters
        λ = FT(5e-7) # m (500 nm)
        energy_per_photon = h * c / λ

        top_index = canopy.hydraulics.n_stem + canopy.hydraulics.n_leaf
        c_co2::FT = canopy.atmos.c_co2(t)
        P::FT = canopy.atmos.P(t)
        u::FT = canopy.atmos.u(t)
        T::FT = canopy.atmos.T(t)
        h::FT = canopy.atmos.h
        q::FT = canopy.atmos.q(t)
        SW_d::FT = canopy.radiation.SW_d(t)
        LW_d::FT = canopy.radiation.LW_d(t)  
        θs::FT  = canopy.radiation.θs(t)      
        # atmos_ts = construct_atmos_ts(canopy.atmos, t, thermo_params)
        # compute VPD
        es = Thermodynamics.saturation_vapor_pressure.(Ref(thermo_params), T, Ref(Thermodynamics.Liquid()))
        ea = @.(q * P / (0.622 + 0.378 * q))
        VPD = es .- ea
        PAR = @.(SW_d / (energy_per_photon * N_a) / 2)

        K = extinction_coeff.(ld, θs)
        APAR = plant_absorbed_ppfd.(PAR, ρ_leaf, K, LAI, Ω)   
        β = moisture_stress.(p.canopy.hydraulics.ψ[top_index], sc, ψc) 
        Jmax = max_electron_transport.(Vcmax25, ΔHJmax,T,To,R)
        J = electron_transport.(APAR, Jmax, θj, ϕ)
        Vcmax = compute_Vcmax.(Vcmax25, T, To, R, ΔHVcmax)
        Γstar = co2_compensation.(Γstar25, ΔHΓstar, T, To, R)
        m = medlyn_term.(g1, VPD)
        ci = intercellular_co2.(c_co2, Γstar, m)
        Aj = light_assimilation.(Ref(mechanism), J, ci, Γstar)
        Kc = MM_Kc.(Kc25, ΔHkc, T, To, R)
        Ko = MM_Ko.(Ko25, ΔHko, T, To, R)
        Ac = rubisco_assimilation.(Ref(mechanism), Vcmax, ci, Γstar, Kc, Ko, oi)
        Rd = dark_respiration.(Vcmax25, β, f, ΔHRd, T, To, R)
        
        p.canopy.radiative_transfer.apar .= APAR
        p.canopy.photosynthesis.An .= net_photosynthesis.(Ac, Aj, Rd, β)
        p.canopy.photosynthesis.GPP .= compute_GPP.(p.canopy.photosynthesis.An, K, LAI, Ω)
        p.canopy.conductance.gs .= medlyn_conductance.(g0, Drel, m, p.canopy.photosynthesis.An, c_co2)
        p.canopy.conductance.medlyn_term .= m
        (transpiration, shf, lhf) = canopy_surface_fluxes(canopy.atmos, canopy, Y, p, t)
        p.canopy.hydraulics.fa[top_index] .= transpiration
        # note, confusingly, that the other aux variables for plant hydraulics are updated in the RHS
    end
    return update_aux!
end

function make_ode_function(canopy::CanopyModel)
    components = canopy_components(canopy)
    rhs_function_list = map(x -> make_rhs(getproperty(canopy, x)), components)
    update_aux! = make_update_aux(canopy)
    function ode_function!(dY, Y, p, t)
        update_aux!(p, Y, t)
        for f! in rhs_function_list
            f!(dY, Y, p, t)
        end
    end
    return ode_function!
end

"""
    canopy_surface_fluxes(atmos::PrescribedAtmosphere{FT},
                          model::CanopyModel,
                          Y,
                          p,
                          t::FT) where {FT}

Computes canopy transpiration 
"""
function canopy_surface_fluxes(atmos::PrescribedAtmosphere{FT},
                               model::CanopyModel,
                               Y,
                               p,
                               t::FT) where {FT}
    # in the long run, we should pass r_sfc to surface_fluxes
    # where it would be handle internally.
    # but it doesn't do that, so we need to hack together something after
    # the fact
    # TODO: adjust latent heat flux for stomatal conductance
    # E_potential = g_ae (q_sat(T_sfc) - q_atmos) 
    # T = g_eff (q_sat(T_sfc) - q_atmos)  < E_potential

    base_lhf, shf, base_transpiration, C_h = surface_fluxes(atmos, model, Y, p, t)
    earth_param_set = model.parameters.earth_param_set
    # here is where we adjust evaporation for stomatal conductance = 1/r_sfc
    r_ae = 1/(C_h * abs(atmos.u(t))) # s/m
    ρ_m = FT(LSMP.ρ_cloud_liq(earth_param_set) / LSMP.molar_mass_water(earth_param_set))

    r_sfc = @. 1/(p.canopy.conductance.gs/ρ_m) # should be s/m
    r_eff = r_ae .+ r_sfc
    transpiration = @. base_transpiration*r_ae/r_eff

    # we also need to correct the LHF
    lhf = @. base_lhf *r_ae/r_eff
    return transpiration, shf, lhf
end


"""
    ClimaLSM.surface_temperature(model::CanopyModel, Y, p, t)

a helper function which returns the surface temperature for the canopy 
model, which is stored in the aux state.
"""
function ClimaLSM.surface_temperature(model::CanopyModel, Y, p, t)
    return model.atmos.T(t) 
end

function ClimaLSM.surface_height(model::CanopyModel)
    return model.parameters.h_c 
end

"""
    ClimaLSM.surface_temperature(model::CanopyModel, Y, p)

a helper function which returns the surface specific humidity for the canopy 
model, which is stored in the aux state.
"""
function ClimaLSM.surface_specific_humidity(model::CanopyModel, Y, p, T_sfc, ρ_sfc)
    thermo_params = LSMP.thermodynamic_parameters(model.parameters.earth_param_set)
    return Thermodynamics.q_vap_saturation_generic.(
        Ref(thermo_params),
        T_sfc,
        ρ_sfc,
        Ref(Thermodynamics.Liquid()),
    )
end

"""
    ClimaLSM.surface_temperature(model::CanopyModel, Y, p)
    
a helper function which computes and returns the surface air density for the canopy 
model.
"""
function ClimaLSM.surface_air_density(
    atmos::PrescribedAtmosphere,
    model::CanopyModel,
    Y,
    p,
    t,
    T_sfc,
)
    thermo_params =
        LSMP.thermodynamic_parameters(model.parameters.earth_param_set)
    ts_in = construct_atmos_ts(atmos, t, thermo_params)
    return compute_ρ_sfc.(Ref(thermo_params), Ref(ts_in), T_sfc)
end

"""
    ClimaLSM.surface_temperature(model::CanopyModel, Y, p)

a helper function which computes and returns the surface evaporative scaling
 factor for the canopy model.
"""
function ClimaLSM.surface_evaporative_scaling(model::CanopyModel{FT}, Y, p) where {FT}
    beta = FT(1.0)
    return beta
end



"""
    DiagnosticTranspiration{FT} <: AbstractTranspiration{FT}

A concrete type used for dispatch when computing the transpiration
from the leaves, in the case where transpiration is computed
diagnostically, as a function of prognostic variables and parameters
"""
struct DiagnosticTranspiration{FT} <: AbstractTranspiration{FT} end
function transpiration(
    transpiration::DiagnosticTranspiration,
    t,
    Y,
    p
)

   @inbounds return p.canopy.hydraulics.fa[length(propertynames(p.canopy.hydraulics.fa))]
end


end
