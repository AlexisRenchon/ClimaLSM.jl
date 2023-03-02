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
export SharedCanopyParameters, CanopyModel
include("./PlantHydraulics.jl")
using .PlantHydraulics
import .PlantHydraulics: transpiration, AbstractTranspiration
include("./canopy_parameterizations.jl")
include("./stomatalconductance.jl")
include("./photosynthesis.jl")
include("./radiation.jl")


"""
    SharedCanopyParameters{FT <: AbstractFloat, PSE}

A place to store shared parameters that are required by all canopy components.
$(DocStringExtensions.FIELDS)
"""
struct SharedCanopyParameters{FT <: AbstractFloat, PSE}
    "Leaf Area Index"
    LAI::FT
    # "Canopy height"
    # h_c::FT
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

function ClimaLSM.make_update_aux(canopy::CanopyModel{FT}) where {FT}
    plant_hydraulics_update_aux! = make_update_aux(canopy.hydraulics)
    function update_aux!(p, Y, t)
        # update the plant hydraulics system
        plant_hydraulics_update_aux!(p,Y,t)
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

        c_co2::FT = canopy.atmos.c_co2(t)
        P::FT = canopy.atmos.P(t)
        u::FT = canopy.atmos.u(t)
        T::FT = canopy.atmos.T(t)
        h::FT = canopy.atmos.h
        q::FT = canopy.atmos.q(t)
        SW_d::FT = canopy.radiation.shortwave_radiation(t)
        LW_d::FT = canopy.radiation.longwave_radiation(t)  
        θs::FT  = canopy.radiation.zenith_angle(t)      
        # atmos_ts = construct_atmos_ts(canopy.atmos, t, thermo_params)
        # compute VPD
        es = Thermodynamics.saturation_vapor_pressure.(Ref(thermo_params), T, Ref(Thermodynamics.Liquid()))
        ea = @.(q * P / (0.622 + 0.378 * q))
        VPD = es .- ea
        PAR = @.(SW_d / (energy_per_photon * N_a) / 2)

        K = extinction_coeff.(ld, θs)
        APAR = plant_absorbed_ppfd.(PAR, ρ_leaf, K, LAI, Ω)   
        β = moisture_stress.(p.canopy.hydraulics.ψ[1], sc, ψc) 
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

        # compute transpiration and store in p so plant hydraulics can access it.
    end
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
    DiagnosticTranspiration{FT} <: AbstractTranspiration{FT}

A concrete type used for dispatch when computing the transpiration
from the leaves, in the case where transpiration is computed
diagnostically, as a function of prognostic variables and parameters
"""
struct DiagnosticTranspiration{FT} <: AbstractTranspiration{FT} end
function transpiration(
    transpiration::DiagnosticTranspiration{FT},
    t::FT,
    Y,
    p
)::FT where {FT}
#    return the correct value!
end

end
