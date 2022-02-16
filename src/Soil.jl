module Soil

using ClimaLSM
using ClimaCore
import ClimaCore: Fields, Operators, Geometry
import ClimaLSM.Domains: coordinates
using ClimaLSM.Configurations: RootSoilConfiguration, AbstractConfiguration
import ClimaLSM:
    AbstractModel, make_update_aux, make_rhs, prognostic_vars, auxiliary_vars
using UnPack
export AbstractSoilModel, RichardsModel, RichardsParameters, FluxBC


abstract type AbstractSoilModel{FT} <: AbstractModel{FT} end
"""
    RichardsModel

A model for simulating the flow of water in a soil column.
"""
struct RichardsModel{FT, PS, D, C, B} <: AbstractSoilModel{FT}
    param_set::PS
    domain::D
    coordinates::C
    configuration::B
    model_name::Symbol
end

function RichardsModel{FT}(; param_set, domain, configuration) where {FT}
    coords = coordinates(domain)
    args = (param_set, domain, coords, configuration)
    RichardsModel{FT, typeof.(args)...}(
        param_set,
        domain,
        coords,
        configuration,
        :soil,
    )
end

coordinates(model::RichardsModel) = model.coordinates

"""
    RichardsParameters{FT <: AbstractFloat}

A struct for storing parameters of the Richard's Eq. Soil Model.
"""
struct RichardsParameters{FT <: AbstractFloat}
    ν::FT
    vg_α::FT
    vg_n::FT
    vg_m::FT
    Ksat::FT
    S_s::FT
    θ_r::FT
end

function volumetric_liquid_fraction(ϑ_l::FT, ν_eff::FT) where {FT}
    if ϑ_l < ν_eff
        θ_l = ϑ_l
    else
        θ_l = ν_eff
    end
    return θ_l
end

function matric_potential(α::FT, n::FT, m::FT, S::FT) where {FT}
    ψ_m = -((S^(-FT(1) / m) - FT(1)) * α^(-n))^(FT(1) / n)
    return ψ_m
end

function effective_saturation(porosity::FT, ϑ_l::FT, θr::FT) where {FT}
    ϑ_l_safe = max(ϑ_l, θr + eps(FT))
    S_l = (ϑ_l_safe - θr) / (porosity - θr)
    return S_l
end

function pressure_head(
    α::FT,
    n::FT,
    m::FT,
    θ_r::FT,
    ϑ_l::FT,
    ν_eff::FT,
    S_s::FT,
) where {FT}
    S_l_eff = effective_saturation(ν_eff, ϑ_l, θ_r)
    if S_l_eff <= FT(1.0)
        ψ = matric_potential(α, n, m, S_l_eff)
    else
        ψ = (ϑ_l - ν_eff) / S_s
    end
    return ψ
end

function hydraulic_conductivity(Ksat::FT, m::FT, S::FT) where {FT}
    if S < FT(1)
        K = sqrt(S) * (FT(1) - (FT(1) - S^(FT(1) / m))^m)^FT(2)
    else
        K = FT(1)
    end
    return K * Ksat
end

abstract type AbstractSoilConfiguration{FT} <: AbstractConfiguration{FT} end
# Eventually would support other types of BC
struct FluxBC{FT} <: AbstractSoilConfiguration{FT}
    top_flux_bc::FT
    bot_flux_bc::FT
end

function compute_boundary_fluxes(be::FluxBC, _...)
    return be.top_flux_bc, be.bot_flux_bc
end

function compute_boundary_fluxes(be::RootSoilConfiguration{FT}, t) where {FT}
    return be.P(t), FT(0.0)
end

function compute_source(be::FluxBC{FT}, _...) where {FT}
    return FT(0.0)
end

function compute_source(be::RootSoilConfiguration{FT}, Y, p) where {FT}
    V_layer = FT(1.0 * 0.15) # hardcoded
    ρm = FT(1e6 / 18) # moles/m^3
    return p.root_extraction ./ ρm ./ V_layer # Field
end


function make_rhs(model::RichardsModel)
    function rhs!(dY, Y, p, t)
        @unpack ν, vg_α, vg_n, vg_m, Ksat, S_s, θ_r = model.param_set
        top_flux_bc, bot_flux_bc =
            compute_boundary_fluxes(model.configuration, t)
        z = model.coordinates
        interpc2f = Operators.InterpolateC2F()
        gradc2f_water = Operators.GradientC2F()
        divf2c_water = Operators.DivergenceF2C(
            top = Operators.SetValue(Geometry.WVector(top_flux_bc)),
            bottom = Operators.SetValue(Geometry.WVector(bot_flux_bc)),
        )
        @. dY.soil.ϑ_l =
            -(divf2c_water(-interpc2f(p.soil.K) * gradc2f_water(p.soil.ψ + z)))
        dY.soil.ϑ_l .= compute_source(model.configuration, Y, p)
    end
    return rhs!
end

prognostic_vars(soil::RichardsModel) = (:ϑ_l,)
auxiliary_vars(soil::RichardsModel) = (:K, :ψ)

function make_update_aux(model::RichardsModel)
    function update_aux!(p, Y, t)
        @unpack ν, vg_α, vg_n, vg_m, Ksat, S_s, θ_r = model.param_set
        @. p.soil.K = hydraulic_conductivity(
            Ksat,
            vg_m,
            effective_saturation(ν, Y.soil.ϑ_l, θ_r),
        )
        @. p.soil.ψ = pressure_head(vg_α, vg_n, vg_m, θ_r, Y.soil.ϑ_l, ν, S_s)
        return update_aux!
    end
end

end
