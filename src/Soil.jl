module Soil

using ClimaLSM
using ClimaCore
import ClimaCore: Fields, Operators, Geometry
import ClimaLSM.Domains: coordinates
using ClimaLSM.ComponentExchanges: LSMExchange, AbstractComponentExchange
import ClimaLSM: AbstractModel, make_update_aux, make_rhs, prognostic_vars, auxiliary_vars
using UnPack
export RichardsModel, RichardsParameters, FluxBC

"""
    RichardsModel

A model for simulating the flow of water in a soil column.
"""
struct RichardsModel{FT, PS, D, C, B} <: AbstractModel{FT}
    param_set::PS
    domain::D
    coordinates::C
    boundary_exchanges::B
    model_name::Symbol
end

function RichardsModel{FT}(; param_set, domain, boundary_exchanges) where {FT}
    coords = coordinates(domain)
    args = (param_set, domain, coords, boundary_exchanges)
    RichardsModel{FT, typeof.(args)...}(
        param_set,
        domain,
        coords,
        boundary_exchanges,
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

struct FluxBC{FT} <: AbstractComponentExchange{FT}
    top_flux_bc::FT
    bot_flux_bc::FT
end

####What if exchange terms always live in p.soil, p.roots?
### and we set p.soil.root_extraction = root source
### p.roots.flow_from_roots = root flow in?
### We'd need an LSM update aux that does both.

function compute_boundary_fluxes(be::FluxBC, _...)
    return be.top_flux_bc, be.bot_flux_bc
end

function compute_boundary_fluxes(be::LSMExchange{FT}, t) where {FT}
    return be.P(t), FT(0.0)
end

function compute_source(be::FluxBC{FT},_...) where {FT}
    return FT(0.0)
end

function compute_source(be::LSMExchange{FT},Y,p) where {FT}
    return p.soil.root_extraction # Field
end

p.soil.top_flux_bc

# p.roots.flow_in_from_roots # single number 

function make_rhs(model::RichardsModel)
    function rhs!(dY, Y, p, t)
        @unpack ν, vg_α, vg_n, vg_m, Ksat, S_s, θ_r = model.param_set
        top_flux_bc, bot_flux_bc = compute_boundary_fluxes(model.boundary_exchanges, t)
        z = model.coordinates
        interpc2f = Operators.InterpolateC2F()
        gradc2f_water = Operators.GradientC2F()
        divf2c_water = Operators.DivergenceF2C(
            top = Operators.SetValue(Geometry.WVector(top_flux_bc)),
            bottom = Operators.SetValue(Geometry.WVector(bot_flux_bc)),
        )
        @. dY.soil.ϑ_l =
            -(divf2c_water(-interpc2f(p.soil.K) * gradc2f_water(p.soil.ψ + z)))
        dY.soil.ϑ_l .+= compute_source(model.boundary_exchanges, Y, p)
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
