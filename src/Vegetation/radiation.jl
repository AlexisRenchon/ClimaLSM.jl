export BeerLambertParameters, BeerLambertModel, AlbedoParameters, AlbedoModel, EmissivityParameters, EmissivityModel

abstract type AbstractRadiationModel{FT} <: AbstractModel{FT} end

"""
    BeerLambertParameters{FT <: AbstractFloat}

The required parameters for the Beer-Lambert radiative transfer model.
$(DocStringExtensions.FIELDS)
"""
struct BeerLambertParameters{FT <: AbstractFloat}
    "Leaf angle distribution function (unitless)"
    ld::FT
    "PAR canopy reflectance (unitless)"
    ρ_leaf::FT
    "Clumping index following Braghiere (2021) (unitless)"
    Ω::FT
end

"""
    function BeerLambertParameters{FT}(;
        ld = FT(0.5),    
        ρ_leaf = FT(0.1),
        Ω = FT(1),
    ) where {FT}

A constructor supplying default values for the BeerLambertParameters struct.
"""
function BeerLambertParameters{FT}(;
    ld = FT(0.5),
    ρ_leaf = FT(0.1),
    Ω = FT(1),
) where {FT}
    return BeerLambertParameters{FT}(ld, ρ_leaf, Ω)
end

struct BeerLambertModel{FT} <: AbstractRadiationModel{FT}
    parameters::BeerLambertParameters{FT}
end

#=
"""
"""
struct AlbedoParameters{FT <: AbstractFloat}
    ρ_leaf_sw::FT,
    α_soil::FT,
    Ω::FT
end

"""
"""
function AlbedoParameters{FT}(;
    ρ_leaf_sw = FT(0.1),
    α_soil = FT(0.1),
    Ω = FT(1),
) where {FT}
    return AlbedoParameters{FT}(ρ_leaf_sw, α_soil, Ω)
end

struct AlbedoModel{FT} <: AbstractRadiationModel{FT}
    parameters::AlbedoParameters{FT}
end

struct EmissivityParameters{FT <: AbstractFloat}
    LAI_max::FT,
    k_ϵ::FT,
    fsnow::FT,
    fground::FT,
    ϵ_soil::FT,
    ϵ_veg::FT,
    ϵ_snow::FT
end

"""
"""
function EmissivityParameters{FT}(;
    LAI_max = FT(9),
    k_ϵ = FT(50),
    fsnow = FT(0.0),
    fground = FT(0.5),
    ϵ_soil = FT(0.96),
    ϵ_veg = FT(0.98),
    ϵ_snow = FT(0.97)
) where {FT}
    return EmissivityParameters{FT}(LAI_max, k_ϵ, fsnow, fground, ϵ_soil, ϵ_veg, ϵ_snow)
end

struct EmissivityModel{FT} <: AbstractRadiationModel{FT}
    parameters::EmissivityParameters{FT}
end
=#

ClimaLSM.name(model::AbstractRadiationModel) = :radiative_transfer
ClimaLSM.auxiliary_vars(model::BeerLambertModel) = (:apar,)
ClimaLSM.auxiliary_types(model::BeerLambertModel{FT}) where {FT} = (FT,)

#=
ClimaLSM.auxiliary_vars(model::AlbedoModel) = (:albedo,)
ClimaLSM.auxiliary_types(model::AlbedoModel{FT}) where {FT} = (FT,)

ClimaLSM.auxiliary_vars(model::EmissivityModel) = (:emissivity,)
ClimaLSM.auxiliary_types(model::EmissivityModel{FT}) where {FT} = (FT,)
=#