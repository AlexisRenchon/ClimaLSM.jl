export FarquharParameters, FarquharModel, C3, C4

abstract type AbstractPhotosynthesisMechanism end
"""
    C3 <: AbstractPhotosynthesisMechanism

Helper struct for dispatching between C3 and C4 photosynthesis.
"""
struct C3 <: AbstractPhotosynthesisMechanism end

"""
    C4 <: AbstractPhotosynthesisMechanism

Helper struct for dispatching between C3 and C4 photosynthesis.
"""
struct C4 <: AbstractPhotosynthesisMechanism end

abstract type AbstractPhotosynthesisModel{FT} <: AbstractCanopyComponent{FT} end
"""
    FarquharParameters{FT<:AbstractFloat}

The required parameters for the Farquhar photosynthesis model.
$(DocStringExtensions.FIELDS)
"""
struct FarquharParameters{FT <: AbstractFloat}
    "Photosynthesis mechanism: C3 or C4"
    mechanism::AbstractPhotosynthesisMechanism
    "Vcmax at 25 °C (mol CO2/m^2/s)"
    Vcmax25::FT
    "Γstar at 25 °C (mol/mol)"
    Γstar25::FT
    "Michaelis-Menten parameter for CO2 at 25 °C (mol/mol)"
    Kc25::FT
    "Michaelis-Menten parameter for O2 at 25 °C (mol/mol)"
    Ko25::FT
    "Energy of activation for CO2 (J/mol)"
    ΔHkc::FT
    "Energy of activation for oxygen (J/mol)"
    ΔHko::FT
    "Energy of activation for Vcmax (J/mol)"
    ΔHVcmax::FT
    "Energy of activation for Γstar (J/mol)"
    ΔHΓstar::FT
    "Energy of activation for Jmax (J/mol)"
    ΔHJmax::FT
    "Energy of activation for Rd (J/mol)"
    ΔHRd::FT
    "Reference temperature equal to 25 degrees Celsius (K)"
    To::FT
    "Intercelluar O2 concentration (mol/mol); taken to be constant"
    oi::FT
    "Quantum yield of photosystem II (Bernacchi, 2003; unitless)"
    ϕ::FT
    "Curvature parameter, a fitting constant to compute J, unitless"
    θj::FT
    "Constant factor appearing the dark respiration term, equal to 0.015."
    f::FT
    "Fitting constant to compute the moisture stress factor (Pa^{-1})"
    sc::FT
    "Fitting constant to compute the moisture stress factor (Pa)"
    ψc::FT
end

"""
    function FarquharParameters{FT}(mechanism::AbstractPhotosynthesisMechanism;
        oi = FT(0.209),# mol/mol
        ϕ = FT(0.6), # unitless
        θj = FT(0.9), # unitless
        f = FT(0.015), # unitless
        sc = FT(5e-6),# Pa
        ψc = FT(-2e6), # Pa
        Vcmax25 = FT(5e-5), # converted from 50 μmol/mol CO2/m^2/s to mol/m^2/s
        Γstar25 = FT(4.275e-5),  # converted from 42.75 μmol/mol to mol/mol
        Kc25 = FT(4.049e-4), # converted from 404.9 μmol/mol to mol/mol
        Ko25 = FT(0.2874), # converted from 278.4 mmol/mol to mol/mol
        To = FT(298.15), # 25 C
        ΔHkc = FT(79430), #J/mol, Table 11.2 Bonan
        ΔHko = FT(36380), #J/mol, Table 11.2 Bonan
        ΔHVcmax = FT(58520), #J/mol, Table 11.2 Bonan
        ΔHΓstar = FT(37830), #J/mol, 11.2 Bonan
        ΔHJmax = FT(43540), # J/mol, 11.2 Bonan
        ΔHRd = FT(43390), # J/mol, 11.2 Bonan
        ) where {FT}

A constructor supplying default values for the FarquharParameters struct.
"""
function FarquharParameters{FT}(
    mechanism::AbstractPhotosynthesisMechanism;
    oi = FT(0.209),
    ϕ = FT(0.6),
    θj = FT(0.9),
    f = FT(0.015),
    sc = FT(5e-6),
    ψc = FT(-2e6),
    Vcmax25 = FT(5e-5),
    Γstar25 = FT(4.275e-5),
    Kc25 = FT(4.049e-4),
    Ko25 = FT(0.2874),
    To = FT(298.15),
    ΔHkc = FT(79430),
    ΔHko = FT(36380),
    ΔHVcmax = FT(58520),
    ΔHΓstar = FT(37830),
    ΔHJmax = FT(43540),
    ΔHRd = FT(46390),
) where {FT}
    return FarquharParameters{FT}(
        mechanism,
        Vcmax25,
        Γstar25,
        Kc25,
        Ko25,
        ΔHkc,
        ΔHko,
        ΔHVcmax,
        ΔHΓstar,
        ΔHJmax,
        ΔHRd,
        To,
        oi,
        ϕ,
        θj,
        f,
        sc,
        ψc,
    )
end

struct FarquharModel{FT} <: AbstractPhotosynthesisModel{FT}
    parameters::FarquharParameters{FT}
end

ClimaLSM.name(model::AbstractPhotosynthesisModel) = :photosynthesis
ClimaLSM.auxiliary_vars(model::FarquharModel) = (:An, :GPP)
ClimaLSM.auxiliary_types(model::FarquharModel{FT}) where {FT} = (FT, FT)
