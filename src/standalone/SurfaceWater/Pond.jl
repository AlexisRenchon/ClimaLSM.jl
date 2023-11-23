module Pond
using ClimaLSM
using ClimaCore
using DocStringExtensions
import ClimaCore: Fields

import ClimaLSM:
    AbstractExpModel,
    make_compute_exp_tendency,
    prognostic_vars,
    name,
    prognostic_types,
    prognostic_domain_names
using ClimaLSM.Domains
export PondModel, PrescribedRunoff, surface_runoff

abstract type AbstractSurfaceWaterModel{FT} <: AbstractExpModel{FT} end
abstract type AbstractSurfaceRunoff end
"""
    PondModel{FT, D, R} <: AbstractSurfaceWaterModel{FT}

A stand-in model for models like the snow or river model. In
standalone mode, a prescribed soil infiltration rate
 and precipitation rate
control the rate of change of the pond height variable `η` via an ODE.
In integrated LSM mode, the infiltration into the soil will be computed
via a different method, and also be applied as a flux boundary condition
for the soil model.
$(DocStringExtensions.FIELDS)
"""
struct PondModel{FT, D, R} <: AbstractSurfaceWaterModel{FT}
    "The domain for the pond model"
    domain::D
    "The runoff model for the pond model"
    runoff::R
end

function PondModel{FT}(;
    domain::ClimaLSM.Domains.AbstractDomain{FT} = ClimaLSM.Domains.Point(
        z_sfc = FT(0),
    ),
    runoff::AbstractSurfaceRunoff,
) where {FT}
    return PondModel{FT, typeof(domain), typeof(runoff)}(domain, runoff)
end


"""
    PrescribedRunoff{FT} <:  AbstractSurfaceRunoff

The required input for driving the simple pond model: precipitation, as a
function of time, soil effective saturation at a depth `Δz` below the surface,
as a function of time, and soil parameters, which affect infiltration.
"""
struct PrescribedRunoff{FT} <: AbstractSurfaceRunoff
    "Time dependent precipitation magnitude, given in m/s. Negative is into the soil"
    precip::Function
    "Time dependent infiltration magnitude, given in m/s. Negative is into the soil."
    infil::Function
end

ClimaLSM.prognostic_vars(model::PondModel) = (:η,)
ClimaLSM.prognostic_types(model::PondModel{FT}) where {FT} = (FT,)
ClimaLSM.prognostic_domain_names(model::PondModel) = (:surface,)
ClimaLSM.name(::AbstractSurfaceWaterModel) = :surface_water

function ClimaLSM.make_compute_exp_tendency(model::PondModel)
    function compute_exp_tendency!(dY, Y, p, t)
        runoff = surface_runoff(model.runoff, Y, p, t)
        @. dY.surface_water.η = runoff
    end
    return compute_exp_tendency!
end

# Runoff > 0 -> into river system.
function surface_runoff(runoff::PrescribedRunoff{FT}, Y, p, t) where {FT}
    return @. -FT((runoff.precip(t) - runoff.infil(t)))
end

end
