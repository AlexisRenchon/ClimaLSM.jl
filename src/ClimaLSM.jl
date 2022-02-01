module ClimaLSM

using UnPack
using ClimaCore
import ClimaCore: Fields
include("Domains.jl")
using .Domains
import .Domains: coordinates
include("ComponentExchanges.jl")
using .ComponentExchanges

export AbstractModel,
    make_rhs,
    make_ode_function,
    make_update_aux,
    initialize_prognostic,
    initialize_auxiliary,
    initialize,
    LandModel,
    prognostic_vars,
    auxiliary_vars
## Default methods for all models - to be in a seperate module at some point.

abstract type AbstractModel{FT <: AbstractFloat} end

"""
   prognostic_vars(m::AbstractModel)

Returns the prognostic variable symbols for the model in the form of a tuple.
"""
prognostic_vars(m::AbstractModel) = ()

"""
   auxiliary_vars(m::AbstractModel)

Returns the auxiliary variable symbols for the model in the form of a tuple.
"""
auxiliary_vars(m::AbstractModel) = ()

"""
    make_rhs(model::AbstractModel)

Return a `rhs!` function that updates state variables.

`rhs!` should be compatible with OrdinaryDiffEq.jl solvers.
"""
function make_rhs(model::AbstractModel)
    function rhs!(dY, Y, p, t) end
    return rhs!
end

"""
    make_update_aux(model::AbstractModel)

Return an `update_aux!` function that updates auxiliary parameters `p`.
"""
function make_update_aux(model::AbstractModel)
    function update_aux!(p, Y, t) end
    return update_aux!
end

"""
    make_ode_function(model::AbstractModel)

Returns an `ode_function` that updates auxiliary variables and
updates the prognostic state.

`ode_function!` should be compatible with OrdinaryDiffEq.jl solvers.
"""
function make_ode_function(model::AbstractModel)
    rhs! = make_rhs(model)
    update_aux! = make_update_aux(model)
    function ode_function!(dY, Y, p, t)
        update_aux!(p, Y, t)
        rhs!(dY, Y, p, t)
    end
    return ode_function!
end

"""
    initialize_prognostic(model::AbstractModel, state::Union{ClimaCore.Fields.Field, Vector{FT}, FT})

Returns a FieldVector of prognostic variables for `model` with the required
structure, with values equal to `similar(state)`. This assumes that all prognostic
variables are defined over the entire domain. 

The default is a model with no prognostic variables, in which case 
the returned FieldVector contains only an empty array.

Adjustments to this - for example because different prognostic variables
have different dimensions - require defining a new method.
"""
function initialize_prognostic(
    model::AbstractModel{FT},
    state::Union{ClimaCore.Fields.Field, Vector{FT}, FT},
) where {FT}
    keys = prognostic_vars(model)
    if length(keys) == 0
        return ClimaCore.Fields.FieldVector(; model.model_name => FT[])
    else
        values = map((x) -> similar(state), keys)
        return ClimaCore.Fields.FieldVector(;
            model.model_name => (; zip(keys, values)...),
        )
    end

end

"""
    initialize_auxiliary(model::AbstractModel,state::Union{ClimaCore.Fields.Field, Vector{FT}, FT})

Returns a FieldVector of auxiliary variables for `model` with the required
structure, with values equal to `similar(state)`. This assumes that all auxiliary
variables are defined over the entire domain. 

The default is a model with no auxiliary variables, in which case 
the returned FieldVector contains only an empty array.

Adjustments to this - for example because different auxiliary variables
have different dimensions - require defining a new method.
"""
function initialize_auxiliary(
    model::AbstractModel{FT},
    state::Union{ClimaCore.Fields.Field, Vector{FT}, FT},
) where {FT}
    keys = auxiliary_vars(model)
    if length(keys) == 0
        return ClimaCore.Fields.FieldVector(; model.model_name => FT[])
    else
        values = map((x) -> similar(state), keys)
        return ClimaCore.Fields.FieldVector(;
            model.model_name => (; zip(keys, values)...),
        )
    end
end

"""
    initialize(model::AbstractModel)

Creates the prognostic and auxiliary states structures, but with unset values
not set; constructs and returns the coordinates for the `model` domain.
We may need to consider this default more as we add diverse components and 
`Simulations`.
"""
function initialize(model::AbstractModel)
    coords = coordinates(model)
    Y = initialize_prognostic(model, coords)
    p = initialize_auxiliary(model, coords)
    return Y, p, coords
end

Domains.coordinates(model::AbstractModel) = Domains.coordinates(model.domain)

#### LandModel Specific
include("Soil.jl")
include("Roots.jl")
using .Roots
"""
    struct LandModel{FT, SM <: AbstractModel{FT}, RM <: AbstractModel{FT}} <: AbstractModel{FT}

A concrete type of `AbstractModel` for use in land surface modeling. Each component model of the
`LandModel` is itself an `AbstractModel`.

If a user wants to run in standalone, would they use this interface?
No, but should it work ?
"""
struct LandModel{FT, SM <: AbstractModel{FT}, RM <: AbstractModel{FT}} <:
       AbstractModel{FT}
    soil::SM
    roots::RM
end

function LandModel{FT}(; domain::NamedTuple,
                       parameters::NamedTuple,
                       ) where {FT}
    
    boundary_exchanges = LSMExchange{FT}()
    soil = Soil.RichardsModel{FT}(;param_set = parameters.soil,
                             domain = domain.soil,
                             boundary_exchanges =  boundary_exchanges)
    roots = Roots.RootsModel{FT}(;param_set = parameters.roots,
                          domain = domain.roots,
                          boundary_exchanges = boundary_exchanges)
    args = (soil, roots)
    return LandModel{FT, typeof.(args)...}(args...)
end


auxiliary_vars(land::LandModel) = (:root_extraction,)

function initialize(land::LandModel)
    Y_soil, p_soil, coords_soil = initialize(land.soil)
    Y_roots, p_roots, coords_roots = initialize(land.roots)
    # Do we just hardwire all of the interactions? This is so hardcoded, how would
    # we do this for other combos of components?
    Y = ClimaCore.Fields.FieldVector(;
                                     soil = Y_soil.soil,
                                     roots = Y_roots.roots,
                                     )
    p = ClimaCore.Fields.FieldVector(;
                                     root_extraction = similar(Y_soil.soil.ϑ_l),
                                     soil = p_soil.soil,
                                     roots = p_roots.roots,
                                     )
    coords =
        ClimaCore.Fields.FieldVector(; soil = coords_soil, roots = coords_roots)
    return Y, p, coords
end

function make_update_aux(land::LandModel)
    interactions_update_aux! = make_interactions_update_aux(land)
    soil_update_aux! = make_update_aux(land.soil)
    roots_update_aux! = make_update_aux(land.roots)
    function update_aux!(p, Y, t)
        interactions_update_aux!(p, Y, t)
        soil_update_aux!(p, Y, t)
        roots_update_aux!(p, Y, t)
    end
    return update_aux!
end

function make_interactions_update_aux(land::LandModel{FT}) where {FT}
    function update_aux!(p, Y, t)
        root_params = land.roots.param_set
        z = coordinates(land.soil)
        z_up = land.roots.domain.compartment_heights[1]
        rhog_MPa = FT(0.0098)
        ψ = FT(0.0)
        @unpack a_root,
        b_root,
        K_max_root_moles =  land.roots.param_set
        ρm = FT(1e6/18) # moles/m^3
        @. p.root_extraction = compute_flow(z, z_up, rhog_MPa * ψ, Y.roots.rwc[1], a_root, b_root, K_max_root_moles) / ρm # m^3/s need to convert to θ̇...
    end
    return update_aux!
end


function make_ode_function(land::LandModel)
    rhs_soil! = make_rhs(land.soil)
    rhs_roots! = make_rhs(land.roots)
    update_aux! = make_update_aux(land)
    function ode_function!(dY, Y, p, t)
        update_aux!(p, Y, t)
        rhs_soil!(dY, Y, p, t)
        rhs_roots!(dY, Y, p, t)
    end
    return ode_function!
end

end # module
