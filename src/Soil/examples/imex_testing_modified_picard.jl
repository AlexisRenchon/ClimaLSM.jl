using ClimaCore
using DiffEqBase
import OrdinaryDiffEq as ODE
import ClimaTimeSteppers as CTS
using StaticArrays

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Soil
using ClimaLSM.Domains: Column
include("./TridiagonalJacobian.jl")
using .TridiagonalJacobian:
    TridiagonalW, make_Wfact, make_implicit_tendency, explicit_tendency!

FT = Float64

is_imex_CTS_algo(::CTS.IMEXAlgorithm) = true
is_imex_CTS_algo(::DiffEqBase.AbstractODEAlgorithm) = false

# is_implicit(::ODE.OrdinaryDiffEqImplicitAlgorithm) = true
# is_implicit(::ODE.OrdinaryDiffEqAdaptiveImplicitAlgorithm) = true
# is_implicit(ode_algo) = is_imex_CTS_algo(ode_algo)

is_rosenbrock(::ODE.Rosenbrock23) = true
is_rosenbrock(::ODE.Rosenbrock32) = true
is_rosenbrock(::DiffEqBase.AbstractODEAlgorithm) = false
use_transform(ode_algo) =
    !(is_imex_CTS_algo(ode_algo) || is_rosenbrock(ode_algo))
stepper = CTS.ARS111()
norm_condition = CTS.MaximumError(Float64(1e-8))
conv_checker = CTS.ConvergenceChecker(; norm_condition)
# backeuler_cn_tableau = CTS.IMEXTableau(;
#     a_exp = @SArray([0 0; 1 0]), # not being used (no exp tendency)
#     a_imp = @SArray([1 0; 0 0]),
#     b_imp = @SArray([1/2, 1/2]),
#     c_imp = @SArray([1, 0])
# )

# ode_algo = CTS.IMEXAlgorithm(
#     backeuler_cn_tableau,
#     CTS.NewtonsMethod(
#         max_iters = 500,
#         update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
#         convergence_checker = conv_checker,
#         verbose = CTS.Verbose()
#     ),
# )
ode_algo = CTS.IMEXAlgorithm(
    stepper,
    CTS.NewtonsMethod(
        max_iters = 1,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
       # convergence_checker = conv_checker,
       # verbose = CTS.Verbose()
    ),
)

#function main(ode_algo, t_end::Float64, dt::Float64; explicit=false)

t_start = FT(0)
t_end = FT(1e6)
dt = FT(1e5) #1000 = 1e3, 10000 = 1e4
explicit = false

# van Genuchten parameters for clay (from Bonan 2019 supplemental program 8.2)
ν = FT(0.495)
K_sat = FT(0.0443 / 3600 / 100) # m/s
vg_n = FT(1.43)
vg_α = FT(0.026 * 100) # inverse meters
vg_m = FT(1) - FT(1) / vg_n
θ_r = FT(0.124)
S_s = FT(1e-3) #inverse meters

zmax = FT(0)
zmin = FT(-1.5)
nelems = 150

soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems);

top_bc = Soil.MoistureStateBC((p, t) -> eltype(t)(ν - 1e-3))
#top_bc = Soil.FluxBC((p, t) -> eltype(t)(-1.23e-7))

# function flux_function(p, t::T) where {T}
#     if t < 1e5
#         return -T(1e-6) # TODO find cutoff time
#     else
#         return T(0)
#     end
# end
# top_bc = Soil.FluxBC(flux_function)
top_bc = Soil.MoistureStateBC((p, t) -> eltype(t)(ν - 1e-3))
flux_in = FT(-1e-7)
#top_bc = Soil.FluxBC((p, t) -> eltype(t)(flux_in))


# bot_bc = Soil.FreeDrainage()
flux_out = FT(0)
bot_bc = Soil.FluxBC((p, t) -> eltype(t)(flux_out))

sources = ()
boundary_fluxes = (; top = (water = top_bc,), bottom = (water = bot_bc,))
params = Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, K_sat, S_s, θ_r)

soil = Soil.RichardsModel{FT}(;
    parameters = params,
    domain = soil_domain,
    boundary_conditions = boundary_fluxes,
    sources = sources,
)

Y, p, coords = initialize(soil)
@. Y.soil.ϑ_l = FT(0.24)

if !explicit
    transform = use_transform(ode_algo)

    W = TridiagonalW(Y, transform)

    implicit_tendency! = make_implicit_tendency(soil)
    Wfact! = make_Wfact(soil)
    Wfact!(W, Y, p, dt, t_start)

    jac_kwargs = if use_transform(ode_algo)
        (; jac_prototype = W, Wfact_t = Wfact!)
    else
        (; jac_prototype = W, Wfact = Wfact!)
    end

    implicit_problem = ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = explicit_tendency!,
            T_imp! = ODEFunction(implicit_tendency!; jac_kwargs...),
            dss! = ClimaLSM.Soil.dss!,
        ),
        Y,
        (t_start, t_end),
        p,
    )

    integrator = init(
        implicit_problem,
        ode_algo;
        dt = dt,
        adaptive = false,
        progress = true,
        saveat = t_start:1:t_end,
    )
else
    implicit_tendency! = make_implicit_tendency(soil)
    explicit_problem = ODE.ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = implicit_tendency!,
            T_imp! = nothing,
            dss! = ClimaLSM.Soil.dss!,
        ),
        Y,
        (t_start, t_end),
        p,
    )
    integrator = init(
        explicit_problem,
        ode_algo;
        dt = dt,
        adaptive = false,
        progress = true,
        saveat = t_start:1000:t_end,
    )
end



plot(parent(integrator.sol.u[end].soil.ϑ_l), parent(coords.z))

# for step in 1:t_end
#     @show step
#     ODE.step!(integrator)
# end
ODE.solve!(integrator)

plot(parent(integrator.sol.u[end].soil.ϑ_l), parent(coords.z))

# calculate water mass balance over entire simulation
mass_end = sum(integrator.sol.u[end].soil.ϑ_l)
mass_start = sum(integrator.sol.u[1].soil.ϑ_l)
t_sim = integrator.sol.t[end] - integrator.sol.t[1]
# flux changes water content every timestep (assumes constant flux_in, flux_out)
mass_change_exp = -(flux_in - flux_out) * t_sim
mass_change_actual = mass_end - mass_start
relerr = abs(mass_change_actual - mass_change_exp) / mass_change_exp * 100 # %

@show relerr
@show integrator.sol.u[end]