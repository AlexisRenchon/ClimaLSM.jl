const ν = FT(0.495);
const Ksat = FT(0.0443 / 3600 / 100); # m/s
const S_s = FT(1e-3); #inverse meters
const vg_n = FT(2.0);
const vg_α = FT(2.6); # inverse meters
const vg_m = FT(1) - FT(1) / vg_n;
const θ_r = FT(0);
const zmax = FT(0);
const zmin = FT(-10);
const nelems = 50;

soil_domain = Column(FT, zlim = (zmin, zmax), nelements = nelems);
top_flux_bc = FT(-1e-2);
bot_flux_bc = FT(0.0);
boundary_fluxes = (top_flux_bc = top_flux_bc, bot_flux_bc = bot_flux_bc)
params = Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, Ksat, S_s, θ_r);

soil = Soil.RichardsModel{FT}(;
    param_set = params,
    domain = soil_domain,
    boundary_exchanges = boundary_fluxes,
)

Y, p, coords = initialize(soil)

# specify ICs
function init_soil!(Ysoil, z, params)
    function hydrostatic_profile(
        z::FT,
        params::RichardsParameters{FT},
    ) where {FT}
        @unpack ν, vg_α, vg_n, vg_m, θ_r = params
        #unsaturated zone only, assumes water table starts at z_∇
        z_∇ = FT(-10)# matches zmin
        S = FT((FT(1) + (vg_α * (z - z_∇))^vg_n)^(-vg_m))
        ϑ_l = S * (ν - θ_r) + θ_r
        return FT(ϑ_l)
    end
    Ysoil.soil.ϑ_l .= hydrostatic_profile.(z, Ref(params))
end

init_soil!(Y, coords, soil.param_set)

soil_ode! = make_ode_function(soil)

t0 = FT(0);
tf = FT(10);
dt = FT(1);

sv1 = SavedValues(FT, ClimaCore.Fields.FieldVector)
sv2 = SavedValues(FT, Vector{FT})                           
cb = SavingCallback((t, u, integrator) -> (integrator.p), sv1)
cb2 = SavingCallback((t, u, integrator) -> ([parent(integrator.p.soil.ψ)[end]]), sv2)
prob = ODEProblem(soil_ode!, Y, (t0, tf), p);
sol = solve(prob, Euler(); dt = dt, callback = CallbackSet(cb,cb2));
parent(sv1.saveval[2].soil.ψ)[end] - parent(sv1.saveval[1].soil.ψ)[end]
sv2.saveval[2] .- sv2.saveval[1]

@test sum(parent(sol.u[end]) .== parent(Y.soil.ϑ_l)) == nelems
# Testing that ψ +z is constant - > hydrostatic equilibrium
@test mean(parent(p.soil.ψ .+ coords)[:] .- (-10.0)) < eps(FT)
# should be at every layer, at each step too:
@test mean(
    sum([
        parent(saved_values.saveval[k].soil.ψ .+ coords)[:] .+ 10.0 for
        k in 1:1:50
    ]),
) < 1e-14
