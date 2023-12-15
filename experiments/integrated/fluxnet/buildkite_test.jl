using ClimaLSM

ARGS = ["US-MOz"]

climalsm_dir = pkgdir(ClimaLSM)

include(joinpath(climalsm_dir, "experiments/integrated/fluxnet/run_fluxnet.jl"))
include(joinpath(climalsm_dir, "experiments/integrated/fluxnet/inputs_dataframe.jl"))
include(joinpath(climalsm_dir, "experiments/integrated/fluxnet/climalsm_output_dataframe.jl")) 

if isdir(joinpath(climalsm_dir, "experiments/integrated/fluxnet/figures"))
    nothing
else
    mkdir(joinpath(climalsm_dir, "experiments/integrated/fluxnet/figures"))
end

include(joinpath(climalsm_dir, "experiments/integrated/fluxnet/makie_plots.jl"))

fig1 = timeseries_fig(inputs, climalsm)
fig2 = fingerprint_fig(inputs, climalsm)
fig3 = diurnals_fig(inputs, climalsm)

names = ["timeseries.pdf", "fingerprint.pdf", "diurnals.pdf"]

[save(joinpath(savedir, name), fig) for (name, fig) in zip(names, [fig1, fig2, fig3])]

