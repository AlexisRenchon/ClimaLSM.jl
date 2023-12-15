#= run these lines for test / development
ARGS = ["US-MOz"]
include("integrated/fluxnet/run_fluxnet.jl")
include("integrated/fluxnet/inputs_dataframe.jl")
include("integrated/fluxnet/climalsm_output_dataframe.jl") 
# TO DO:
# change unit from SI to what we want (e.g., umol m-2 s-1)
=# 

# Make all sort of plots with data and model output
# 1. Time series (e.g., C fluxes, h2o fluxes, energy fluxes, met drivers)
# 2. Seasonal pattern (i.e., monthly average and std)
# 3. Diurnal pattern (i.e., hourly average and std)
# 4. Response curves (e.g., NEE vs. PAR with VPD color and SWC brightness...)  
# 5. Energy conservation (i.e., Rn - G vs. L + H)
# 6. Water budget (i.e., cumulative ET vs. P)
# 7. Data quality plots (e.g., NEE vs. u* by T and SWC bins)
# 8. Fingerprint plots (showing both seasonality and diurnal pattern)
# 9. Wavelet coherence
# to do in another script: animations

# TO DO:  
# Script for plot utilities, Ed started this with plot_utils.jl

# TO DO:
# make white or dark background figures
# publication style and presentation style (bigger font etc.)

using ClimaLSM 
climalsm_dir = pkgdir(ClimaLSM)
savedir = joinpath(climalsm_dir, "experiments", "integrated", "fluxnet/figures/") 
using CairoMakie # Draw vector graphics to SVG or PDF. High quality plots! 
using LaTeXStrings # To have latex labels
using PlotUtils: optimize_ticks

# drivers will be in drivers from Ed PR

# 1. Time series of GPP, ET and SW_OUT
function timeseries_fig(inputs, climalsm) # will run for any inputs or climalsm output of FLUXNET sites
    # create an empty figure
    fig = Figure(size = (1000, 1000)) # note: do not load Plots.jl in this branch (it is loading in plot_utils)
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    # create empty axis, with a specific layout
    ax_C = Axis(fig[1, 1], ylabel = L"\text{GPP} \, (\mu\text{mol m}^{-2} \, \text{s}^{-1})") # C fluxes
    ax_W = Axis(fig[2, 1], ylabel = "ET (mm)") # h2o fluxes
    ax_SWOUT = Axis(fig[3, 1], ylabel = L"\text{SW OUT} \, (\text{W} \, \text{m}^{-2})") # shortwave out 
    # ax_T = Axis(fig[4, 1]) # air, canopy, and soil temperature

    # for time series, Makie should allow DateTime type soon (but not yet)
    # so the 2 lines of code below are a trick to be able to use DateTime - will be removed later
    dateticks = optimize_ticks(climalsm.DateTime[1], climalsm.DateTime[end])[1][2:end-1] # first and last are weirdly placed

    # add plots into axis ax_C
    p_GPP_m = CairoMakie.scatter!(ax_C, datetime2unix.(climalsm.DateTime), climalsm.GPP .* 1e6, color = :blue)
    p_GPP_d = CairoMakie.scatter!(ax_C, datetime2unix.(inputs.DateTime[index_t_start:index_t_end]), inputs.GPP[index_t_start:index_t_end] .*1e6, color = :black) 

    # ax_W
    p_ET_m = CairoMakie.scatter!(ax_W, datetime2unix.(climalsm.DateTime), (climalsm.vapor_flux .* 1e3 .* 24 .* 3600) .+ (climalsm.transpiration .* 1e3 .* 24 .* 3600), color = :blue) # not sure about units
    p_ET_d = CairoMakie.scatter!(ax_W, datetime2unix.(inputs.DateTime[index_t_start:index_t_end]), inputs.LE[index_t_start:index_t_end] ./ (LSMP.LH_v0(earth_param_set) * 1000) .* (1e3 * 24 * 3600), color = :black) # not sure units

    # ax_SW_OUT
    p_SWOUT_m = CairoMakie.scatter!(ax_SWOUT, datetime2unix.(climalsm.DateTime), climalsm.SW_out, color = :blue)
    p_SWOUT_d = CairoMakie.scatter!(ax_SWOUT, datetime2unix.(inputs.DateTime[index_t_start:index_t_end]), FT.(inputs.SW_OUT[index_t_start:index_t_end]), color = :black) 

    # xticks
    ax_C.xticks[] = (datetime2unix.(dateticks) , Dates.format.(dateticks, "mm/dd"))
    ax_W.xticks[] = (datetime2unix.(dateticks) , Dates.format.(dateticks, "mm/dd"))
    ax_SWOUT.xticks[] = (datetime2unix.(dateticks) , Dates.format.(dateticks, "mm/dd"))

    axislegend(ax_C, [p_GPP_d, p_GPP_m], ["Observations", "ClimaLSM"], "", position = :rt, orientation = :horizontal)

    CairoMakie.ylims!(ax_C, (0, 40))
    CairoMakie.ylims!(ax_W, (0, 30))
    CairoMakie.ylims!(ax_SWOUT, (0, 200))

    [CairoMakie.xlims!(axes, (datetime2unix(climalsm.DateTime[1]), datetime2unix(climalsm.DateTime[end]))) for axes in [ax_C, ax_W, ax_SWOUT]]

    fig
    return fig
end

#= to test. These will be in another script though.
fig = timeseries_fig(inputs, climalsm) 
save(joinpath(savedir, "timeseries.pdf"), fig)
=#


# 2. Fingerprint plot
function fingerprint_fig(inputs, climalsm)
    fig = Figure(size = (1000, 1000))
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    # create empty axis, with a specific layout
    ax_C = Axis(fig[1, 1], ylabel = "Hour of the day", xlabel = "Date", title = L"\text{GPP} \, (\mu\text{mol m}^{-2} \, \text{s}^{-1})") # C fluxes
    ax_W = Axis(fig[2, 1]) # h2o fluxes
    ax_P = Axis(fig[3, 1]) # Precip & soil moisture
    ax_T = Axis(fig[4, 1]) # air, canopy, and soil temperature

    # for time series, Makie should allow DateTime type soon (but not yet)
    # so the 2 lines of code below are a trick to be able to use DateTime - will be removed later
    dateticks = optimize_ticks(inputs.DateTime[1], inputs.DateTime[end])[1][2:end-1] # first and last are weirdly placed

    # Fingerprint plot
    hm_GPP = CairoMakie.heatmap!(ax_C, datetime2unix.(DateTime.(Date.(inputs.DateTime))), hour.(inputs.DateTime) .+ (minute.(inputs.DateTime) ./ 60), inputs.GPP .* 1e6)
    Colorbar(fig[1, 2], hm_GPP)

    ax_C.xticks[] = (datetime2unix.(dateticks) , Dates.format.(dateticks, "mm/dd"))

    fig
    return fig
end

#=
fig = fingerprint_fig(inputs, climalsm)
save(joinpath(savedir, "fingerprint.pdf"), fig)
=#

# 3. Diurnals, with quantiles, for C, h2o, energy
using Statistics

function diurnal(datetime, data)
    hourlyquantile = [quantile(data[hour.(datetime) .== h], [0.25, 0.5, 0.75]) for h = 0:1:23]
    return hourlyquantile 
end

function diurnals_fig(inputs, climalsm)
    fig = Figure(size = (1000, 1000))
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    ax_C = Axis(fig[1, 1], ylabel = L"\text{CO}_{2} \, (\mu\text{mol m}^{-2} \, \text{s}^{-1})") # C fluxes
    ax_W = Axis(fig[2, 1], ylabel = L"\text{H}_{2}\text{O} \, \text{(mm)}") # h2o fluxes
    ax_E = Axis(fig[3, 1], ylabel = L"\text{Radiation} \, (\text{W} \, \text{m}^{-2})") # shortwave out 

    # to do: make a diurnalplot function to do this line and band, then loop over all the variables - should be 3 lines, not 20 lines...
    GPP_diurnal = diurnal(climalsm.DateTime, climalsm.GPP)
    GPP_diurnal_p = CairoMakie.lines!(ax_C, 0.5:1:23.5, getindex.(GPP_diurnal[1:24], 2) .* 1e6, color = :black)
    GPP_diurnal_q = CairoMakie.band!(ax_C, 0.5:1:23.5, getindex.(GPP_diurnal[1:24], 1) .* 1e6, getindex.(GPP_diurnal[1:24], 3) .* 1e6, color = (:green, 0.3))

    ET_diurnal = diurnal(climalsm.DateTime, climalsm.transpiration)
    ET_diurnal_p = CairoMakie.lines!(ax_W, 0.5:1:23.5, getindex.(ET_diurnal[1:24], 2) .* 1e3 .* 24 .* 3600, color = :black)
    ET_diurnal_q = CairoMakie.band!(ax_W, 0.5:1:23.5, getindex.(ET_diurnal[1:24], 1) .* 1e3 .* 24 .* 3600, getindex.(ET_diurnal[1:24], 3) .* 1e3 .* 24 .* 3600, color = (:blue, 0.3))

    [CairoMakie.xlims!(axes, (0, 24)) for axes in [ax_C, ax_W, ax_E]]

    fig
    return fig
end

#=
fig = diurnals_fig(inputs, climalsm)
save(joinpath(savedir, "diurnals.pdf"), fig)
=#


