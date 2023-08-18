using JSServe
import JSServe.TailwindDashboard as D
using WGLMakie

include("ozark.jl")
include("ozark_met_drivers_FLUXNET.jl"); GPP = GPP .* 1e6;
function init()
    include("ozark_domain.jl")
    include("ozark_parameters.jl")
    include("ozark_simulation.jl")
end

init()

# fig = Figure(); display(fig) # should be in ClimaLSM_dashboard, but bug for some reason
Rn = SW_IN .- SW_OUT .+ LW_IN .- LW_OUT

function ClimaLSM_dashboard(button, g1_input)
    fig = Figure(resolution = (1600, 800))
    ax, ax2 = Axis(fig[1,1], xlabel = "Day of the year", ylabel = "GPP (μmol m⁻² s⁻¹)"), Axis(fig[2,1], xlabel = "Day of the year", ylabel = "Energy (W m⁻²)")
    ylims!(ax, (-2, 20))
    ylims!(ax2, (-150, 900))
    linkxaxes!(ax, ax2)

    n = N_days*24+1
    daily = collect(range(1, N_days + 1, n))

    g1 = g1_input.value
    out = Observable(ClimaLSM_ozark(g1[]))
    GPP_model, H_model, L_model = Observable(out[][:GPP_model]), Observable(out[][:H_model]), Observable(out[][:L_model])

    pdata_GPP, pdata_H, pdata_L = @lift(Vec2f.(daily, $GPP_model)), @lift(Vec2f.(daily, $H_model)), @lift(Vec2f.(daily, $L_model))
    plt_GPP_mod, plt_H_mod, plt_L_mod = lines!(ax, pdata_GPP), lines!(ax2, pdata_H), lines!(ax2, pdata_L)

    data_daily, data_hh = collect(range(1, 1+N_days, N_days*48+1)), 1*48:1*48+N_days*48

    plt_GPP_obs = lines!(ax, data_daily, GPP[data_hh])
    plt_H_obs, plt_LE_obs, plt_G_obs, plt_Rn_obs = lines!(ax2, data_daily, FT.(H[data_hh])), lines!(ax2, data_daily, FT.(LE[data_hh])), lines!(ax2, data_daily, FT.(G[data_hh])), lines!(ax2, data_daily, FT.(Rn[data_hh]))

    axislegend(ax, [plt_GPP_obs, plt_GPP_mod], ["GPP model", "GPP data"])
    axislegend(ax2, [plt_Rn_obs, plt_H_obs, plt_LE_obs, plt_G_obs, plt_H_mod, plt_L_mod], ["Rn data", "H data", "L data", "G data", "L model", "H model"])

    on(button) do click
        g1[] = g1_input.value[]
        init()
        out[] = @lift(ClimaLSM_ozark($g1))[]
        GPP_model[], H_model[], L_model[] = out[][:GPP_model], out[][:H_model], out[][:L_model]
    end

    fig
end

App() do
    button = D.Button("Run ClimaLSM")
    g1_input = D.NumberInput(141.0)
    fig = ClimaLSM_dashboard(button, g1_input)
    return DOM.div(
        D.Card(button),
        D.Card(D.FlexRow("g1: ", g1_input)),
        D.Card(fig)
    )
end

