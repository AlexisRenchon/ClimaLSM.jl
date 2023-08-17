using JSServe
import JSServe.TailwindDashboard as D
using WGLMakie

include("ozark.jl")
include("ozark_met_drivers_FLUXNET.jl"); GPP = GPP .* 1e6;
function init() # not sure why needed to rerun each time before using ClimaLSM_ozark(), otherwise return Any[]
  include("ozark_domain.jl");
  include("ozark_parameters.jl");
  include("ozark_simulation.jl");
end
init();

# fig = Figure(); display(fig) # should be in ClimaLSM_dashboard, but bug for some reason

Rn = SW_IN .- SW_OUT .+ LW_IN .- LW_OUT 

function ClimaLSM_dashboard(button, g1_input)
  
  # Create figure and axes 
  fig = Figure(resolution = (1600, 800))
  ax = Axis(fig[1,1], xlabel = "Day of the year", ylabel = "GPP (μmol m⁻² s⁻¹)"); ylims!(ax, (0, 30))
  ax2 = Axis(fig[2,1], xlabel = "Day of the year", ylabel = "Energy (W m⁻²)"); ylims!(ax2, (-150, 900))
  linkxaxes!(ax, ax2)

  # Create model Observables
  daily = collect(range(120, 140, 481)) # could be output from ClimaLSM_ozark()
  GPP_model = Observable(ones(481))
  H_model = Observable(ones(481))
  L_model = Observable(ones(481))
  
  # Plot model
  pdata_GPP = @lift(Vec2f.(daily, $GPP_model))
  pdata_H = @lift(Vec2f.(daily, $H_model))
  pdata_L = @lift(Vec2f.(daily, $L_model))
  plt_GPP_mod = lines!(ax, pdata_GPP)

  plt_H_mod = lines!(ax2, pdata_H)
  plt_L_mod = lines!(ax2, pdata_L)

  # Get parameters value from dashboard
  g1 = g1_input.value 

  # data
  data_daily = collect(range(120, 120+N_days, N_days*48+1))
  data_hh = 120*48:120*48+N_days*48 # data is half-hourly, model is hourly?
  
  # ax 1
  plt_GPP_obs = lines!(ax, data_daily, GPP[data_hh])
  
  # ax 2
  plt_H_obs = lines!(ax2, data_daily, FT.(H[data_hh]))
  plt_LE_obs = lines!(ax2, data_daily, FT.(LE[data_hh]))
  plt_G_obs = lines!(ax2, data_daily, FT.(G[data_hh]))
  plt_Rn_obs = lines!(ax2, data_daily, FT.(Rn[data_hh]))

  # legend
  axislegend(ax, [plt_GPP_obs, plt_GPP_mod], ["GPP model", "GPP data"])
  axislegend(ax2, [plt_Rn_obs, plt_H_obs, plt_LE_obs, plt_G_obs, plt_H_mod, plt_L_mod], ["Rn data", "H data", "L data", "G data", "L model", "H model"])

  # Action when clicking button "Run ClimaLSM"
  on(button) do click
    g1[] = g1_input.value[]
    init(); GPP_model[] = @lift(ClimaLSM_ozark($g1)[:GPP_model])[]
    init(); H_model[] = @lift(ClimaLSM_ozark($g1)[:H_model])[]
    init(); L_model[] = @lift(ClimaLSM_ozark($g1)[:L_model])[]
    # autolimits!(ax)
  end

  fig
  return fig
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



# TO DO: add a card or column right to GPP, put energy flux in it 
# top row: diurnal
# bottom row: seasonal
