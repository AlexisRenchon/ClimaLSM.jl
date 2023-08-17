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

function ClimaLSM_dashboard(button, g1_input)
  fig = Figure()
  ax = Axis(fig[1,1], xlabel = "Day of the year", ylabel = "GPP (μmol m⁻² s⁻¹)")
  daily = collect(range(120, 140, 481)) # could be output from ClimaLSM_ozark()
  GPP_model = Observable(ones(481))
  pdata = @lift(Vec2f.(daily, $GPP_model))
  plt = lines!(ax, pdata)
  g1 = g1_input.value 

  # GPP data
  data_daily = collect(range(120, 120+N_days, N_days*48+1))
  data_hh = 120*48:120*48+N_days*48 # data is half-hourly, model is hourly?
  plt_GPP = lines!(ax, data_daily, GPP[data_hh])

  # legend
  axislegend(ax, [plt, plt_GPP], ["GPP model", "GPP data"])

  # Action when clicking button "Run ClimaLSM"
  on(button) do click
    g1[] = g1_input.value[]
    init();
    GPP_model[] = @lift(ClimaLSM_ozark($g1))[]
    autolimits!(ax)
  end

  fig
  return fig
end

App() do
    button = D.Button("Run ClimaLSM")
    g1_input = D.NumberInput(141.0)
    fig = ClimaLSM_dashboard(button, g1_input)
    return DOM.div(
                   D.Card(button), D.Card(D.FlexRow("g1: ", g1_input)), D.Card(fig)
                  )
end


# TO DO: add a card or column right to GPP, put energy flux in it 
# top row: diurnal
# bottom row: seasonal
