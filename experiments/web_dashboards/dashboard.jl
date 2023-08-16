using JSServe
import JSServe.TailwindDashboard as D
using WGLMakie

include("ozark.jl")
include("ozark_met_drivers_FLUXNET.jl");
include("ozark_domain.jl");
include("ozark_parameters.jl");
include("ozark_simulation.jl");

GPP = ClimaLSM_ozark(141.0) .* 1e6 # test - works on 1st call, but subsequent calls return Any[]

fig = Figure(); display(fig) # should be in ClimaLSM_dashboard, but bug for some reason

function ClimaLSM_dashboard(button, g1_input)
  ax = Axis(fig[1,1], xlabel = "day", ylabel = "GPP")
  daily = collect(range(120, 140, 481)) # could be output from ClimaLSM_ozark()
  #GPP = Observable(ones(481))
  #pdata = @lift(Vec2f.(daily, $GPP))
  #plt = lines!(ax, pdata)

  g1 = g1_input.value  
  GPP = @lift(ClimaLSM_ozark($g1)) # empty array...
  pdata = @lift(Vec2f.(daily, $GPP))
  plt = lines!(ax, pdata)

  on(button) do click
    g1 = g1_input.value
    # sv[] = 
    # probably need to work with g1 as an observable
    # GPP[] = [parent(sv.saveval[k].canopy.photosynthesis.GPP)[1] for k in 1:length(sv.saveval)]   
    # GPP[] = ClimaLSM_ozark(g1)
    
    GPP = @lift(ClimaLSM_ozark($g1))
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



