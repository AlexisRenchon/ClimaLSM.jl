using JSServe
import JSServe.TailwindDashboard as D
using WGLMakie
WGLMakie.activate!()

using DiffEqCallbacks
import OrdinaryDiffEq as ODE
import ClimaTimeSteppers as CTS
using ClimaCore
import CLIMAParameters as CP
using Statistics
using Dates
using Insolation

using ClimaLSM
using ClimaLSM.Domains: LSMSingleColumnDomain
using ClimaLSM.Soil
using ClimaLSM.Canopy
using ClimaLSM.Canopy.PlantHydraulics
import ClimaLSM
import ClimaLSM.Parameters as LSMP

using ArtifactWrappers
using DelimitedFiles
using Dierckx
using Thermodynamics
using Dates

using ArtifactWrappers
using DelimitedFiles
using ClimaLSM
using Statistics



fig = Figure(); display(fig)

function ClimaLSM_dashboard(button)
  ax = Axis(fig[1,1], xlabel = "day", ylabel = "GPP")
  daily = collect(range(120, 240, 2881))
  GPP = Observable(ones(2881))
  pdata = @lift(Vec2f.(daily, $GPP))
  plt = lines!(ax, pdata)

  on(button) do click
    # include("ozark.jl")
    # GPP[] = [parent(sv.saveval[k].canopy.photosynthesis.GPP)[1] for k in 1:length(sv.saveval)]   
    GPP[] = GPP[] .+ 3
    autolimits!(ax)
  end

  fig
  return fig
end

App() do
    button = D.Button("Run ClimaLSM")
    fig = ClimaLSM_dashboard(button)
    return DOM.div(
                   D.Card(button), D.Card(fig)
                  )
end






# test with GLMakie
# using GLMakie # works
using WGLMakie, JSServe # ok, works too
fig = Figure(); display(fig)
ax = Axis(fig[1,1])
x = [1,2]
y = [4,4]
data = Observable(Vec2f.(x,y))
scatter!(ax, data)
data[] = Vec2f.(x,[6,6])

























