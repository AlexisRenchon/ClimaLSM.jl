# These plotting functions are meant to be used with Makie
# Can be used on FLUXNET data, as well as ClimaLSM single site output

#= to test functions
using GLMakie
fig = Figure(); ax = Axis(fig[1,1])
x = collect(1:50); y = rand(50);
=# 

# Simple time series
function timeseries!(ax, x, y)
  data = Vec2f.(x,y)
  plt = lines!(ax, data)
  return plt
end

# Seasonality pattern (mean + std, or median + quartiles)

# Diurnal pattern (mean + std, or median + quartiles)

# Energy balance

# Water cumulative

# Light response

# R vs. T & θ

# NEE, Ft Fs vs. u* 

# Light saturated NEE & ET vs. D, by SWC quantiles (color dark to light)

# Other? ...

