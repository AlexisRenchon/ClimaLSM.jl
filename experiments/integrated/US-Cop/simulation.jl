global t0 
global tf 
global dt
global n 
global N_days=round((tf .- t0)./FT(3600 * 24))

saveat = Array(t0:(n * dt[1]):tf)
timestepper = CTS.RK4()
# Set up timestepper
ode_algo = CTS.ExplicitAlgorithm(timestepper)



