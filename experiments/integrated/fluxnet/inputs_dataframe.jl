#= run those lines for testing / development
ARGS = ["US-MOz"]
include("integrated/fluxnet/setup.jl")
include("integrated/fluxnet/run_fluxnet.jl")
=#

using DataFrames
variables_name = (
    "DateTime",
    "RECO",
    "TA",
    "VPD",
    "PA",
    "P",
    "WS",
    "LW_IN",
    "SW_IN",
    "CO2",
    "SWC",
    "TS",
    "GPP",
    "LE",
    "H",
    "G",
    "SW_OUT",
) # I imagine this may not work for all fluxnet sites...

variables = [
    vec(mat) for mat in (
        LOCAL_DATETIME,
        drivers.RECO.values,
        drivers.TA.values,
        drivers.VPD.values,
        drivers.PA.values,
        drivers.P.values,
        drivers.WS.values,
        drivers.LW_IN.values,
        drivers.SW_IN.values,
        drivers.CO2.values,
        drivers.SWC.values,
        drivers.TS.values,
        drivers.GPP.values,
        drivers.LE.values,
        drivers.H.values,
        drivers.G.values,
        drivers.SW_OUT.values,
    )
]

inputs = DataFrame([
    variables_name[i] => variables[i] for i in 1:length(variables_name)
])
