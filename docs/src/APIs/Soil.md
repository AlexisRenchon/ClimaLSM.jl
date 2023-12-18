# Soil Models

```@meta
CurrentModule = ClimaLSM.Soil
```
## Soil Models

```@docs
ClimaLSM.Soil.AbstractSoilModel
ClimaLSM.Soil.RichardsModel
ClimaLSM.Soil.EnergyHydrology
```
## Soil Parameter Structs

```@docs
ClimaLSM.Soil.RichardsParameters
ClimaLSM.Soil.EnergyHydrologyParameters
```

## Soil Hydrology Parameterizations

```@docs
ClimaLSM.Soil.volumetric_liquid_fraction
ClimaLSM.Soil.pressure_head
ClimaLSM.Soil.hydraulic_conductivity
ClimaLSM.Soil.impedance_factor
ClimaLSM.Soil.viscosity_factor
ClimaLSM.Soil.effective_saturation
ClimaLSM.Soil.matric_potential
ClimaLSM.Soil.dψdϑ
ClimaLSM.Soil.inverse_matric_potential
ClimaLSM.Soil.AbstractSoilHydrologyClosure
ClimaLSM.Soil.vanGenuchten
ClimaLSM.Soil.BrooksCorey
```

## Soil Heat Parameterizations

```@docs
ClimaLSM.Soil.volumetric_heat_capacity
ClimaLSM.Soil.κ_solid
ClimaLSM.Soil.κ_sat_frozen
ClimaLSM.Soil.κ_sat_unfrozen
ClimaLSM.Soil.κ_sat
ClimaLSM.Soil.κ_dry
ClimaLSM.Soil.kersten_number
ClimaLSM.Soil.relative_saturation
ClimaLSM.Soil.volumetric_internal_energy
ClimaLSM.Soil.volumetric_internal_energy_liq
ClimaLSM.Soil.temperature_from_ρe_int
ClimaLSM.Soil.thermal_conductivity
ClimaLSM.Soil.phase_change_source
ClimaLSM.Soil.thermal_time
```

## Soil Surface Parameterizations

```@docs
ClimaLSM.soil.soil_resistance
ClimaLSM.Soil.dry_soil_layer_thickness
ClimaLSM.Soil.soil_tortuosity
```

## Soil Runoff Types and Methods

```@docs
ClimaLSM.Soil.NoRunoff
ClimaLSM.Soil.subsurface_runoff_source
ClimaLSM.Soil.soil_surface_infiltration
```

## Soil BC Methods and Types

```@docs
ClimaLSM.Soil.AbstractSoilBC
ClimaLSM.Soil.MoistureStateBC
ClimaLSM.Soil.FluxBC
ClimaLSM.Soil.TemperatureStateBC
ClimaLSM.Soil.FreeDrainage
ClimaLSM.Soil.RichardsAtmosDrivenFluxBC
ClimaLSM.Soil.AtmosDrivenFluxBC
ClimaLSM.Soil.boundary_vars
ClimaLSM.Soil.boundary_var_domain_names
ClimaLSM.Soil.boundary_var_types
```

## Soil Source Types

```@docs
ClimaLSM.Soil.AbstractSoilSource
ClimaLSM.Soil.PhaseChange
ClimaLSM.Soil.RootExtraction
```

## Soil Jacobian Structures

```@docs
ClimaLSM.Soil.RichardsTridiagonalW
```