# Astronomical Physics Units Library

A Rust library providing strongly-typed physical units for astronomical calculations with built-in unit conversion.

## Overview

This library implements a type-safe way to represent common astronomical measurements and handle unit conversions. All types use `f64` precision for scientific computing accuracy.

### Basic Physical Quantities
- **Length**: AU, Earth radii, Jupiter radii, solar radii, km, m, cm, microns
- **Mass**: Solar masses, Earth masses, Jupiter masses, grams
- **Time**: Years, million years (Myr), days, hours, seconds
- **Temperature**: Kelvin, Celsius, Fahrenheit
- **Velocity**: AU/year, cm/sec

### Derived Quantities (for planet formation simulations)
- **MassRate**: Mass flow rates (M☉/year, M⊕/Myr) for accretion and disk evolution
- **SurfaceDensity**: Disk surface densities (g/cm², kg/m²) with power-law scaling
- **VolumeDensity**: Bulk densities (g/cm³) with material constants and composition tools

## Features

- **Type-safe arithmetic**: Compile-time prevention of dimensional analysis errors
- **f64 precision**: Suitable for scientific computing
- **Zero-cost abstractions**: No runtime overhead compared to raw floats
- **Scientifically accurate**: Conversion factors from authoritative sources
- **Comprehensive tests**: Full unit test and doctest coverage
- **Standard operators**: Implements `+`, `-`, `*`, `/` for supported types
- **Domain-specific helpers**: Physics-aware methods like `integrate()`, `power_law_scaling()`, `weighted_average()`

## Usage Examples

### Length Conversions

```rust
use units::Length;

// Create lengths in various units
let orbit = Length::from_au(1.0);
let earth_radius = Length::from_earth_radii(1.0);
let jupiter_radius = Length::from_jupiter_radii(1.0);
let stellar_radius = Length::from_solar_radii(1.0);

// Convert between units
let cm_value = orbit.to_cm();           // 1.496e13 cm
let km_value = orbit.to_km();           // 1.496e8 km
let earth_r = orbit.to_earth_radii();   // ~23481 R⊕
```

### Mass Conversions

```rust
use units::Mass;

// Create masses in different units
let star = Mass::from_solar_masses(1.0);
let planet = Mass::from_earth_masses(317.8);  // Jupiter
let rock = Mass::from_grams(1e24);

// Convert between units
let grams = star.to_grams();
let earth_masses = planet.to_earth_masses();
```

### Time Operations

```rust
use units::Time;

// Create time spans
let orbital_period = Time::from_days(365.25);
let disk_lifetime = Time::from_myr(3.0);  // 3 million years

// Convert and combine
let years = orbital_period.to_years();
let total = disk_lifetime + Time::from_years(1e6);
```

### Temperature

```rust
use units::Temperature;

// Create temperatures
let disk_temp = Temperature::from_kelvin(280.0);
let earth_surface = Temperature::from_celsius(15.0);

// Convert between units
let celsius = disk_temp.to_celsius();
let fahrenheit = disk_temp.to_fahrenheit();

// Physical constants
let freezing = Temperature::water_freezing();  // 273.15 K
let boiling = Temperature::water_boiling();    // 373.15 K
```

### Velocity

```rust
use units::{Velocity, circular_orbital_velocity};

// Create velocities
let v = Velocity::from_au_per_year(6.28);  // ~Earth's orbital velocity
let cm_s = v.to_cm_per_sec();

// Calculate circular orbital velocity
let earth_v = circular_orbital_velocity(1.0, 1.0);  // M=1 M☉, r=1 AU
```

### Mass Rate (Accretion Fluxes)

```rust
use units::{MassRate, Time};

// Pebble accretion flux in planet formation
let pebble_flux = MassRate::from_earth_masses_per_myr(60.0);

// Integrate over time to get total accreted mass
let duration = Time::from_myr(1.0);
let accreted_mass = pebble_flux.integrate(duration);
// Result: ~60 Earth masses

// Stellar wind mass loss
let solar_wind = MassRate::from_solar_masses_per_year(1e-14);
```

### Surface Density (Disk Profiles)

```rust
use units::SurfaceDensity;

// Minimum Mass Solar Nebula at 1 AU
let sigma_1au = SurfaceDensity::from_grams_per_cm2(1700.0);

// Apply power-law scaling: Σ(r) = Σ₀ (r/r₀)^(-1.5)
let sigma_5au = sigma_1au.power_law_scaling(5.0, 1.5);
// Result: ~152 g/cm² at Jupiter's distance

// Convert to SI units
let kg_per_m2 = sigma_1au.to_kg_per_m2();  // 17000 kg/m²
```

### Volume Density (Planetary Composition)

```rust
use units::VolumeDensity;

// Earth's bulk density
let earth_density = VolumeDensity::from_grams_per_cm3(5.5);

// Use material constants
let iron_core = VolumeDensity::iron();           // 7.9 g/cm³
let rock_mantle = VolumeDensity::silicate_rock(); // 3.3 g/cm³
let ice = VolumeDensity::water_ice();            // 0.92 g/cm³

// Calculate weighted average for Earth-like planet
let components = vec![
    (VolumeDensity::iron(), 0.32),           // 32% iron core
    (VolumeDensity::silicate_rock(), 0.68),  // 68% rocky mantle
];
let bulk_density = VolumeDensity::weighted_average(&components).unwrap();
// Result: ~5.1 g/cm³ (uncompressed)
```

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
units = { path = "crates/units" }
```

## Requirements

- Rust 1.56 or later
- The `approx` crate for floating-point comparisons in tests
- The `serde` crate for serialization support
