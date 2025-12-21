# ParticleBin Implementation Plan

## Overview

This document outlines the implementation of `SizeDistribution` and `ParticleBin` - the statistical representation of particle populations in protoplanetary disks.

**Goal:** Efficiently represent large populations of solid particles (dust, pebbles, boulders) without tracking individuals.

**Dependencies:**
- `Particle` struct (completed)
- `DiskModel` trait (completed)

---

## Design

### Why Statistical Representation?

A protoplanetary disk contains ~10^15+ dust grains. We can't track them individually. Instead, we represent particles statistically:

- **Size distribution**: How many particles of each size?
- **Spatial binning**: Where are they in the disk?
- **Bulk properties**: Total mass, velocity dispersion, scale height

This enables:
- Efficient computation (O(n_bins) not O(n_particles))
- Coagulation/fragmentation via Smoluchowski equation
- Transition to discrete bodies when gravitational instability triggers

### Key Types

```
┌─────────────────────────────────────────────────────┐
│ ParticleBin                                         │
│                                                     │
│  radial_center: Length     (where in disk)          │
│  radial_width: Length      (annulus width)          │
│  size_distribution: SizeDistribution                │
│  material_density: Density (internal ρ_m)           │
│  surface_density: SurfaceDensity (total Σ_p)        │
│  scale_height: Length      (vertical extent)        │
│                                                     │
└─────────────────────────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────┐
│ SizeDistribution                                    │
│                                                     │
│  PowerLaw { s_min, s_max, exponent, total_mass }    │
│  Binned { bin_edges, mass_per_bin }                 │
│  Monodisperse { size, total_mass }                  │
│                                                     │
└─────────────────────────────────────────────────────┘
```

---

## Part 1: SizeDistribution

### 1.1 Enum Variants

```rust
pub enum SizeDistribution {
    /// Power-law: dn/ds ∝ s^(-q)
    /// MRN distribution has q = 3.5
    PowerLaw {
        s_min: f64,        // cm, minimum size
        s_max: f64,        // cm, maximum size
        exponent: f64,     // q in s^(-q)
        total_mass: f64,   // g, integrated mass
        material_density: f64,  // g/cm³
    },

    /// Discretized into logarithmic mass bins
    /// More flexible, required for coagulation
    Binned {
        bin_edges: Vec<f64>,    // cm, N+1 edges for N bins
        mass_per_bin: Vec<f64>, // g, mass in each bin
        material_density: f64,  // g/cm³
    },

    /// Single size (delta function)
    /// Useful for testing and simple models
    Monodisperse {
        size: f64,         // cm
        total_mass: f64,   // g
        material_density: f64,  // g/cm³
    },
}
```

### 1.2 Core Methods

```rust
impl SizeDistribution {
    // === Constructors ===

    /// MRN-like power law (q=3.5) from μm to mm
    fn mrn(total_mass: f64, material_density: f64) -> Self;

    /// Power law with custom parameters
    fn power_law(s_min: f64, s_max: f64, exponent: f64,
                 total_mass: f64, material_density: f64) -> Self;

    /// Single size
    fn monodisperse(size: f64, total_mass: f64, material_density: f64) -> Self;

    /// Convert to binned representation (for coagulation)
    fn to_binned(&self, n_bins: usize) -> Self;

    // === Queries ===

    /// Total mass in the distribution
    fn total_mass(&self) -> f64;

    /// Mass-weighted mean size: ∫ s × m(s) ds / M_tot
    fn mean_size(&self) -> f64;

    /// Size at which half the mass is in smaller particles
    fn median_size(&self) -> f64;

    /// Maximum particle size
    fn max_size(&self) -> f64;

    /// Minimum particle size
    fn min_size(&self) -> f64;

    /// Number density at size s: dn/ds
    fn number_density_at(&self, s: f64) -> f64;

    /// Mass density at size s: dm/ds = (4π/3) ρ_m s³ dn/ds
    fn mass_density_at(&self, s: f64) -> f64;

    // === Weighted averages ===

    /// Compute ∫ f(s) × m(s) ds / M_tot for any function f
    fn mass_weighted_average<F: Fn(f64) -> f64>(&self, f: F) -> f64;
}
```

### 1.3 Physics Notes

**Power-law distributions:**
- MRN (Mathis-Rumpl-Nordsieck): q = 3.5, used for ISM dust
- Collisional equilibrium: q ≈ 3.5 (Dohnanyi 1969)
- Most mass at large sizes when q < 4
- Most surface area at small sizes when q > 3

**Mass in a power-law:**
For dn/ds ∝ s^(-q) with q ≠ 4:
```
M_tot = (4π/3) ρ_m ∫ s³ × s^(-q) ds
      = (4π/3) ρ_m × C × [s^(4-q)]_{s_min}^{s_max} / (4-q)
```

---

## Part 2: ParticleBin

### 2.1 Structure

```rust
pub struct ParticleBin {
    /// Center of radial annulus
    radial_center: f64,  // cm

    /// Width of radial annulus
    radial_width: f64,   // cm

    /// Size distribution of particles
    size_distribution: SizeDistribution,

    /// Total surface density of solids
    surface_density: f64,  // g/cm²

    /// Vertical scale height of particle layer
    scale_height: f64,  // cm

    /// Velocity dispersion (random velocities)
    velocity_dispersion: f64,  // cm/s
}
```

### 2.2 Core Methods

```rust
impl ParticleBin {
    // === Constructors ===

    /// Create from disk conditions at radius r
    /// Uses standard dust-to-gas ratio and MRN distribution
    fn from_disk<D: DiskModel>(disk: &D, r: Length, width: Length) -> Self;

    /// Create with custom parameters
    fn new(radial_center: Length, radial_width: Length,
           size_distribution: SizeDistribution,
           surface_density: SurfaceDensity,
           scale_height: Length) -> Self;

    // === Geometric queries ===

    /// Inner edge of annulus
    fn inner_radius(&self) -> Length;

    /// Outer edge of annulus
    fn outer_radius(&self) -> Length;

    /// Area of annulus: π(r_out² - r_in²)
    fn area(&self) -> f64;  // cm²

    // === Mass queries ===

    /// Total mass in this bin: Σ × Area
    fn total_mass(&self) -> Mass;

    /// Midplane density: Σ / (√2π × h)
    fn midplane_density(&self) -> Density;

    // === Aerodynamic properties ===

    /// Mean stopping time (mass-weighted over size distribution)
    fn mean_stopping_time<D: DiskModel>(&self, disk: &D) -> Time;

    /// Mean Stokes number
    fn mean_stokes_number<D: DiskModel>(&self, disk: &D) -> f64;

    /// Mean radial drift velocity
    fn mean_drift_velocity<D: DiskModel>(&self, disk: &D) -> Velocity;

    // === Evolution (future) ===

    /// Update scale height toward equilibrium
    fn relax_scale_height<D: DiskModel>(&mut self, disk: &D, dt: Time);

    /// Apply radial drift (move mass between bins)
    fn apply_drift(&mut self, dt: Time) -> MassFlux;
}
```

### 2.3 Initial Conditions

Standard assumptions for disk initialization:

| Parameter | Value | Notes |
|-----------|-------|-------|
| Dust-to-gas ratio | 0.01 | Solar metallicity |
| Initial size range | 0.1 μm - 1 μm | ISM grain sizes |
| Size exponent | 3.5 | MRN distribution |
| Material density | 3.0 g/cm³ | Silicate |
| Initial scale height | h_gas | Well-mixed with gas |

```rust
impl ParticleBin {
    fn from_disk<D: DiskModel>(disk: &D, r: Length, width: Length) -> Self {
        let dust_to_gas = 0.01;
        let sigma_gas = disk.surface_density(r);
        let sigma_dust = sigma_gas * dust_to_gas;

        let size_dist = SizeDistribution::mrn(
            sigma_dust.to_grams_per_cm2() * area,
            3.0,  // silicate density
        );

        Self {
            radial_center: r.to_cm(),
            radial_width: width.to_cm(),
            size_distribution: size_dist,
            surface_density: sigma_dust.to_grams_per_cm2(),
            scale_height: disk.scale_height(r).to_cm(),
            velocity_dispersion: disk.sound_speed(r).to_cm_per_sec() * 0.01,
        }
    }
}
```

---

## Part 3: Tests

### 3.1 SizeDistribution Tests

```rust
#[test]
fn power_law_mass_conservation() {
    let dist = SizeDistribution::power_law(
        1e-4, 0.1,  // 1 μm to 1 mm
        3.5,        // MRN exponent
        1000.0,     // 1 kg total
        3.0,        // silicate
    );
    assert_relative_eq!(dist.total_mass(), 1000.0, epsilon = 1e-10);
}

#[test]
fn to_binned_preserves_mass() {
    let power_law = SizeDistribution::mrn(1000.0, 3.0);
    let binned = power_law.to_binned(50);
    assert_relative_eq!(
        power_law.total_mass(),
        binned.total_mass(),
        epsilon = 0.01  // 1% tolerance
    );
}

#[test]
fn mean_size_between_bounds() {
    let dist = SizeDistribution::power_law(1e-4, 0.1, 3.5, 1000.0, 3.0);
    let mean = dist.mean_size();
    assert!(mean > 1e-4 && mean < 0.1);
}

#[test]
fn monodisperse_mean_equals_size() {
    let dist = SizeDistribution::monodisperse(0.01, 1000.0, 3.0);
    assert_relative_eq!(dist.mean_size(), 0.01, epsilon = 1e-10);
}
```

### 3.2 ParticleBin Tests

```rust
#[test]
fn from_disk_has_correct_dust_to_gas() {
    let disk = GasDisk::mmsn();
    let bin = ParticleBin::from_disk(&disk, Length::from_au(1.0), Length::from_au(0.1));

    let sigma_gas = disk.surface_density(Length::from_au(1.0));
    let sigma_dust = bin.surface_density();

    let ratio = sigma_dust.to_grams_per_cm2() / sigma_gas.to_grams_per_cm2();
    assert_relative_eq!(ratio, 0.01, epsilon = 1e-10);
}

#[test]
fn total_mass_equals_surface_density_times_area() {
    let disk = GasDisk::mmsn();
    let bin = ParticleBin::from_disk(&disk, Length::from_au(1.0), Length::from_au(0.1));

    let expected = bin.surface_density().to_grams_per_cm2() * bin.area();
    let actual = bin.total_mass().to_grams();

    assert_relative_eq!(actual, expected, epsilon = 1e-10);
}

#[test]
fn mean_stokes_number_reasonable() {
    let disk = GasDisk::mmsn();
    let bin = ParticleBin::from_disk(&disk, Length::from_au(1.0), Length::from_au(0.1));

    let tau = bin.mean_stokes_number(&disk);

    // μm grains should have τ << 1
    assert!(tau < 1e-3, "τ = {} too large for μm grains", tau);
}
```

---

## File Structure

```
src/particles/
├── mod.rs                 # Module exports
├── particle.rs            # Single particle (existing)
├── particle_test.rs       # Single particle tests (existing)
├── size_distribution.rs   # SizeDistribution enum
├── size_distribution_test.rs
├── particle_bin.rs        # ParticleBin struct
└── particle_bin_test.rs
```

---

## Implementation Order

### PR 1: SizeDistribution (~150 lines)

1. Create `size_distribution.rs` with enum and variants
2. Implement constructors: `power_law()`, `monodisperse()`, `mrn()`
3. Implement queries: `total_mass()`, `mean_size()`, `max_size()`, `min_size()`
4. Add tests for mass conservation and bounds
5. Defer `to_binned()` and `mass_weighted_average()` to PR 2

### PR 2: ParticleBin (~150 lines)

1. Create `particle_bin.rs` with struct
2. Implement `from_disk()` constructor
3. Implement geometric queries: `area()`, `inner_radius()`, `outer_radius()`
4. Implement mass queries: `total_mass()`, `midplane_density()`
5. Add tests for dust-to-gas ratio and mass calculations

### PR 3: Aerodynamic Integration (~100 lines)

1. Add `to_binned()` for SizeDistribution
2. Add `mass_weighted_average()`
3. Implement `mean_stopping_time()`, `mean_stokes_number()`, `mean_drift_velocity()`
4. Tests comparing to single-particle calculations

---

## Future Work (Not This PR)

- **Coagulation kernel**: Collision rates between size bins
- **Smoluchowski solver**: Evolve size distribution over time
- **Radial transport**: Move mass between spatial bins
- **Fragmentation**: Break large particles into small ones
- **Gravitational instability**: Trigger planetesimal formation

---

## References

- Birnstiel et al. (2010) - "A simple model for the evolution of the dust population in protoplanetary disks"
- Mathis, Rumpl & Nordsieck (1977) - "The size distribution of interstellar grains" (MRN)
- Dohnanyi (1969) - "Collisional model of asteroids and their debris"
