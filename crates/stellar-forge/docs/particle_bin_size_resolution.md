# ParticleBin Size Resolution

## The Problem

The current `ParticleBin` implementation uses mean particle size as a representative for all aerodynamic calculations:

```rust
// Current implementation in particle_bin.rs
pub fn mean_stopping_time<D: DiskModel>(&self, disk: &D) -> Time {
    // For now, use the mean size as representative
    // A more accurate implementation would integrate over the distribution
    let mean_s = self.mean_size();
    let particle = Particle::new(mean_s, rho_m);
    particle.stopping_time(disk, r)
}
```

This loses critical physics. Different-sized particles have different Stokes numbers, leading to different:

| Property | Consequence |
|----------|-------------|
| Radial drift velocity | τ_s = 1 particles drift fastest (~50 m/s) |
| Settling rate | Larger particles settle to thinner midplane layers |
| Turbulent coupling | Different sizes couple to different eddy scales |

For coagulation, the *relative velocity* between size bins drives collision rates. Homogenizing to a mean erases this differential motion entirely.

---

## Option A: Size-Resolved Queries

Keep the current structure but add methods that operate over the full size distribution.

```rust
impl ParticleBin {
    /// Drift velocity for each size in the distribution.
    pub fn drift_velocities<D: DiskModel>(&self, disk: &D) -> Vec<(Length, Velocity)> {
        let r = self.radial_center();
        let rho_m = self.material_density();
        
        self.size_distribution
            .size_bins()
            .map(|s| {
                let particle = Particle::new(s, rho_m);
                (s, particle.radial_drift_velocity(disk, r))
            })
            .collect()
    }
    
    /// Relative velocity between particles of sizes s1 and s2.
    pub fn relative_velocity<D: DiskModel>(
        &self, 
        disk: &D, 
        s1: Length, 
        s2: Length,
    ) -> Velocity {
        let r = self.radial_center();
        let rho_m = self.material_density();
        
        let p1 = Particle::new(s1, rho_m);
        let p2 = Particle::new(s2, rho_m);
        
        // Differential radial drift
        let dv_r = (p1.radial_drift_velocity(disk, r).to_cm_per_sec()
                  - p2.radial_drift_velocity(disk, r).to_cm_per_sec()).abs();
        
        // Differential azimuthal drift  
        let dv_phi = (p1.azimuthal_drift_velocity(disk, r).to_cm_per_sec()
                    - p2.azimuthal_drift_velocity(disk, r).to_cm_per_sec()).abs();
        
        // Turbulent relative velocity (Ormel & Cuzzi 2007)
        let dv_turb = self.turbulent_relative_velocity(disk, s1, s2);
        
        // Quadrature sum
        Velocity::from_cm_per_sec(
            (dv_r.powi(2) + dv_phi.powi(2) + dv_turb.powi(2)).sqrt()
        )
    }
}
```

Requires `SizeDistribution` to expose its size bins:

```rust
impl SizeDistribution {
    /// Iterator over representative sizes.
    pub fn size_bins(&self) -> impl Iterator<Item = Length> + '_ {
        match self {
            Self::Monodisperse { size, .. } => vec![*size].into_iter(),
            Self::Binned { bin_edges, .. } => {
                bin_edges.windows(2)
                    .map(|w| Length::from_cm((w[0] * w[1]).sqrt()))
                    .collect::<Vec<_>>()
                    .into_iter()
            }
            Self::PowerLaw { s_min, s_max, .. } => {
                log_spaced(*s_min, *s_max, 30).into_iter()
            }
        }
    }
}
```

**Pros:**
- Minimal structural change
- Preserves existing API

**Cons:**
- PowerLaw case requires sampling or lazy conversion
- Inconsistent iteration patterns across variants

---

## Option B: Always Use Binned Representation

Convert to `Binned` on construction, making size-resolved operations uniform.

```rust
impl ParticleBin {
    pub fn from_disk<D: DiskModel>(disk: &D, r: Length, width: Length) -> Self {
        // ... existing code ...
        
        let size_dist = SizeDistribution::mrn(
            Mass::from_grams(total_mass),
            Density::from_grams_per_cm3(SILICATE_DENSITY),
        ).to_binned(30);  // Always convert
        
        // ...
    }
}
```

Iteration becomes clean and uniform:

```rust
impl ParticleBin {
    pub fn drift_velocities<D: DiskModel>(
        &self, 
        disk: &D,
    ) -> impl Iterator<Item = (Length, Velocity)> + '_ {
        let r = self.radial_center();
        let rho_m = self.material_density();
        
        self.size_distribution
            .bins()  // Always available
            .iter()
            .map(move |(size, _mass)| {
                let particle = Particle::new(*size, rho_m);
                (*size, particle.radial_drift_velocity(disk, r))
            })
    }
}
```

**Pros:**
- Uniform representation across all cases
- Clean iteration with no match arms
- Mass already tracked per bin
- `to_binned()` preserves mass exactly

**Cons:**
- Loses analytical form
- Slightly more memory (bin arrays vs. parameters)

---

## Option C: Explicit 2D Structure

Separate radial and size dimensions into distinct types:

```rust
/// A radial zone containing size-resolved particle populations.
pub struct RadialZone {
    radial_center: f64,
    radial_width: f64,
    size_bins: Vec<SizeBin>,
}

/// A single size bin within a radial zone.
pub struct SizeBin {
    size: f64,
    mass: f64,
    scale_height: f64,  // Can differ per size (settling)
    velocity_dispersion: f64,
}
```

**Pros:**
- Conceptually cleanest separation
- Per-size-bin scale heights (important for differential settling)
- Natural structure for coagulation matrix operations

**Cons:**
- Larger refactor
- May be overkill if per-size-bin scale heights aren't needed yet

---

## Recommendation

Follow the incremental development approach:

### Step 1: Convert to Binned (Now)

Modify `ParticleBin::from_disk` to always use binned representation:

```rust
let size_dist = SizeDistribution::mrn(...).to_binned(30);
```

This is approximately 5 lines changed. The `to_binned()` method already preserves mass exactly, so no fidelity is lost.

### Step 2: Add Size-Resolved Methods (Next PR)

Add `drift_velocities()` and `relative_velocity()` methods that iterate over the binned distribution.

### Step 3: Refactor to Option C (When Needed)

If coagulation implementation reveals that per-size-bin scale heights are necessary (likely, since large grains settle faster), refactor to the explicit 2D structure at that point.

---

## Rationale

The key insight is that `to_binned()` is already implemented and conserves mass exactly. We're not losing physical fidelity—only representation flexibility that isn't being used.

For the coagulation physics ahead (collision kernels, fragmentation thresholds), binned representation is what the algorithms naturally operate on. Converting early simplifies all downstream code.

The explicit 2D structure (Option C) is probably the right long-term architecture, but following the "no forward references" rule, we defer that refactor until the physics demands it.
