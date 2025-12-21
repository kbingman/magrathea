# DiskModel Trait & GridDisk Implementation Plan

## Overview

This plan introduces a trait abstraction for protoplanetary disk models, enabling both fast analytical disks (for statistical generation) and evolving grid-based disks (for emergent simulation) to share physics calculations.

**Goals:**
1. Abstract common disk physics into a reusable trait
2. Preserve existing `GasDisk` functionality unchanged
3. Enable future `GridDisk` with viscous evolution
4. Maintain idiomatic Rust: traits, pattern matching, minimal mutation

**Estimated Effort:** 4-6 focused PRs over 1-2 weeks

---

## Part 1: DiskModel Trait

### 1.1 Design Rationale

The trait separates "what defines a disk" from "what we can compute from it":

| Category | Methods | Rationale |
|----------|---------|-----------|
| **Required** | `surface_density`, `temperature`, `stellar_mass`, `alpha`, `inner_radius`, `outer_radius` | Minimal set that varies between representations |
| **Override-able** | `pressure_gradient_log` | Numerical default works for any disk; analytical override for power-law |
| **Derived** | Everything else (15+ methods) | Physics identical regardless of representation |

This means a new disk type only needs to implement 6 methods to get full functionality.

### 1.2 Trait Definition

```rust
pub trait DiskModel {
    // === Required ===
    fn surface_density(&self, r: Length) -> SurfaceDensity;
    fn temperature(&self, r: Length) -> Temperature;
    fn stellar_mass(&self) -> Mass;
    fn alpha(&self) -> f64;
    fn inner_radius(&self) -> Length;
    fn outer_radius(&self) -> Length;

    // === Override-able (has default) ===
    fn pressure_gradient_log(&self, r: Length) -> f64 {
        // Numerical differentiation via central difference
    }

    // === Derived (default implementations) ===
    fn orbital_frequency(&self, r: Length) -> AngularVelocity { /* Ω_K = √(GM/r³) */ }
    fn keplerian_velocity(&self, r: Length) -> Velocity { /* v_K = √(GM/r) */ }
    fn orbital_period(&self, r: Length) -> Time { /* 2π/Ω */ }
    fn sound_speed(&self, r: Length) -> Velocity { /* √(kT/μm_p) */ }
    fn thermal_velocity(&self, r: Length) -> Velocity { /* √(8/π) × c_s */ }
    fn scale_height(&self, r: Length) -> Length { /* c_s / Ω */ }
    fn aspect_ratio(&self, r: Length) -> f64 { /* h/r */ }
    fn midplane_density(&self, r: Length) -> Density { /* Σ / (√2π h) */ }
    fn pressure(&self, r: Length) -> Pressure { /* ρ c_s² */ }
    fn pressure_gradient_parameter(&self, r: Length) -> f64 { /* η */ }
    fn sub_keplerian_velocity(&self, r: Length) -> Velocity { /* η v_K */ }
    fn viscosity(&self, r: Length) -> f64 { /* α c_s h */ }
    fn viscous_timescale(&self, r: Length) -> Time { /* r²/ν */ }
    fn mean_free_path(&self, r: Length) -> Length { /* μ m_p / (σ ρ) */ }
    fn is_valid_radius(&self, r: Length) -> bool { /* bounds check */ }
}
```

### 1.3 Separate Mass Trait

Total mass calculation requires knowing the profile form:
- Power-law: closed-form integral
- Tabulated: numerical integration

```rust
pub trait DiskMass: DiskModel {
    fn total_mass(&self) -> Mass;
}
```

### 1.4 Implementation for GasDisk

```rust
impl DiskModel for GasDisk {
    fn surface_density(&self, r: Length) -> SurfaceDensity {
        // Existing implementation unchanged
        let ratio = r.to_cm() / self.r_0.to_cm();
        SurfaceDensity::from_grams_per_cm2(
            self.sigma_0.to_grams_per_cm2() * ratio.powf(-self.sigma_exponent)
        )
    }

    fn temperature(&self, r: Length) -> Temperature {
        // Existing implementation unchanged
    }

    fn stellar_mass(&self) -> Mass { self.stellar_mass }
    fn alpha(&self) -> f64 { self.alpha }
    fn inner_radius(&self) -> Length { self.inner_radius }
    fn outer_radius(&self) -> Length { self.outer_radius }

    // OVERRIDE: Analytical expression for power-law
    fn pressure_gradient_log(&self, _r: Length) -> f64 {
        // d ln P / d ln r = -(p + (3+q)/2)
        -(self.sigma_exponent + (3.0 + self.temp_exponent) / 2.0)
    }
}

impl DiskMass for GasDisk {
    fn total_mass(&self) -> Mass {
        // Existing analytical integration
    }
}
```

### 1.5 PR Sequence for Part 1

#### PR 1.1: Trait Definition
**Files:** `src/disk/disk_model.rs`

**Contents:**
- `DiskModel` trait with all method signatures
- Default implementations for derived quantities
- `DiskMass` trait
- Unit tests using a minimal `TestDisk` struct

**Test:**
```rust
#[test]
fn trait_defaults_produce_reasonable_values() {
    let disk = TestDisk::mmsn();
    let eta = disk.pressure_gradient_parameter(Length::from_au(1.0));
    assert!(eta > 0.001 && eta < 0.01);
}
```

#### PR 1.2: GasDisk Implementation
**Files:** `src/disk/gas_disk.rs` (modify)

**Contents:**
- `impl DiskModel for GasDisk`
- `impl DiskMass for GasDisk`
- Analytical `pressure_gradient_log` override

**Test:**
```rust
#[test]
fn gas_disk_trait_matches_direct_methods() {
    let disk = GasDisk::mmsn();
    
    // Trait method vs existing method should be identical
    let eta_trait = DiskModel::pressure_gradient_parameter(&disk, Length::from_au(1.0));
    let eta_direct = disk.pressure_gradient_parameter(Length::from_au(1.0));
    
    assert!((eta_trait - eta_direct).abs() < 1e-15);
}
```

#### PR 1.3: Deprecation & Cleanup (Optional)
**Decision Point:** Keep both trait methods and inherent methods, or migrate callers to trait?

**Option A: Keep Both**
- Trait for generic code (`fn analyze<D: DiskModel>(disk: &D)`)
- Inherent methods for direct usage (`disk.pressure(r)`)
- No breaking changes

**Option B: Migrate to Trait**
- Remove inherent methods that duplicate trait
- All callers use trait methods
- Cleaner API, but requires `use DiskModel` at call sites

**Recommendation:** Option A initially. Revisit after GridDisk exists.

---

## Part 2: GridDisk Implementation

### 2.1 Design Overview

`GridDisk` stores surface density on a discrete radial grid and supports:
- Interpolated property queries
- Viscous evolution via diffusion equation
- Local modifications (gap opening, accretion)

```rust
pub struct GridDisk {
    /// Radial grid points (logarithmically spaced)
    radii: Vec<f64>,          // cm
    
    /// Surface density at each grid point
    sigma: Vec<f64>,          // g/cm²
    
    /// Temperature profile (analytical for now, could be gridded later)
    temp_0: f64,              // K at r_0
    temp_exponent: f64,       // T ∝ r^(-q)
    r_0: f64,                 // cm, reference radius
    
    /// Disk parameters
    stellar_mass: f64,        // g
    alpha: f64,
}
```

### 2.2 Grid Design Decisions

**Logarithmic Spacing:**
Disk properties vary over orders of magnitude in radius. Logarithmic spacing gives uniform resolution in `ln(r)`:

```rust
fn log_spaced_grid(r_min: f64, r_max: f64, n: usize) -> Vec<f64> {
    let log_min = r_min.ln();
    let log_max = r_max.ln();
    (0..n)
        .map(|i| {
            let frac = i as f64 / (n - 1) as f64;
            (log_min + frac * (log_max - log_min)).exp()
        })
        .collect()
}
```

**Typical Grid Size:**
- Minimum viable: 50 points (fast but coarse)
- Recommended: 100-200 points (good balance)
- High fidelity: 500+ points (for validation)

**Interpolation:**
Log-linear interpolation (linear in log-log space) matches power-law profiles exactly:

```rust
fn interpolate_log(&self, r: f64) -> f64 {
    let log_r = r.ln();
    let log_radii: Vec<f64> = self.radii.iter().map(|&x| x.ln()).collect();
    
    // Find bracketing indices
    let i = log_radii.partition_point(|&x| x < log_r).saturating_sub(1);
    let i = i.min(self.radii.len() - 2);
    
    // Linear interpolation in log-log space
    let t = (log_r - log_radii[i]) / (log_radii[i + 1] - log_radii[i]);
    let log_sigma_i = self.sigma[i].ln();
    let log_sigma_ip1 = self.sigma[i + 1].ln();
    
    (log_sigma_i + t * (log_sigma_ip1 - log_sigma_i)).exp()
}
```

### 2.3 Viscous Evolution

The surface density evolves according to:

```
∂Σ/∂t = (3/r) ∂/∂r [r^(1/2) ∂/∂r (ν Σ r^(1/2))]
```

**Discretization (explicit scheme):**

```rust
impl GridDisk {
    pub fn evolve_viscous(&mut self, dt: f64) {
        let n = self.radii.len();
        let mut d_sigma = vec![0.0; n];
        
        for i in 1..n-1 {
            let r = self.radii[i];
            let r_m = self.radii[i - 1];
            let r_p = self.radii[i + 1];
            
            let nu = self.viscosity_at(i);
            let nu_m = self.viscosity_at(i - 1);
            let nu_p = self.viscosity_at(i + 1);
            
            // Compute F = ν Σ r^(1/2) at cell faces
            let f_m = 0.5 * (nu_m * self.sigma[i-1] * r_m.sqrt() 
                          + nu * self.sigma[i] * r.sqrt());
            let f_p = 0.5 * (nu * self.sigma[i] * r.sqrt() 
                          + nu_p * self.sigma[i+1] * r_p.sqrt());
            
            // Compute ∂F/∂r at cell faces
            let df_dr_m = (nu * self.sigma[i] * r.sqrt() 
                         - nu_m * self.sigma[i-1] * r_m.sqrt()) / (r - r_m);
            let df_dr_p = (nu_p * self.sigma[i+1] * r_p.sqrt() 
                         - nu * self.sigma[i] * r.sqrt()) / (r_p - r);
            
            // ∂Σ/∂t = (3/r) ∂/∂r [r^(1/2) ∂F/∂r]
            let r_face_m = (r * r_m).sqrt();
            let r_face_p = (r * r_p).sqrt();
            
            let term_m = r_face_m.sqrt() * df_dr_m;
            let term_p = r_face_p.sqrt() * df_dr_p;
            
            d_sigma[i] = (3.0 / r) * (term_p - term_m) / (r_face_p - r_face_m);
        }
        
        // Apply update
        for i in 1..n-1 {
            self.sigma[i] += d_sigma[i] * dt;
            self.sigma[i] = self.sigma[i].max(1e-30); // floor to avoid negatives
        }
        
        // Boundary conditions: zero-torque at inner, zero-flux at outer
        self.sigma[0] = self.sigma[1];
        self.sigma[n-1] = self.sigma[n-2];
    }
    
    fn viscosity_at(&self, i: usize) -> f64 {
        let r = self.radii[i];
        let c_s = self.sound_speed_at(i);
        let h = c_s / self.orbital_frequency_at(i);
        self.alpha * c_s * h
    }
}
```

**Timestep Constraint (CFL condition):**

```rust
fn max_timestep(&self) -> f64 {
    self.radii.windows(2)
        .enumerate()
        .map(|(i, w)| {
            let dr = w[1] - w[0];
            let nu = self.viscosity_at(i);
            0.5 * dr * dr / nu  // CFL factor of 0.5
        })
        .fold(f64::INFINITY, f64::min)
}
```

### 2.4 Validation Strategy

**Test 1: Self-Similar Solution**

A disk with `Σ ∝ r^(-1)` evolving under constant α has an analytical self-similar solution. The surface density profile maintains its shape while spreading:

```rust
#[test]
fn viscous_evolution_self_similar() {
    let mut disk = GridDisk::power_law(
        stellar_mass: Mass::solar(),
        sigma_0: 1700.0,
        sigma_exponent: 1.0,
        n_radii: 200,
    );
    
    let t_visc_1au = disk.viscous_timescale(Length::from_au(1.0));
    
    // Evolve for 0.1 viscous times
    let dt = disk.max_timestep();
    let n_steps = (0.1 * t_visc_1au.to_seconds() / dt) as usize;
    
    for _ in 0..n_steps {
        disk.evolve_viscous(dt);
    }
    
    // Check that inner disk has depleted, outer disk has gained mass
    let sigma_0p5au_initial = 1700.0 / 0.5;  // Power-law at 0.5 AU
    let sigma_0p5au_final = disk.surface_density(Length::from_au(0.5));
    
    assert!(sigma_0p5au_final.to_grams_per_cm2() < sigma_0p5au_initial);
}
```

**Test 2: Mass Conservation**

```rust
#[test]
fn viscous_evolution_conserves_mass() {
    let mut disk = GridDisk::mmsn(n_radii: 200);
    
    let mass_initial = disk.total_mass();
    
    let dt = disk.max_timestep();
    for _ in 0..1000 {
        disk.evolve_viscous(dt);
    }
    
    let mass_final = disk.total_mass();
    let mass_lost = (mass_initial - mass_final) / mass_initial;
    
    // Should lose < 1% (boundary losses only)
    assert!(mass_lost.abs() < 0.01, "Lost {:.2}% of mass", mass_lost * 100.0);
}
```

**Test 3: Comparison with Armitage Figure 3.3**

The textbook has a specific figure showing viscous spreading. We can reproduce it:

```rust
#[test]
fn matches_armitage_figure_3_3() {
    // Initial ring at 1 AU
    let mut disk = GridDisk::ring(
        r_ring: Length::from_au(1.0),
        sigma_ring: 1000.0,  // g/cm²
        width: Length::from_au(0.1),
        n_radii: 500,
    );
    
    // Evolve to t = 0.1, 0.3, 1.0 viscous times
    // Compare profile shapes to textbook figure
}
```

### 2.5 PR Sequence for Part 2

#### PR 2.1: GridDisk Structure & Construction
**Files:** `src/disk/grid_disk.rs`

**Contents:**
- `GridDisk` struct definition
- Constructors: `new()`, `from_power_law()`, `from_gas_disk()`
- Logarithmic grid generation
- Log-linear interpolation

**Test:**
```rust
#[test]
fn from_gas_disk_matches_original() {
    let power_law = GasDisk::mmsn();
    let grid = GridDisk::from_gas_disk(&power_law, 200);
    
    for r_au in [0.5, 1.0, 5.0, 10.0, 50.0] {
        let r = Length::from_au(r_au);
        let sigma_pl = power_law.surface_density(r);
        let sigma_grid = grid.surface_density(r);
        
        let rel_err = (sigma_pl.to_grams_per_cm2() - sigma_grid.to_grams_per_cm2()).abs()
                    / sigma_pl.to_grams_per_cm2();
        assert!(rel_err < 0.01, "Error at {} AU: {:.2}%", r_au, rel_err * 100.0);
    }
}
```

#### PR 2.2: DiskModel Implementation for GridDisk
**Files:** `src/disk/grid_disk.rs`

**Contents:**
- `impl DiskModel for GridDisk`
- Uses numerical `pressure_gradient_log` default
- `impl DiskMass for GridDisk` with trapezoidal integration

**Test:**
```rust
#[test]
fn grid_disk_eta_reasonable() {
    let disk = GridDisk::from_gas_disk(&GasDisk::mmsn(), 200);
    let eta = disk.pressure_gradient_parameter(Length::from_au(1.0));
    assert!(eta > 0.001 && eta < 0.01);
}
```

#### PR 2.3: Viscous Evolution (Explicit)
**Files:** `src/disk/grid_disk.rs`

**Contents:**
- `evolve_viscous(&mut self, dt: f64)`
- `max_timestep(&self) -> f64`
- CFL condition

**Test:**
```rust
#[test]
fn viscous_evolution_depletes_inner_disk() {
    // ... as shown in validation strategy
}
```

#### PR 2.4: Mass Conservation & Boundary Conditions
**Files:** `src/disk/grid_disk.rs`

**Contents:**
- Proper boundary condition handling
- Mass accounting methods
- Angular momentum tracking (optional)

**Test:**
```rust
#[test]
fn mass_conservation_over_long_evolution() {
    // ... as shown in validation strategy
}
```

#### PR 2.5: Photoevaporation (Optional Enhancement)
**Files:** `src/disk/grid_disk.rs`

**Contents:**
- `PhotoevaporationModel` enum
- `apply_photoevaporation(&mut self, dt: f64)`
- Disk lifetime behavior

**Test:**
```rust
#[test]
fn photoevaporation_clears_disk_in_few_myr() {
    let mut disk = GridDisk::mmsn(200);
    disk.set_photoevaporation(PhotoevaporationModel::EUV { 
        rate: 1e-9  // M_sun / yr
    });
    
    // Evolve for 10 Myr
    while disk.age() < Time::from_myr(10.0) {
        disk.step();
    }
    
    assert!(disk.total_mass() < Mass::from_solar_masses(1e-4));
}
```

---

## Part 3: Integration & Usage

### 3.1 Generic Functions

With the trait in place, physics code can be generic:

```rust
/// Calculate stopping time for a particle in any disk
fn stopping_time<D: DiskModel>(
    disk: &D, 
    r: Length, 
    particle_radius: Length,
    particle_density: Density,
) -> Time {
    let rho_gas = disk.midplane_density(r);
    let c_s = disk.sound_speed(r);
    let mfp = disk.mean_free_path(r);
    
    // Epstein vs Stokes regime...
}

/// Calculate radial drift velocity for a particle
fn radial_drift<D: DiskModel>(
    disk: &D,
    r: Length,
    stokes_number: f64,
) -> Velocity {
    let eta = disk.pressure_gradient_parameter(r);
    let v_k = disk.keplerian_velocity(r);
    
    // Weidenschilling (1977) formula
    let v_r = -2.0 * eta * v_k * stokes_number / (1.0 + stokes_number.powi(2));
    Velocity::from_cm_per_sec(v_r.to_cm_per_sec())
}
```

### 3.2 Simulation Context

```rust
pub struct EmergentSimulation<D: DiskModel> {
    disk: D,
    particles: Vec<ParticleBin>,
    bodies: Vec<DiscreteBody>,
    time: Time,
}

// For evolving simulations
impl EmergentSimulation<GridDisk> {
    fn step(&mut self) {
        let dt = self.adaptive_timestep();
        self.disk.evolve_viscous(dt.to_seconds());
        // ... rest of simulation
    }
}

// For static snapshots (fast generation)
impl EmergentSimulation<GasDisk> {
    fn evaluate_at_time(&self, t: Time) -> SystemSnapshot {
        // No disk evolution, just evaluate planet formation
    }
}
```

### 3.3 Module Structure

```
src/disk/
├── mod.rs
├── constants.rs          # Physical constants
├── disk_model.rs         # DiskModel trait + DiskMass trait
├── gas_disk.rs           # Power-law disk (existing, + trait impl)
├── gas_disk_test.rs      # Existing tests
├── grid_disk.rs          # Tabulated disk with evolution
└── grid_disk_test.rs     # Grid disk tests
```

---

## Part 4: Timeline & Milestones

### Week 1: Trait Foundation

| Day | PR | Deliverable |
|-----|-----|-------------|
| 1-2 | 1.1 | `DiskModel` trait with defaults, unit tests |
| 3-4 | 1.2 | `impl DiskModel for GasDisk`, verify existing tests pass |
| 5 | — | Review, cleanup, merge |

**Milestone Demo:** Print table comparing trait methods vs direct methods for MMSN disk at various radii.

### Week 2: GridDisk Core

| Day | PR | Deliverable |
|-----|-----|-------------|
| 1-2 | 2.1 | `GridDisk` construction, interpolation |
| 3-4 | 2.2 | `impl DiskModel for GridDisk` |
| 5 | — | Review, cleanup, merge |

**Milestone Demo:** Plot `GridDisk` Σ(r) overlaid on `GasDisk` Σ(r), showing interpolation accuracy.

### Week 3: Viscous Evolution

| Day | PR | Deliverable |
|-----|-----|-------------|
| 1-2 | 2.3 | `evolve_viscous()`, explicit scheme |
| 3-4 | 2.4 | Mass conservation, boundary conditions |
| 5 | — | Review, cleanup, merge |

**Milestone Demo:** Animate disk surface density evolution over 1 Myr (matches Armitage Fig 3.3).

### Week 4 (Optional): Enhancements

| Day | PR | Deliverable |
|-----|-----|-------------|
| 1-2 | 2.5 | Photoevaporation |
| 3-4 | — | Implicit solver for stiff evolution (if needed) |
| 5 | — | Documentation, examples |

**Milestone Demo:** Plot disk mass vs time showing characteristic 1-10 Myr disk lifetime.

---

## Appendix A: Physical Constants

For reference, the constants used in derived quantities:

| Constant | Symbol | Value | Units |
|----------|--------|-------|-------|
| Gravitational constant | G | 6.674×10⁻⁸ | cm³/(g·s²) |
| Boltzmann constant | k_B | 1.381×10⁻¹⁶ | erg/K |
| Proton mass | m_p | 1.673×10⁻²⁴ | g |
| Mean molecular weight | μ | 2.34 | — |
| H₂ cross section | σ_mol | 2×10⁻¹⁵ | cm² |

---

## Appendix B: Key Equations

### Derived Quantities

| Quantity | Symbol | Formula |
|----------|--------|---------|
| Keplerian frequency | Ω_K | √(GM_*/r³) |
| Keplerian velocity | v_K | √(GM_*/r) |
| Sound speed | c_s | √(k_B T / μ m_p) |
| Scale height | h | c_s / Ω_K |
| Aspect ratio | h/r | c_s / v_K |
| Midplane density | ρ | Σ / (√2π h) |
| Pressure | P | ρ c_s² |
| Pressure gradient | η | -(h/r)² × ½ × d ln P / d ln r |
| Viscosity | ν | α c_s h |

### Power-Law Profiles

For Σ ∝ r^(-p) and T ∝ r^(-q):

| Quantity | Radial Dependence |
|----------|-------------------|
| c_s | r^(-q/2) |
| h | r^((3-q)/2) |
| h/r | r^((1-q)/2) |
| ρ | r^(-(p + (3-q)/2)) |
| P | r^(-(p + (3+q)/2)) |
| d ln P / d ln r | -(p + (3+q)/2) = const |
| η | r^(1-q) |

For MMSN (p=1, q=0.5): d ln P / d ln r = -2.75

---

## Appendix C: Alternative Approaches Considered

### Why Not Enum-Based Profiles?

```rust
// Rejected approach
enum SigmaProfile {
    PowerLaw { ... },
    Tabulated { ... },
}
```

**Problems:**
- Every method has a match statement
- `evolve()` only valid for `Tabulated` — awkward API
- Mixing mutable/immutable semantics

### Why Not Just GridDisk?

Could skip the trait and always use `GridDisk`, initializing from power-law.

**Problems:**
- Unnecessary overhead for statistical generation
- Interpolation errors compound over many queries
- Power-law evaluation is trivial; why add complexity?

The trait approach lets each representation optimize for its use case.
