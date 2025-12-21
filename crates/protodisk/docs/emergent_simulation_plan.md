# Emergent Planet Formation Simulation

## Design Document & Development Plan

---

# Part 1: Design Document

## 1. Introduction & Philosophy

### 1.1 What "Emergent" Means Here

An emergent simulation differs fundamentally from both statistical models and hybrid approaches. Rather than prescribing outcomes (like a statistical model) or mixing rules with physics (like a hybrid), an emergent simulation defines **local rules** and lets global structure arise from their interactions.

The key insight from Armitage's textbook is that planet formation involves a cascade of coupled physical processes, each operating on different scales. Rather than modeling each stage separately and stitching them together, an emergent approach simulates the underlying interactions and lets the standard formation pathways (runaway growth, oligarchic growth, core accretion, disk instability) emerge naturally from the physics.

### 1.2 Design Goals

Given the stated priorities (believable exoplanets, fast generation, Rust idiomaticity), the simulation should:

1. **Produce diverse, plausible systems** without hand-tuning per-regime parameters
2. **Run quickly** by using simplified but physically-grounded local rules
3. **Exhibit emergent behaviors** like runaway growth, gap clearing, and resonant capture without explicitly coding them
4. **Maintain reproducibility** via seeded RNG infrastructure (ChaChaRng)

---

## 2. Core Architecture

### 2.1 Fundamental Entities

The simulation tracks three categories of entities, each defined by traits:

```rust
trait MassCarrier {
    fn mass(&self) -> Mass;
    fn position(&self) -> OrbitalPosition;  // (a, e, i) or (r, φ, z)
    fn velocity_dispersion(&self) -> Velocity;
}

trait GasInteracting {
    fn friction_time(&self, gas_density: Density) -> Duration;
    fn aerodynamic_regime(&self) -> DragRegime;  // Epstein | Stokes
}

trait GravitationallyInteracting {
    fn hill_radius(&self, stellar_mass: Mass) -> Length;
    fn escape_velocity(&self) -> Velocity;
    fn gravitational_focusing_factor(&self, approach_velocity: Velocity) -> f64;
}
```

**Entity Types:**

| Entity | Represented As | Key Properties |
|--------|----------------|----------------|
| Gas Disk | Continuous field (radial bins) | Σ(r), T(r), h/r, α |
| Particle Population | Statistical bin (per radial zone) | Size distribution, surface density, velocity dispersion |
| Discrete Bodies | Individual objects (N ≤ ~1000) | Mass, orbital elements, composition |
| Protoplanet | Discrete body + envelope | Core mass, envelope mass, accretion state |

The critical insight: **transition between representations is part of the physics**. A particle population "spawns" a discrete body when local conditions trigger gravitational instability. A discrete body "absorbs" nearby particle populations through accretion.

### 2.2 State Representation

```rust
struct SimulationState {
    time: Duration,
    star: StellarObject,
    gas_disk: GasDisk,
    particle_bins: Vec<ParticleBin>,     // Radial bins of size-distributed solids
    discrete_bodies: Vec<DiscreteBody>,   // Protoplanets, embryos
    rng: ChaChaRng,
}

struct GasDisk {
    surface_density: RadialProfile<SurfaceDensity>,
    temperature: RadialProfile<Temperature>,
    aspect_ratio: RadialProfile<f64>,
    alpha: f64,  // Shakura-Sunyaev viscosity parameter
    photoevaporation_rate: MassRate,
}

struct ParticleBin {
    radial_position: Length,
    width: Length,
    size_distribution: SizeDistribution,  // Power-law or evolved
    surface_density: SurfaceDensity,
    velocity_dispersion: Velocity,
    scale_height: Length,
}

struct DiscreteBody {
    orbital_elements: OrbitalElements,
    core_mass: Mass,
    envelope_mass: Mass,
    composition: Composition,
    envelope_state: EnvelopeState,  // None | Hydrostatic | Runaway
}
```

---

## 3. Local Interaction Rules

### 3.1 Particle-Gas Interactions

Aerodynamic coupling is the dominant physics for sub-planetesimal bodies. Two key processes:

**Vertical Settling:**

Particles settle toward the midplane on a timescale controlled by friction:

```
t_settle = (2/π) × (Σ_gas / ρ_m × s) × exp(-z²/2h²)
```

Rather than solving this exactly, define a settling efficiency per timestep:

```rust
fn settling_factor(&self, dt: Duration, gas: &GasDisk) -> f64 {
    let t_settle = self.settling_timescale(gas);
    1.0 - (-dt / t_settle).exp()
}
```

**Radial Drift:**

Peak drift occurs at τ_fric × Ω ≈ 1, with drift velocity:

```
v_r ≈ -2η × v_K × τ_s / (1 + τ_s²)
```

where η ≈ (h/r)² × |d ln P / d ln r| / 2.

This naturally produces the "meter-size barrier" without explicit coding—bodies at the peak drift size lose material to the star unless growth outpaces drift.

### 3.2 Particle-Particle Interactions (Coagulation)

The collision timescale between particles depends on their relative velocities:

```rust
fn collision_timescale(
    &self,
    other: &ParticleBin,
    turbulent_alpha: f64,
) -> Duration {
    let relative_velocity = self.relative_velocity(other, turbulent_alpha);
    let cross_section = self.collision_cross_section(other);
    let number_density = self.number_density();
    
    Duration::from_secs_f64(1.0 / (number_density * cross_section * relative_velocity))
}
```

**Collision Outcomes:**

Rather than deterministically growing, sample outcomes:

```rust
enum CollisionOutcome {
    Sticking { mass_ratio: f64 },      // Full or partial sticking
    Bouncing,                           // No mass transfer
    Fragmentation { fragment_dist: SizeDistribution },
    Erosion { mass_lost: Mass },
}

fn collision_outcome(
    &self,
    impact_velocity: Velocity,
    size_ratio: f64,
    rng: &mut ChaChaRng,
) -> CollisionOutcome {
    let v_frag = self.fragmentation_velocity();  // ~1-10 m/s for silicates
    let v_bounce = self.bouncing_velocity();     // ~0.01-0.1 m/s
    
    match impact_velocity {
        v if v < v_bounce => CollisionOutcome::Sticking { mass_ratio: 1.0 },
        v if v < v_frag => {
            // Probabilistic bouncing/sticking
            let stick_prob = (v_frag - v) / (v_frag - v_bounce);
            if rng.gen::<f64>() < stick_prob {
                CollisionOutcome::Sticking { mass_ratio: 0.5 + 0.5 * rng.gen::<f64>() }
            } else {
                CollisionOutcome::Bouncing
            }
        }
        _ => CollisionOutcome::Fragmentation { ... }
    }
}
```

This produces emergent behavior: in turbulent regions, relative velocities exceed fragmentation thresholds, naturally limiting growth. In quiescent zones (dead zones, outer disk), growth proceeds efficiently.

### 3.3 Gravitational Instability of Particle Layers

The Goldreich-Ward mechanism and streaming instability emerge from local density conditions. Track Toomre Q for the particle layer:

```
Q_particles = σ × Ω / (π × G × Σ_particles)
```

When Q < Q_crit (≈ 1-2), the particle layer fragments:

```rust
fn check_gravitational_instability(&self, gas: &GasDisk) -> Option<PlanetesimalFormation> {
    let q = self.toomre_q();
    let richardson = self.richardson_number(gas);  // Must be > 0.25 to avoid shear instability
    
    match (q < Q_CRIT, richardson > RI_CRIT) {
        (true, true) => {
            // Instability proceeds: form planetesimals
            let fragment_mass = self.jeans_mass();
            Some(PlanetesimalFormation {
                characteristic_mass: fragment_mass,
                formation_timescale: self.freefall_time(),
            })
        }
        (true, false) => None,  // Self-excited turbulence prevents collapse
        (false, _) => None,     // Stable against self-gravity
    }
}
```

This naturally produces planetesimal formation only where conditions allow—outside turbulent zones, at pressure bumps, in settled particle layers—without explicit zone logic.

### 3.4 Protoplanet Growth

Once discrete bodies exist, track their interactions with both particles and other bodies.

**Gravitational Focusing:**

```rust
fn accretion_rate(&self, particle_bin: &ParticleBin) -> MassRate {
    let sigma = particle_bin.velocity_dispersion;
    let v_esc = self.escape_velocity();
    
    let focusing_factor = 1.0 + (v_esc / sigma).powi(2);
    let geometric_rate = particle_bin.surface_density 
        * self.physical_radius().powi(2) 
        * self.orbital_frequency();
    
    geometric_rate * focusing_factor
}
```

This single formula produces three regimes naturally:

- **Low mass bodies** (v_esc ≪ σ): geometric accretion, slow growth
- **Intermediate bodies** (v_esc ~ σ): transition regime
- **Massive bodies** (v_esc ≫ σ): runaway growth, F_g ∝ M^(4/3)

**Dynamical Friction:**

The most massive bodies cool (reduce eccentricity/inclination) while heating the surrounding population:

```rust
fn dynamical_friction_rate(&self, particle_bin: &ParticleBin) -> (f64, f64) {
    // Returns de/dt and di/dt
    let mass_ratio = self.mass / particle_bin.typical_mass();
    let coulomb_log = (mass_ratio).ln();
    
    // Bodies try to reach equipartition: m × σ² = M × σ_M²
    // Massive bodies have lower velocity dispersion
    ...
}
```

### 3.5 Gas Envelope Capture

The critical core mass emerges from envelope physics. Rather than prescribing it, solve for envelope structure:

```rust
fn envelope_state(&self, gas_disk: &GasDisk, accretion_rate: MassRate) -> EnvelopeState {
    let r_hill = self.hill_radius(gas_disk.stellar_mass);
    let r_acc = self.bondi_radius(gas_disk.sound_speed());
    let r_out = r_hill.min(r_acc);
    
    // Solve for envelope mass given boundary conditions
    let envelope = self.hydrostatic_envelope(r_out, gas_disk, accretion_rate);
    
    match envelope {
        Ok(env) if env.mass < self.core_mass => EnvelopeState::Hydrostatic(env),
        Ok(_) => EnvelopeState::Runaway,  // M_env > M_core: runaway begins
        Err(_) => EnvelopeState::None,    // No hydrostatic solution exists
    }
}
```

The critical core mass (typically ~5-20 M⊕) emerges from the physics of radiative/convective transport, opacity, and luminosity—not from a lookup table.

### 3.6 Gas Disk Evolution

The disk itself evolves through:

- **Viscous accretion:** Mass flows inward, angular momentum outward
- **Photoevaporation:** UV/X-ray heating drives mass loss
- **Planet interaction:** Gap opening, Type I/II migration

```rust
fn evolve_gas(&mut self, dt: Duration, bodies: &[DiscreteBody]) {
    // Viscous evolution (self-similar solution or explicit diffusion)
    self.viscous_spreading(dt);
    
    // Photoevaporative mass loss
    self.photoevaporate(dt);
    
    // Gap opening by massive planets
    for body in bodies {
        if let Some(gap_profile) = body.gap_opening_criterion(self) {
            self.apply_gap(body.semi_major_axis, gap_profile);
        }
    }
}
```

Gap opening emerges from the balance between tidal torques and viscous refilling:

```rust
fn gap_opening_criterion(&self, disk: &GasDisk) -> Option<GapProfile> {
    let thermal = self.mass / disk.stellar_mass > disk.aspect_ratio.cubed();
    let viscous = self.mass / disk.stellar_mass > (40.0 * disk.alpha).sqrt() 
                  * disk.aspect_ratio.powi(5);
    
    if thermal && viscous {
        Some(GapProfile::Deep)
    } else if thermal || viscous {
        Some(GapProfile::Partial)
    } else {
        None
    }
}
```

---

## 4. Emergent Phenomena

The beauty of this approach is that complex behaviors emerge without explicit coding:

### 4.1 Runaway Growth → Oligarchic Growth Transition

Initially, all bodies grow slowly. The first body to reach v_esc > σ enters runaway growth and rapidly dominates its feeding zone. But as it grows:

- Viscous stirring heats the planetesimals (σ increases)
- Gravitational focusing decreases (F_g drops)
- Growth rate slows

Simultaneously, other embryos grow in their own feeding zones. The result: an oligarchic regime where multiple similar-mass bodies grow in parallel. This transition emerges from the velocity evolution equations alone.

### 4.2 Snow Line Effects

The condensation physics produces a natural surface density enhancement at the snow line. In the emergent simulation, this produces:

- Faster particle growth (more material)
- Earlier gravitational instability threshold
- Higher isolation masses
- First locations to form giant planet cores

No explicit "is this the snow line?" logic is needed—the physics produces the effect.

### 4.3 Planet Traps

Pressure maxima in the disk naturally trap drifting particles (the radial drift equation shows v_r → 0 where dP/dr = 0). Locations like:

- Snow line (opacity transition)
- Dead zone edges (viscosity transition)
- Gap edges (surface density transition)

These become preferential sites for planetesimal and planet formation, producing architectural features without prescribed formation zones.

### 4.4 Resonant Capture

Bodies undergoing migration naturally approach mean-motion resonances with other planets. The resonant interaction modifies eccentricities and can lead to capture. This produces:

- Resonant chains (like Trappist-1)
- Kozai-Lidov oscillations in inclined systems
- System stability limits

---

## 5. Implementation Strategy

### 5.1 Time Integration

Use adaptive timesteps based on the fastest timescales:

```rust
fn adaptive_timestep(&self) -> Duration {
    let dynamical = self.bodies.iter()
        .map(|b| b.orbital_period())
        .min();
    
    let drift = self.particle_bins.iter()
        .map(|p| p.radial_drift_timescale())
        .min();
    
    let collision = self.particle_bins.iter()
        .map(|p| p.collision_timescale())
        .min();
    
    let disk_evolution = self.gas_disk.viscous_timescale();
    
    [dynamical / 20.0, drift / 10.0, collision / 5.0, disk_evolution / 100.0]
        .iter()
        .copied()
        .min()
}
```

### 5.2 Statistical vs. Direct Treatment

The key to performance is knowing when to use statistical methods:

| Scale | Treatment | Justification |
|-------|-----------|---------------|
| Dust (μm-mm) | Continuous distribution | N ~ 10^15+, central limit applies |
| Pebbles (cm-m) | Binned population | N ~ 10^10, Monte Carlo sampling |
| Planetesimals (1-100 km) | Super-particles | N ~ 10^6, statistical weights |
| Embryos (>1000 km) | Direct N-body | N < 100, individual tracking |
| Planets | Direct N-body | N < 20, high accuracy |

Transitions between representations:

```rust
enum BodyRepresentation {
    Statistical(ParticleBin),
    SuperParticle { weight: f64, representative: DiscreteBody },
    Direct(DiscreteBody),
}

fn promote_to_direct(&mut self, bin: &mut ParticleBin) -> Option<DiscreteBody> {
    // When a single body dominates the bin mass
    let largest = bin.largest_body_mass();
    if largest > bin.total_mass() * 0.1 {
        Some(bin.extract_largest())
    } else {
        None
    }
}
```

### 5.3 Spatial Decomposition

Use dynamic formation zones:

```rust
struct SpatialZone {
    inner_edge: Length,
    outer_edge: Length,
    particle_bins: Vec<ParticleBin>,
    discrete_bodies: Vec<DiscreteBody>,
}

impl SpatialZone {
    fn should_merge_with(&self, other: &SpatialZone) -> bool {
        // Merge if particle populations are well-mixed
        let diffusion_timescale = self.radial_diffusion_time(other);
        let age = self.simulation_time();
        diffusion_timescale < age * 0.1
    }
    
    fn should_split(&self) -> Option<Length> {
        // Split at density discontinuities
        self.find_pressure_maximum()
    }
}
```

---

## 6. Validation Strategy

### 6.1 Unit Tests: Physical Limits

Each interaction rule should reproduce known limits:

```rust
#[test]
fn runaway_growth_scaling() {
    // Verify dM/dt ∝ M^(4/3) in strong focusing limit
    let body = DiscreteBody::new(Mass::earth_masses(0.1), ...);
    let rate1 = body.accretion_rate(low_dispersion_bin);
    
    let bigger = body.with_mass(body.mass * 2.0);
    let rate2 = bigger.accretion_rate(low_dispersion_bin);
    
    let expected_ratio = 2.0_f64.powf(4.0/3.0);
    assert!((rate2 / rate1 - expected_ratio).abs() < 0.1);
}

#[test]
fn critical_core_mass_emergence() {
    // Verify hydrostatic solutions cease to exist above M_crit
    for core_mass in (1..30).map(Mass::earth_masses) {
        let envelope = hydrostatic_envelope(core_mass, disk_conditions);
        if core_mass < Mass::earth_masses(15.0) {
            assert!(envelope.is_ok());
        }
        // The transition should happen somewhere around 10-20 M_earth
    }
}
```

### 6.2 Integration Tests: Reference Systems

Use reference seeds for known-good outcomes:

```rust
#[test]
fn solar_system_analog() {
    let system = simulate_with_seed(SOLAR_SEED, MMSN_DISK);
    
    // Should produce:
    // - Terrestrial planets in inner region
    // - Gas giants near/beyond snow line
    // - Ice giants in outer region
    // - Reasonable mass ratios
    
    assert!(system.has_gas_giant_between(3.0.au(), 10.0.au()));
    assert!(system.terrestrial_planets_count() >= 2);
    assert!(!system.has_hot_jupiter());
}

#[test]
fn hot_jupiter_from_migration() {
    let system = simulate_with_seed(HOT_JUPITER_SEED, HIGH_MASS_DISK);
    
    // High-mass disk should produce massive planet that migrates
    assert!(system.has_planet_with(
        |p| p.mass > Mass::jupiter_masses(0.3) && p.semi_major_axis < 0.1.au()
    ));
}
```

### 6.3 Statistical Validation

Run population synthesis to verify distributions:

```rust
fn validate_population(n_systems: usize) -> PopulationStats {
    let systems: Vec<_> = (0..n_systems)
        .map(|i| simulate_with_seed(i as u64, random_disk_conditions(i)))
        .collect();
    
    PopulationStats {
        planet_occurrence_vs_period: histogram(...),
        mass_distribution: histogram(...),
        eccentricity_distribution: histogram(...),
        metallicity_correlation: regression(...),
    }
}

#[test]
fn matches_kepler_statistics() {
    let stats = validate_population(10_000);
    
    // Super-Earths should be common
    assert!(stats.occurrence_rate(Mass::earth_masses(1.0)..Mass::earth_masses(4.0)) > 0.3);
    
    // Hot Jupiters should be rare
    assert!(stats.occurrence_rate_hot_jupiter() < 0.02);
    
    // Giant planets should correlate with metallicity
    assert!(stats.metallicity_correlation_giants() > 0.3);
}
```

---

## 7. Performance Considerations

### 7.1 Algorithmic Complexity

| Operation | Naive | Optimized |
|-----------|-------|-----------|
| Pairwise particle interactions | O(N²) | O(N log N) with spatial hashing |
| Gravitational N-body | O(N²) | O(N log N) with tree codes |
| Disk-body interactions | O(N × M) | O(N) with precomputed profiles |
| Gas disk evolution | O(M²) | O(M) with implicit methods |

### 7.2 Parallelization Opportunities

```rust
// Particle bins evolve independently except at boundaries
use rayon::prelude::*;

self.particle_bins.par_iter_mut()
    .for_each(|bin| bin.evolve_internal(dt));

// Then handle boundary fluxes (must be sequential)
self.handle_radial_fluxes(dt);
```

### 7.3 Early Termination

For fast generation, implement "boring" detection:

```rust
fn is_system_settled(&self) -> bool {
    // No more gas
    if self.gas_disk.total_mass() < Mass::jupiter_masses(0.01) {
        // And orbits are stable
        self.discrete_bodies.iter().all(|b| b.eccentricity < 0.1)
            && self.is_hill_stable()
    } else {
        false
    }
}
```

---

## 8. Module Structure

```
src/
├── emergent/
│   ├── mod.rs
│   ├── state.rs           # SimulationState, integration loop
│   ├── gas_disk.rs        # GasDisk, viscous evolution, photoevaporation
│   ├── particles/
│   │   ├── mod.rs
│   │   ├── bin.rs         # ParticleBin, size distributions
│   │   ├── aerodynamics.rs # Drag, drift, settling
│   │   └── coagulation.rs  # Collision outcomes, growth
│   ├── bodies/
│   │   ├── mod.rs
│   │   ├── discrete.rs     # DiscreteBody, orbital elements
│   │   ├── accretion.rs    # Gravitational focusing, feeding zones
│   │   └── envelope.rs     # Gas envelope physics, critical mass
│   ├── interactions/
│   │   ├── mod.rs
│   │   ├── gravitational.rs # N-body, dynamical friction, stirring
│   │   ├── migration.rs     # Type I/II migration, resonances
│   │   └── scattering.rs    # Close encounters, ejection
│   └── transitions.rs      # Representation promotions/demotions
├── validation/
│   ├── physics_validator.rs
│   ├── reference_systems.rs
│   └── population_stats.rs
```

---

## 9. Open Questions & Design Decisions

### 9.1 Initial Conditions

**Option A: Start from gas disk only**
- Pro: Most physically complete
- Con: Dust settling phase is fast but computationally tedious

**Option B: Start with particle population at t ~ 10^4 yr**
- Pro: Skips trivial early evolution
- Con: Initial size distribution is a free parameter

**Recommendation:** Option B for speed, with configurable initial distributions that can be validated against Option A runs.

### 9.2 Envelope Resolution

**Option A: Full 1D envelope integration**
- Pro: Accurate critical masses
- Con: Expensive per timestep

**Option B: Tabulated envelope solutions**
- Pro: Fast lookup
- Con: Interpolation errors, complex parameter space

**Option C: Analytic approximations (Stevenson 1982)**
- Pro: Fast, captures scaling
- Con: 30-50% errors in critical mass

**Recommendation:** Option C for exploration, Option B for production.

### 9.3 N-body Threshold

When should bodies transition to direct integration?

**Conservative:** Always track individually above isolation mass (~10^-2 M⊕)

**Aggressive:** Only track gas giants individually, use super-particles for embryos

**Recommendation:** Configurable threshold based on generation speed requirements.

---

## 10. Summary

This emergent simulation architecture:

1. **Defines local physics** (drag, collisions, gravitational interactions) rather than global outcomes
2. **Uses adaptive representations** (statistical → super-particle → direct) based on local conditions  
3. **Lets formation pathways emerge** from physical interactions rather than prescribed rules
4. **Maintains existing patterns** (trait-based architecture, functional methods, seeded RNG)
5. **Supports both quick generation and detailed investigation** via configurable fidelity

The key difference from a hybrid approach: instead of partitioning the disk into zones with different formation rules, the zones and their properties emerge from the physics. Instead of tuning parameters for each planet type, the diversity arises from the sensitivity to initial conditions and stochastic events during evolution.

---

# Part 2: Development Plan

## Overview

This plan follows a bottom-up approach: build and validate foundational physics modules first, then compose them into increasingly complex interactions, and finally integrate into a complete simulation loop. Each phase produces usable, testable artifacts.

**Estimated Timeline:** 12-16 weeks for a fully functional system, with useful intermediate milestones throughout.

---

## Phase 1: Foundations (Weeks 1-3)

### 1.1 Gas Disk Infrastructure

**Goal:** A self-consistent gas disk that evolves viscously and responds to thermal structure.

#### Week 1: Static Disk Profiles

Start with existing infrastructure, ensuring the interfaces support evolution:

```rust
// Core trait for radial profiles
trait RadialProfile<T> {
    fn at(&self, r: Length) -> T;
    fn gradient_at(&self, r: Length) -> T;  // d/dr
    fn integrate(&self, r_in: Length, r_out: Length) -> T;
}

// Disk structure
struct GasDisk {
    surface_density: Box<dyn RadialProfile<SurfaceDensity>>,
    temperature: Box<dyn RadialProfile<Temperature>>,
    // Derived quantities cached and invalidated on mutation
    cache: DiskCache,
}

impl GasDisk {
    fn sound_speed(&self, r: Length) -> Velocity { ... }
    fn scale_height(&self, r: Length) -> Length { ... }
    fn aspect_ratio(&self, r: Length) -> f64 { ... }
    fn midplane_density(&self, r: Length) -> Density { ... }
    fn pressure(&self, r: Length) -> Pressure { ... }
    fn pressure_gradient_parameter(&self, r: Length) -> f64 { ... }  // η
}
```

**Deliverables:**
- `GasDisk` with configurable initial profiles (power-law, self-similar, MMSN)
- Unit tests verifying hydrostatic equilibrium, known profile shapes
- Integration with existing stellar/disk initialization

#### Week 2: Viscous Evolution

Implement the diffusion equation for surface density evolution:

```rust
impl GasDisk {
    fn evolve_viscous(&mut self, dt: Duration) {
        // ∂Σ/∂t = (3/r) ∂/∂r [r^(1/2) ∂/∂r (ν Σ r^(1/2))]
        // Discretize on radial grid, solve implicitly for stability
    }
    
    fn viscosity(&self, r: Length) -> Viscosity {
        // α-prescription: ν = α × c_s × h
        self.alpha * self.sound_speed(r) * self.scale_height(r)
    }
}
```

**Deliverables:**
- Viscous spreading that matches self-similar solutions
- Configurable α parameter (spatially varying for dead zones later)
- Conservation tests: total angular momentum, mass accounting

#### Week 3: Disk Dispersal

Add photoevaporation and verify disk lifetime behavior:

```rust
enum PhotoevaporationModel {
    None,
    EUV { phi: PhotonFlux },           // ~10^41 photons/s
    Xray { luminosity: Luminosity },   // ~10^30 erg/s
    Combined { euv: PhotonFlux, xray: Luminosity },
}

impl GasDisk {
    fn photoevaporation_rate(&self, r: Length) -> MassRate {
        match &self.photoevaporation {
            PhotoevaporationModel::EUV { phi } => {
                // Mass loss beyond critical radius r_g
                let r_g = self.gravitational_radius();
                if r > r_g { /* EUV wind formula */ } else { MassRate::zero() }
            }
            // ...
        }
    }
    
    fn evolve(&mut self, dt: Duration) {
        self.evolve_viscous(dt);
        self.apply_photoevaporation(dt);
        self.invalidate_cache();
    }
}
```

**Deliverables:**
- Disk lifetimes of 1-10 Myr depending on parameters
- Two-timescale dispersal (slow viscous decline → rapid inside-out clearing)
- Validation against observed disk fraction vs. cluster age

**Phase 1 Milestone:** Stand-alone gas disk evolution that can be visualized and validated independently.

---

## Phase 2: Particle Dynamics (Weeks 4-6)

### 2.1 Single-Particle Aerodynamics

**Goal:** Correct drag physics across all size regimes.

#### Week 4: Drag Laws and Stopping Times

```rust
#[derive(Clone, Copy)]
enum DragRegime {
    Epstein,      // s < 9λ/4 (mean free path)
    Stokes,       // Re < 1
    Intermediate, // 1 < Re < 800
    Quadratic,    // Re > 800
}

struct Particle {
    size: Length,
    material_density: Density,
}

impl Particle {
    fn drag_regime(&self, gas: &GasDisk, r: Length) -> DragRegime {
        let mfp = gas.mean_free_path(r);
        let reynolds = self.reynolds_number(gas, r);
        
        match (self.size < 9.0 * mfp / 4.0, reynolds) {
            (true, _) => DragRegime::Epstein,
            (false, re) if re < 1.0 => DragRegime::Stokes,
            (false, re) if re < 800.0 => DragRegime::Intermediate,
            _ => DragRegime::Quadratic,
        }
    }
    
    fn stopping_time(&self, gas: &GasDisk, r: Length) -> Duration {
        let rho_g = gas.midplane_density(r);
        let v_th = gas.thermal_velocity(r);
        
        match self.drag_regime(gas, r) {
            DragRegime::Epstein => {
                // t_s = ρ_m × s / (ρ_g × v_th)
                Duration::from_secs_f64(
                    (self.material_density * self.size) / (rho_g * v_th)
                )
            }
            DragRegime::Stokes => {
                // t_s = 2 ρ_m s² / (9 μ)  where μ is dynamic viscosity
                ...
            }
            // ...
        }
    }
    
    fn stokes_number(&self, gas: &GasDisk, r: Length) -> f64 {
        // τ_s = t_s × Ω
        self.stopping_time(gas, r) * gas.orbital_frequency(r)
    }
}
```

**Deliverables:**
- Correct regime transitions
- Stokes number profiles showing τ_s = 1 at characteristic sizes
- Validation: known analytic limits

#### Week 5: Radial Drift

```rust
impl Particle {
    fn radial_drift_velocity(&self, gas: &GasDisk, r: Length) -> Velocity {
        let tau = self.stokes_number(gas, r);
        let eta = gas.pressure_gradient_parameter(r);
        let v_k = gas.keplerian_velocity(r);
        
        // v_r = -2 η v_K τ_s / (1 + τ_s²)
        -2.0 * eta * v_k * tau / (1.0 + tau.powi(2))
    }
    
    fn azimuthal_drift_velocity(&self, gas: &GasDisk, r: Length) -> Velocity {
        let tau = self.stokes_number(gas, r);
        let eta = gas.pressure_gradient_parameter(r);
        let v_k = gas.keplerian_velocity(r);
        
        // Δv_φ = -η v_K / (1 + τ_s²)
        -eta * v_k / (1.0 + tau.powi(2))
    }
}
```

**Deliverables:**
- Peak drift velocity at τ_s = 1
- Drift timescales matching textbook values (~100 yr for meter-sized at 1 AU)
- Visualization: drift velocity vs. size and radius

#### Week 6: Vertical Settling

```rust
impl Particle {
    fn settling_velocity(&self, gas: &GasDisk, r: Length, z: Length) -> Velocity {
        let tau = self.stokes_number(gas, r);
        let omega = gas.orbital_frequency(r);
        
        // v_z = -τ_s × Ω × z  (for τ_s << 1)
        // More complex for larger τ_s
        if tau < 0.1 {
            -tau * omega * z
        } else {
            // Full solution accounting for gas density stratification
            ...
        }
    }
    
    fn equilibrium_scale_height(&self, gas: &GasDisk, r: Length, alpha: f64) -> Length {
        let tau = self.stokes_number(gas, r);
        let h_g = gas.scale_height(r);
        
        // h_d / h_g = sqrt(α / (α + τ_s))
        h_g * (alpha / (alpha + tau)).sqrt()
    }
}
```

**Deliverables:**
- Settling timescales
- Equilibrium particle scale heights vs. turbulent α
- Dust-to-gas ratio enhancement at midplane

**Phase 2 Milestone:** Single-particle trajectories in evolving gas disk, visualizable as particle tracks showing drift and settling.

---

## Phase 3: Particle Populations (Weeks 7-9)

### 3.1 Statistical Representation

**Goal:** Efficient treatment of large particle populations with size distributions.

#### Week 7: Particle Bins

```rust
struct SizeDistribution {
    variant: SizeDistributionVariant,
}

enum SizeDistributionVariant {
    /// Power law: n(s) ∝ s^(-q), typically q ≈ 3.5 (MRN)
    PowerLaw { 
        s_min: Length, 
        s_max: Length, 
        exponent: f64,
        normalization: f64,
    },
    /// Discretized into logarithmic bins
    Binned {
        bin_edges: Vec<Length>,
        number_density: Vec<f64>,
    },
    /// Single size (useful for testing)
    Monodisperse { size: Length },
}

impl SizeDistribution {
    fn total_mass(&self, material_density: Density) -> Mass { ... }
    fn mean_size(&self) -> Length { ... }
    fn size_weighted_mean(&self, weight: impl Fn(Length) -> f64) -> f64 { ... }
    fn evolve(&mut self, collision_kernel: &CollisionKernel, dt: Duration) { ... }
}

struct ParticleBin {
    radial_center: Length,
    radial_width: Length,
    size_distribution: SizeDistribution,
    surface_density: SurfaceDensity,
    velocity_dispersion: Velocity,
    scale_height: Length,
}
```

**Deliverables:**
- ParticleBin construction from disk initial conditions
- Mass-conserving operations (split, merge, transfer)
- Size distribution evolution framework

#### Week 8: Coagulation Physics

```rust
struct CollisionKernel {
    material: MaterialProperties,
}

impl CollisionKernel {
    fn collision_rate(
        &self,
        size1: Length,
        size2: Length,
        relative_velocity: Velocity,
        number_density: f64,
    ) -> f64 {
        let cross_section = PI * (size1 + size2).powi(2);
        number_density * cross_section * relative_velocity
    }
    
    fn relative_velocity(
        &self,
        size1: Length,
        size2: Length,
        gas: &GasDisk,
        r: Length,
        alpha: f64,
    ) -> Velocity {
        // Sum of contributions:
        let v_brownian = self.brownian_velocity(size1, size2, gas.temperature(r));
        let v_settling = self.differential_settling_velocity(size1, size2, gas, r);
        let v_drift = self.differential_drift_velocity(size1, size2, gas, r);
        let v_turbulent = self.turbulent_velocity(size1, size2, gas, r, alpha);
        
        // Add in quadrature
        (v_brownian.powi(2) + v_settling.powi(2) + v_drift.powi(2) + v_turbulent.powi(2)).sqrt()
    }
    
    fn collision_outcome(
        &self,
        size1: Length,
        size2: Length,
        impact_velocity: Velocity,
        rng: &mut ChaChaRng,
    ) -> CollisionOutcome {
        let v_frag = self.fragmentation_velocity(size1, size2);
        let v_bounce = self.bouncing_velocity(size1, size2);
        // Probabilistic outcomes based on velocity thresholds
        ...
    }
}

enum CollisionOutcome {
    PerfectSticking,
    PartialSticking { efficiency: f64 },
    Bouncing,
    MassTransfer { from_larger: bool, fraction: f64 },
    Fragmentation { power_law_exponent: f64 },
    Erosion { mass_loss_fraction: f64 },
}
```

**Deliverables:**
- Collision rates matching textbook estimates
- Outcome probabilities as function of size and velocity
- Unit tests for limiting cases (all-sticking, all-fragmentation)

#### Week 9: Size Distribution Evolution

Implement the Smoluchowski coagulation equation:

```rust
impl SizeDistribution {
    fn evolve_coagulation(
        &mut self,
        kernel: &CollisionKernel,
        gas: &GasDisk,
        r: Length,
        alpha: f64,
        dt: Duration,
    ) {
        match &mut self.variant {
            SizeDistributionVariant::Binned { bin_edges, number_density } => {
                // dn_k/dt = (1/2) Σ_{i+j=k} A_ij n_i n_j - n_k Σ_i A_ik n_i
                let mut dn = vec![0.0; number_density.len()];
                
                for i in 0..number_density.len() {
                    for j in 0..number_density.len() {
                        let v_rel = kernel.relative_velocity(
                            self.bin_center(i), self.bin_center(j), gas, r, alpha
                        );
                        let rate = kernel.collision_rate(
                            self.bin_center(i), self.bin_center(j), v_rel, 1.0
                        );
                        // Gain and loss terms...
                    }
                }
                // Apply changes with mass conservation
            }
            _ => { /* Convert to binned first */ }
        }
    }
}
```

**Deliverables:**
- Coagulation evolution matching known solutions (constant kernel, product kernel)
- Runaway growth emergence for product kernel
- Fragmentation equilibrium when velocities exceed threshold

**Phase 3 Milestone:** Particle populations that evolve their size distributions in an evolving gas disk, showing growth from dust to pebbles and the emergence of growth barriers.

---

## Phase 4: Gravitational Dynamics (Weeks 10-12)

### 4.1 Particle Layer Instability

**Goal:** Planetesimal formation via gravitational instability.

#### Week 10: Stability Criteria

```rust
impl ParticleBin {
    fn toomre_q(&self, stellar_mass: Mass, r: Length) -> f64 {
        let omega = orbital_frequency(stellar_mass, r);
        // Q = σ Ω / (π G Σ)
        self.velocity_dispersion * omega / (PI * G * self.surface_density)
    }
    
    fn richardson_number(&self, gas: &GasDisk, r: Length) -> f64 {
        // Ri = N² / (dv/dz)²
        // Simplified calculation...
        let particle_density = self.midplane_density();
        let gas_density = gas.midplane_density(r);
        let rho_ratio = particle_density / gas_density;
        ...
    }
    
    fn is_gravitationally_unstable(&self, gas: &GasDisk, r: Length, stellar_mass: Mass) -> bool {
        let q = self.toomre_q(stellar_mass, r);
        let ri = self.richardson_number(gas, r);
        
        // Must satisfy both: Q < Q_crit AND Ri > Ri_crit
        q < Q_CRIT && ri > RI_CRIT
    }
}
```

**Deliverables:**
- Stability maps (Q, Ri) as function of disk parameters
- Identification of conditions allowing gravitational collapse
- Connection to streaming instability literature

#### Week 11: Planetesimal Formation

```rust
struct PlanetesimalFormationEvent {
    location: Length,
    total_mass: Mass,
    characteristic_mass: Mass,
    size_distribution: SizeDistribution,
    formation_timescale: Duration,
}

impl ParticleBin {
    fn attempt_planetesimal_formation(
        &mut self,
        gas: &GasDisk,
        stellar_mass: Mass,
        rng: &mut ChaChaRng,
    ) -> Option<PlanetesimalFormationEvent> {
        if !self.is_gravitationally_unstable(gas, self.radial_center, stellar_mass) {
            return None;
        }
        
        // Jeans mass in the particle layer
        let sigma = self.velocity_dispersion;
        let surface_density = self.surface_density;
        let lambda_j = sigma.powi(2) / (G * surface_density);
        let m_j = PI * lambda_j.powi(2) * surface_density;
        
        // Formation efficiency
        let efficiency = rng.gen_range(0.1..0.5);
        let mass_to_convert = self.total_mass() * efficiency;
        
        // Remove mass from particle bin
        self.surface_density -= mass_to_convert / self.area();
        
        Some(PlanetesimalFormationEvent {
            location: self.radial_center,
            total_mass: mass_to_convert,
            characteristic_mass: m_j,
            size_distribution: SizeDistribution::power_law(-2.5, m_j * 0.01, m_j * 10.0),
            formation_timescale: self.freefall_time(),
        })
    }
}
```

**Deliverables:**
- Planetesimal formation at appropriate locations
- Characteristic sizes of 1-100 km
- Mass transfer from particle bins to discrete body population

#### Week 12: Discrete Body Dynamics

```rust
struct DiscreteBody {
    id: BodyId,
    mass: Mass,
    physical_radius: Length,
    orbital_elements: OrbitalElements,
    composition: Composition,
}

struct OrbitalElements {
    semi_major_axis: Length,
    eccentricity: f64,
    inclination: Angle,
    longitude_of_ascending_node: Angle,
    argument_of_periapsis: Angle,
    mean_anomaly: Angle,
}

impl DiscreteBody {
    fn hill_radius(&self, stellar_mass: Mass) -> Length {
        self.orbital_elements.semi_major_axis 
            * (self.mass / (3.0 * stellar_mass)).powf(1.0/3.0)
    }
    
    fn escape_velocity(&self) -> Velocity {
        (2.0 * G * self.mass / self.physical_radius).sqrt()
    }
    
    fn gravitational_focusing_factor(&self, velocity_dispersion: Velocity) -> f64 {
        let v_esc = self.escape_velocity();
        1.0 + (v_esc / velocity_dispersion).powi(2)
    }
}

struct NBodyIntegrator {
    bodies: Vec<DiscreteBody>,
    stellar_mass: Mass,
}

impl NBodyIntegrator {
    fn step(&mut self, dt: Duration) {
        // Leapfrog or higher-order symplectic integrator
        self.kick(dt / 2.0);
        self.drift(dt);
        self.kick(dt / 2.0);
    }
}
```

**Deliverables:**
- Stable orbital integration over 10^6+ orbits
- Energy conservation tests
- Close encounter detection

**Phase 4 Milestone:** Complete particle-to-planetesimal pathway, with discrete bodies forming from unstable particle layers and evolving under mutual gravity.

---

## Phase 5: Protoplanet Growth (Weeks 13-14)

### 5.1 Accretion from Particle Disk

#### Week 13: Growth Physics

```rust
impl DiscreteBody {
    fn accretion_rate_from_particles(
        &self,
        particle_bin: &ParticleBin,
        stellar_mass: Mass,
    ) -> MassRate {
        let sigma = particle_bin.velocity_dispersion;
        let sigma_p = particle_bin.surface_density;
        let omega = orbital_frequency(stellar_mass, self.orbital_elements.semi_major_axis);
        
        let f_g = self.gravitational_focusing_factor(sigma);
        
        // dM/dt = √3/2 × Σ_p × Ω × R² × F_g
        (3.0_f64.sqrt() / 2.0) * sigma_p * omega * self.physical_radius.powi(2) * f_g
    }
    
    fn feeding_zone_width(&self, stellar_mass: Mass) -> Length {
        2.5 * self.hill_radius(stellar_mass)
    }
    
    fn isolation_mass(&self, surface_density: SurfaceDensity, stellar_mass: Mass) -> Mass {
        let a = self.orbital_elements.semi_major_axis;
        let delta_a = self.feeding_zone_width(stellar_mass);
        2.0 * PI * a * delta_a * surface_density
    }
}

impl ParticleBin {
    fn viscous_stirring_rate(&self, protoplanet: &DiscreteBody, stellar_mass: Mass) -> f64 {
        // dσ²/dt from gravitational scattering
        ...
    }
    
    fn dynamical_friction_on(&self, protoplanet: &DiscreteBody) -> (f64, f64) {
        // de/dt, di/dt for the protoplanet
        ...
    }
}
```

**Deliverables:**
- Growth rates matching Armitage equations
- Runaway growth emergence (dM/dt ∝ M^{4/3})
- Transition to oligarchic growth when stirring dominates

#### Week 14: Regime Transitions

```rust
enum GrowthRegime {
    Geometric,
    Runaway,
    Oligarchic,
    Isolation,
}

impl DiscreteBody {
    fn growth_regime(
        &self,
        particle_bin: &ParticleBin,
        stellar_mass: Mass,
    ) -> GrowthRegime {
        let f_g = self.gravitational_focusing_factor(particle_bin.velocity_dispersion);
        let m_iso = self.isolation_mass(particle_bin.surface_density, stellar_mass);
        let stirring_dominated = self.viscous_stirring_dominates(particle_bin, stellar_mass);
        
        match (f_g, self.mass / m_iso, stirring_dominated) {
            (fg, _, _) if fg < 2.0 => GrowthRegime::Geometric,
            (_, ratio, _) if ratio > 0.8 => GrowthRegime::Isolation,
            (_, _, true) => GrowthRegime::Oligarchic,
            _ => GrowthRegime::Runaway,
        }
    }
}
```

**Deliverables:**
- Automatic regime detection
- Smooth transitions between growth modes
- Validation against N-body results from literature

**Phase 5 Milestone:** Protoplanet growth from planetesimals showing emergent runaway and oligarchic behavior.

---

## Phase 6: Gas Giant Formation (Weeks 15-16)

### 6.1 Envelope Physics

#### Week 15: Hydrostatic Envelopes

```rust
struct GasEnvelope {
    core_mass: Mass,
    envelope_mass: Mass,
    outer_radius: Length,
    structure: EnvelopeStructure,
}

enum EnvelopeStructure {
    Isothermal { temperature: Temperature },
    Radiative { opacity: Opacity, luminosity: Luminosity },
    Tabulated { profile: EnvelopeProfile },
}

impl GasEnvelope {
    fn hydrostatic_solution(
        core_mass: Mass,
        disk: &GasDisk,
        r: Length,
        accretion_rate: MassRate,
        opacity: Opacity,
    ) -> Result<Self, NoSolutionError> {
        let r_hill = hill_radius(core_mass, disk.stellar_mass, r);
        let r_bondi = bondi_radius(core_mass, disk.sound_speed(r));
        let r_out = r_hill.min(r_bondi);
        
        let luminosity = G * core_mass * accretion_rate / core_radius(core_mass);
        // Integrate structure equations...
        ...
    }
}

impl DiscreteBody {
    fn critical_core_mass(
        &self,
        disk: &GasDisk,
        accretion_rate: MassRate,
        opacity: Opacity,
    ) -> Mass {
        // M_crit ~ 7 × (Ṁ / 10^-7 M_⊕/yr)^0.25 × (κ / 1 cm²/g)^0.25 M_⊕
        let m_dot_normalized = accretion_rate / MassRate::earth_masses_per_year(1e-7);
        let kappa_normalized = opacity / Opacity::new(1.0);
        
        Mass::earth_masses(7.0) 
            * m_dot_normalized.powf(0.25) 
            * kappa_normalized.powf(0.25)
    }
}
```

**Deliverables:**
- Envelope mass as function of core mass
- Critical core mass emergence (~10-20 M⊕)
- Dependence on opacity and accretion rate

#### Week 16: Runaway Gas Accretion

```rust
enum EnvelopeState {
    None,
    Hydrostatic(GasEnvelope),
    Runaway { accretion_rate: MassRate, gap_opening: bool },
    Final(GasEnvelope),
}

impl DiscreteBody {
    fn evolve_envelope(
        &mut self,
        disk: &GasDisk,
        r: Length,
        planetesimal_accretion_rate: MassRate,
        dt: Duration,
    ) {
        match &mut self.envelope_state {
            EnvelopeState::None => {
                if self.mass > self.minimum_envelope_mass(disk, r) {
                    self.envelope_state = EnvelopeState::Hydrostatic(
                        GasEnvelope::hydrostatic_solution(...).unwrap()
                    );
                }
            }
            EnvelopeState::Hydrostatic(env) => {
                let m_crit = self.critical_core_mass(disk, planetesimal_accretion_rate, ...);
                if self.core_mass > m_crit {
                    self.envelope_state = EnvelopeState::Runaway {
                        accretion_rate: self.supply_limited_accretion_rate(disk, r),
                        gap_opening: self.opens_gap(disk, r),
                    };
                } else {
                    let growth_rate = env.kelvin_helmholtz_rate(self.core_mass);
                    env.envelope_mass += growth_rate * dt;
                }
            }
            EnvelopeState::Runaway { accretion_rate, gap_opening } => {
                let dm = *accretion_rate * dt;
                self.envelope_mass += dm;
                *gap_opening = self.opens_gap(disk, r);
                
                if disk.local_mass(r) < self.envelope_mass * 0.1 || *gap_opening {
                    *accretion_rate = self.gap_limited_rate(disk, r);
                }
            }
            EnvelopeState::Final(_) => {}
        }
    }
}
```

**Deliverables:**
- Complete giant planet formation pathway
- Formation timescales of 1-10 Myr
- Final masses depending on disk properties and timing

**Phase 6 Milestone:** Giant planet formation via core accretion, with critical core mass, runaway accretion, and gap opening all emerging from the physics.

---

## Phase 7: Integration & Migration (Weeks 17-18)

### 7.1 Full Simulation Loop

#### Week 17: System Integration

```rust
struct EmergentSimulation {
    state: SimulationState,
    config: SimulationConfig,
    history: SimulationHistory,
}

impl EmergentSimulation {
    fn step(&mut self) {
        let dt = self.adaptive_timestep();
        
        // 1. Gas disk evolution
        self.state.gas_disk.evolve(dt);
        
        // 2. Particle aerodynamics (drift, settling)
        for bin in &mut self.state.particle_bins {
            bin.evolve_aerodynamics(&self.state.gas_disk, dt);
        }
        
        // 3. Particle coagulation/fragmentation
        for bin in &mut self.state.particle_bins {
            bin.evolve_coagulation(&self.state.gas_disk, dt, &mut self.state.rng);
        }
        
        // 4. Check for planetesimal formation
        let mut new_bodies = Vec::new();
        for bin in &mut self.state.particle_bins {
            if let Some(event) = bin.attempt_planetesimal_formation(
                &self.state.gas_disk, 
                self.state.star.mass(),
                &mut self.state.rng,
            ) {
                new_bodies.extend(self.spawn_planetesimals(event));
            }
        }
        self.state.discrete_bodies.extend(new_bodies);
        
        // 5. N-body integration of discrete bodies
        self.integrate_orbits(dt);
        
        // 6. Body-particle interactions (accretion, stirring)
        for body in &mut self.state.discrete_bodies {
            let local_bin = self.find_local_bin(body);
            body.accrete_from(local_bin, dt);
            local_bin.apply_stirring_from(body, dt);
        }
        
        // 7. Envelope evolution
        for body in &mut self.state.discrete_bodies {
            body.evolve_envelope(&self.state.gas_disk, body.semi_major_axis(), dt);
        }
        
        // 8. Migration
        for body in &mut self.state.discrete_bodies {
            body.apply_migration(&self.state.gas_disk, dt);
        }
        
        // 9. Body-body interactions (collisions, scattering)
        self.handle_close_encounters();
        
        // 10. Bookkeeping
        self.state.time += dt;
        self.record_history();
    }
    
    fn run_until(&mut self, condition: impl Fn(&SimulationState) -> bool) {
        while !condition(&self.state) {
            self.step();
        }
    }
}
```

#### Week 18: Migration Physics

```rust
impl DiscreteBody {
    fn migration_torque(&self, disk: &GasDisk) -> Torque {
        let r = self.orbital_elements.semi_major_axis;
        
        if self.opens_gap(disk, r) {
            self.type_ii_torque(disk)
        } else {
            self.type_i_torque(disk)
        }
    }
    
    fn type_i_torque(&self, disk: &GasDisk) -> Torque {
        // Tanaka et al. (2002)
        ...
    }
    
    fn type_ii_torque(&self, disk: &GasDisk) -> Torque {
        // Locked to viscous evolution of gap edge
        ...
    }
    
    fn opens_gap(&self, disk: &GasDisk, r: Length) -> bool {
        let h_r = disk.aspect_ratio(r);
        let q = self.mass / disk.stellar_mass;
        let alpha = disk.alpha;
        
        let thermal = q > h_r.powi(3);
        let viscous = q > 40.0 * alpha * h_r.powi(5);
        
        thermal && viscous
    }
}
```

**Phase 7 Milestone:** Complete simulation loop with all physics integrated.

---

## Phase 8: Validation & Tuning (Weeks 19-20)

### 8.1 Reference System Validation

```rust
#[cfg(test)]
mod validation {
    const SOLAR_ANALOG_SEED: u64 = 0x534F4C4152;
    
    #[test]
    fn solar_system_analog() {
        let config = SimulationConfig {
            disk: DiskConfig::mmsn(),
            star: StellarConfig::solar(),
            ..Default::default()
        };
        
        let mut sim = EmergentSimulation::new(config, SOLAR_ANALOG_SEED);
        sim.run_until(|s| s.gas_disk.is_dispersed() && s.is_dynamically_settled());
        
        let system = sim.final_system();
        
        let terrestrials = system.planets_in_range(0.3.au(), 2.0.au());
        assert!(terrestrials.len() >= 2);
        assert!(terrestrials.iter().all(|p| p.mass < Mass::earth_masses(3.0)));
        
        let giants = system.planets_with(|p| p.mass > Mass::jupiter_masses(0.1));
        assert!(!giants.is_empty());
        assert!(giants.iter().any(|p| p.semi_major_axis > 3.0.au()));
        
        let hot_jupiters = system.planets_with(|p| 
            p.mass > Mass::jupiter_masses(0.3) && p.semi_major_axis < 0.1.au()
        );
        assert!(hot_jupiters.is_empty());
    }
}
```

### 8.2 Population Statistics

```rust
fn validate_population_statistics(n_systems: usize) {
    let systems: Vec<_> = (0..n_systems)
        .into_par_iter()
        .map(|i| {
            let disk = DiskConfig::from_distribution(i as u64);
            let star = StellarConfig::from_distribution(i as u64);
            let mut sim = EmergentSimulation::new(
                SimulationConfig { disk, star, ..Default::default() },
                i as u64,
            );
            sim.run_to_completion();
            sim.final_system()
        })
        .collect();
    
    let stats = PopulationStatistics::from_systems(&systems);
    
    println!("Super-Earth occurrence (P < 100d): {:.1}%", 
             stats.occurrence_rate(1.0..4.0, 0.0..0.3) * 100.0);
    println!("Hot Jupiter occurrence: {:.2}%",
             stats.hot_jupiter_rate() * 100.0);
    println!("Giant planet metallicity correlation: r = {:.2}",
             stats.metallicity_correlation_giants());
}
```

---

## Deliverables Summary

| Phase | Weeks | Key Deliverable |
|-------|-------|-----------------|
| 1 | 1-3 | Evolving gas disk with viscosity and photoevaporation |
| 2 | 4-6 | Single-particle aerodynamics (drag, drift, settling) |
| 3 | 7-9 | Particle populations with coagulation/fragmentation |
| 4 | 10-12 | Gravitational instability and discrete body dynamics |
| 5 | 13-14 | Protoplanet growth with runaway/oligarchic regimes |
| 6 | 15-16 | Gas envelope physics and giant planet formation |
| 7 | 17-18 | Full simulation loop with migration |
| 8 | 19-20 | Validation suite and population statistics |

---

## Risk Mitigation

**Performance Risk:** Simulation too slow for rapid generation
- Mitigation: Aggressive use of statistical representations, early termination, parallelization
- Fallback: Configurable fidelity levels (fast/balanced/accurate)

**Complexity Risk:** Too many interacting systems to debug
- Mitigation: Each phase produces independently testable module
- Fallback: Can freeze certain physics for debugging others

**Physics Risk:** Emergent outcomes don't match observations
- Mitigation: Extensive unit tests against known analytic solutions
- Fallback: Identify which local rules need adjustment rather than adding global corrections

**Scope Risk:** 20 weeks is optimistic
- Mitigation: Phases 1-5 produce useful intermediate products
- Fallback: Defer migration (Phase 7) as optional enhancement

---

## Success Criteria

The simulation is successful if:

1. **Diversity:** Produces qualitatively different system architectures (Solar-like, compact multi-planet, hot Jupiter, ice giant dominated) from different initial conditions
2. **Plausibility:** All systems pass PhysicsValidator checks
3. **Emergence:** Formation pathways arise from local rules without explicit regime-switching logic
4. **Performance:** Can generate 100+ systems per minute at "fast" fidelity
5. **Reproducibility:** Identical seeds produce identical systems

---

# Part 3: Getting Started

## First Step

Start with **Phase 1, Week 1: Static Disk Profiles**, but scoped tightly. Before building anything new, audit existing code.

The actual first step:

```rust
// Create a single file: src/emergent/gas_disk.rs
// Goal: Define the interface, wrap existing disk, add what's missing
```

Specifically, write these three things:

1. **The `GasDisk` trait** that captures what the emergent simulation needs from a disk (surface density, temperature, sound speed, pressure gradient parameter η). Don't implement yet—just the trait.

2. **An adapter** that wraps existing disk structures to implement this trait. This validates current code is sufficient or identifies gaps.

3. **One unit test** that constructs a disk and verifies `pressure_gradient_parameter()` returns sensible values (~0.002-0.005 for typical disks).

That's it. One file, one trait, one adapter, one test. Maybe 100-150 lines.

## Keeping It Under Control

Three rules:

### Rule 1: Each PR adds exactly one capability with one test

Not "implement particle aerodynamics." Instead:

- PR 1: `stopping_time()` for Epstein drag + test against analytic formula
- PR 2: `stopping_time()` for Stokes drag + test at regime boundary
- PR 3: `radial_drift_velocity()` + test that peak occurs at τ_s = 1

If a PR touches more than two files (excluding tests), it's too big.

### Rule 2: No forward references

Don't build infrastructure for things not yet reached. When in Phase 2 (particle dynamics), don't add hooks for "future envelope physics." Build only what the current phase requires, refactor later.

### Rule 3: Milestone demos, not milestone documents

At each phase milestone, run something and see output:

- Phase 1: Plot disk surface density evolving over 10 Myr
- Phase 2: Animate a particle's radial drift trajectory
- Phase 3: Plot size distribution evolution showing growth stalling at fragmentation barrier
- Phase 4: Show planetesimals spawning at pressure maximum

If it can't be visualized, it's not understood well enough yet.
