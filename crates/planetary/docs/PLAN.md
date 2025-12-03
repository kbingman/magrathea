# Statistical Planetary System Generator

## Project Overview

This document describes a fast, physics-informed planetary system generator designed to complement detailed formation simulations. The approach mirrors the stellar generation system: just as the Kroupa IMF provides a statistically-valid way to sample stellar masses, Kepler/TESS occurrence rates provide a foundation for sampling planetary populations.

### Design Philosophy

**Speed over simulation**: This is explicitly *not* a physics simulation. It's a statistical sampler that produces believable planetary systems in ~1ms, suitable for:

- Procedural generation for games/visualization
- Monte Carlo population studies
- Initial conditions for N-body integrations
- Quick "what if" explorations

**Physics-informed, not physics-derived**: The distributions and correlations are grounded in observational data and theoretical understanding, but planets aren't "formed" - they're sampled from occurrence rate distributions and filtered for plausibility.

**Complementary to emergent simulation**: The existing emergent simulation project (outlined in `emergent_simulation_plan.md`) takes the opposite approach: define local physics rules and let planetary architectures emerge. Both approaches have value:

| Aspect | Statistical Generator | Emergent Simulation |
|--------|----------------------|---------------------|
| Speed | ~1ms/system | ~minutes to hours |
| Diversity | High (stochastic sampling) | High (sensitive to ICs) |
| Physics fidelity | Low (correlations only) | High (actual dynamics) |
| Use case | Population generation | Understanding formation |

---

## Theoretical Foundation

### Two-Tier Classification System

The generator uses a physics-based taxonomy developed from Wolfgang et al. (2016) and related work. This avoids problematic nomenclature like "Super-Earth" (which conflates mass and composition) in favor of a two-tier system:

#### Tier 1: PlanetClass (Physical Regime)

Classification by what pressure support mechanism dominates - purely a function of mass:

| Class | Mass Range | Physical Regime | Key Behavior |
|-------|------------|-----------------|--------------|
| **Rocky** | M < 2 M⊕ | Self-compression | ρ increases with M, R ∝ M^0.29 |
| **Transitional** | 2–5 M⊕ | Envelope retention threshold | Stellar flux determines rocky vs. gaseous |
| **Volatile** | 5–160 M⊕ | Envelope accumulation | ρ *decreases* with M, R ∝ M^0.56 |
| **Giant** | > 160 M⊕ | Electron degeneracy | R nearly constant, ρ ∝ M |

The **Transitional** class is particularly interesting - it straddles the Fulton gap (radius valley). Whether a 3 M⊕ planet becomes a rocky Super-Terran or a puffy Mini-Neptune depends on stellar irradiation and photoevaporation history.

#### Tier 2: PlanetType (Observable Expression)

How a planet of a given class manifests based on temperature, composition, and stellar environment:

**Rocky expressions**: SubEarth, Barren, Lava, Desert, Frozen, Oceanic, Terran, Eyeball, Carbon, Iron

**Transitional expressions**: SuperTerran, MiniNeptune, WaterWorld, Hycean, Chthonian

**Volatile expressions**: SubNeptune, WarmNeptune, IceGiant

**Giant expressions**: GasGiant, HotJupiter, PuffySaturn

### Key Physical Correlations

The generator encodes several well-established correlations:

#### Giant-Metallicity Correlation

```
P(giant) ∝ 10^(2 × [Fe/H])
```

Stars with [Fe/H] = +0.3 are ~4× more likely to host giants than [Fe/H] = -0.3 stars. This reflects the core accretion requirement for solid material.

#### M-Dwarf Compact Systems

M dwarfs preferentially host compact multi-planet systems (TRAPPIST-1 style):
- Lower disk masses → smaller planets
- Closer-in snow lines → more ice-rich compositions
- Longer disk lifetimes → more time for migration and packing

#### Radius Valley / Photoevaporation

Transitional planets (2-5 M⊕) near their stars get stripped of H/He envelopes:

```
F > F_threshold → SuperTerran (rocky)
F < F_threshold → MiniNeptune (retained envelope)
```

The threshold scales with mass: more massive cores retain envelopes against stronger irradiation.

#### Stability Constraints

Generated systems must satisfy the Hill stability criterion:

```
Δa > 8-10 × R_Hill,mutual
```

This prevents dynamically unstable configurations from being output.

---

## Architecture

### Crate Structure

```
planetary/
├── Cargo.toml
└── src/
    ├── lib.rs              # Crate root, re-exports
    ├── planet_class.rs     # PlanetClass enum (4 regimes)
    ├── planet_type.rs      # PlanetType enum (21 expressions)
    ├── planet.rs           # Planet struct, M-R relations, factories
    ├── composition.rs      # Composition (iron/silicate/water/H-He)
    ├── sampling.rs         # Occurrence rates, period/mass sampling
    ├── system.rs           # PlanetarySystem, SystemArchitecture
    └── generation.rs       # Main pipeline: generate_planetary_system()
```

### Dependencies

```toml
[dependencies]
units = { path = "../units" }      # Your units crate
serde = { version = "1.0", features = ["derive"] }
rand = "0.8"
rand_chacha = "0.3"
```

### Core Types

```rust
// Primary classification - mass regime
pub enum PlanetClass {
    Rocky,       // M < 2 M⊕
    Transitional, // 2-5 M⊕
    Volatile,    // 5-160 M⊕
    Giant,       // > 160 M⊕
}

// Observable expression (21 variants)
pub enum PlanetType {
    SubEarth,
    Barren,
    Lava { tidally_locked: bool, surface_temp_k: f64 },
    // ... etc
}

// Bulk composition as mass fractions
pub struct Composition {
    pub iron: f64,
    pub silicate: f64,
    pub water: f64,
    pub h_he_gas: f64,  // Named to match existing codebase
}

// Complete planet
pub struct Planet {
    pub mass: Mass,
    pub radius: Length,
    pub semi_major_axis: Length,
    pub eccentricity: f64,
    pub inclination: f64,
    pub class: PlanetClass,
    pub planet_type: PlanetType,
    pub composition: Composition,
    pub equilibrium_temp: Temperature,
    pub surface_gravity: f64,
    pub escape_velocity: Velocity,
}

// System architecture classification
pub enum SystemArchitecture {
    CompactMulti,    // TRAPPIST-1 style
    Mixed,           // Solar System style
    GiantDominated,  // Hot/cold Jupiter systems
    Sparse,          // Few or no planets
}

// Complete system
pub struct PlanetarySystem {
    pub stellar_mass: f64,
    pub stellar_luminosity: f64,
    pub stellar_temperature: f64,
    pub stellar_metallicity: f64,
    pub spectral_type: String,
    pub planets: Vec<Planet>,
    pub architecture: SystemArchitecture,
}
```

### Generation Pipeline

```
1. Sample SystemArchitecture from stellar properties
       ↓
2. Based on architecture, determine planet count and zones
       ↓
3. For each planet:
   a. Sample mass from occurrence-weighted distribution
   b. Classify → PlanetClass
   c. Sample radius from M-R relation (with scatter)
   d. Sample orbital period, convert to semi-major axis
   e. Sample eccentricity, inclination
   f. Sample composition based on position vs. snow line
   g. Compute equilibrium temperature, incident flux
   h. Determine PlanetType from class + environment
       ↓
4. Sort planets by semi-major axis
       ↓
5. Check Hill stability (reject and retry if unstable)
       ↓
6. Return PlanetarySystem
```

---

## Usage

### Basic Generation

```rust
use rand::SeedableRng;
use rand_chacha::ChaChaRng;
use planetary::generate_planetary_system;

let mut rng = ChaChaRng::seed_from_u64(42);

// Sun-like star
let system = generate_planetary_system(
    &mut rng,
    1.0,    // stellar mass (M☉)
    1.0,    // stellar luminosity (L☉)  
    5778.0, // stellar temperature (K)
    0.0,    // stellar metallicity [Fe/H]
    "G",    // spectral type
);

println!("Architecture: {}", system.architecture);
println!("Planets: {}", system.planets.len());

for planet in &system.planets {
    println!(
        "  {}: {:.2} M⊕ at {:.2} AU ({})",
        planet.planet_type,
        planet.mass.to_earth_masses(),
        planet.semi_major_axis.to_au(),
        planet.class,
    );
}
```

### Convenience Functions

```rust
// Solar System analog
let system = solar_system_analog(&mut rng);

// M dwarf system (0.3 M☉)
let system = m_dwarf_system(&mut rng, 0.3);

// Force a hot Jupiter system
let system = hot_jupiter_system(&mut rng);
```

### Custom Configuration

```rust
use planetary::{generate_planetary_system_with_config, GenerationConfig};

let config = GenerationConfig {
    max_stability_attempts: 100,
    min_hill_separation: 10.0,  // Stricter stability
    enforce_stability: true,
};

let system = generate_planetary_system_with_config(
    &mut rng, 1.0, 1.0, 5778.0, 0.0, "G", &config
);
```

---

## Next Steps

### Phase 1: Core Refinement (Current)

- [x] Basic crate structure
- [x] PlanetClass and PlanetType enums
- [x] Planet struct with M-R relations
- [x] Composition system
- [x] Occurrence rate sampling
- [x] System architecture classification
- [x] Stability filtering
- [ ] **Integration with stellar_forge crate** (convenience methods)
- [ ] **Validation against Kepler statistics** (occurrence rates, period ratios)
- [ ] **Unit tests for edge cases** (very low/high metallicity, extreme masses)

### Phase 2: Enhanced Physics

- [ ] **Photoevaporation logic**: Explicit modeling of envelope stripping based on XUV flux and age
- [ ] **Resonant chain generation**: MMR period ratios for compact systems
- [ ] **Orbital element correlations**: e-i coupling, secular architecture
- [ ] **Tidal effects**: Circularization, spin-orbit coupling for close-in planets
- [ ] **Atmospheric escape rates**: Mass loss for hot Jupiters

### Phase 3: System Completeness

- [ ] **Moon systems**: Apply same classification to satellites (see plan doc)
- [ ] **Debris belts**: Asteroid belt / Kuiper belt analogs
- [ ] **Ring systems**: Probability and extent for giant planets
- [ ] **Binary star support**: S-type and P-type orbits

### Phase 4: Validation & Tuning

- [ ] **Population synthesis**: Generate 10^5 systems, compare to Kepler/TESS
- [ ] **Period ratio distribution**: Should show excess near, but not at, MMRs
- [ ] **Radius distribution**: Should reproduce Fulton gap
- [ ] **Giant planet eccentricities**: Should be higher than small planets
- [ ] **Hot Jupiter lonely phenomenon**: Systems with hot Jupiters lack nearby companions

### Phase 5: Integration

- [ ] **WASM bindings**: For web visualization
- [ ] **Stellar integration**: `generate_for_star(&MainSequenceStar)`
- [ ] **Emergent simulation handoff**: Use statistical generator for initial conditions

---

## Integration with Stellar Forge

The `planetary` crate is designed to integrate with the existing stellar generation code. A future integration might look like:

```rust
// In planetary crate, with `stellar` feature enabled
use stellar_forge::stellar::MainSequenceStar;

impl PlanetarySystem {
    pub fn generate_for_star(rng: &mut ChaChaRng, star: &MainSequenceStar) -> Self {
        generate_planetary_system(
            rng,
            star.mass.to_solar_masses(),
            star.luminosity,
            star.temperature.to_kelvin(),
            star.metallicity,
            &star.spectral_type.to_string(),
        )
    }
}

// Usage
let star = sample_main_sequence_star(&mut rng);
let system = PlanetarySystem::generate_for_star(&mut rng, &star);
```

---

## Relationship to Emergent Simulation

This statistical generator and the emergent simulation serve different purposes:

| Use Case | Recommended Approach |
|----------|---------------------|
| Generate 10,000 systems for population study | Statistical generator |
| Understand how Jupiter's migration affected Mars | Emergent simulation |
| Procedural generation for a game | Statistical generator |
| Study resonant chain formation | Emergent simulation |
| Quick visualization demo | Statistical generator |
| Validate against Armitage textbook | Emergent simulation |

The two can work together:

1. **Statistical → Emergent**: Use statistical generator to create initial conditions, then run emergent simulation to evolve and refine
2. **Emergent → Statistical**: Calibrate statistical occurrence rates against emergent simulation outputs
3. **Parallel validation**: Both should produce similar population statistics if tuned correctly

---

## References

### Observational Foundations

- **Fulton et al. (2017)** — "The California-Kepler Survey. III. A Gap in the Radius Distribution of Small Planets" — Discovery of the radius valley
- **Petigura et al. (2018)** — "The California-Kepler Survey. IV. Metal-rich Stars Host a Greater Diversity of Planets" — Occurrence rates by spectral type and metallicity
- **Zhu & Dong (2021)** — "Exoplanet Statistics and Theoretical Implications" — Architecture statistics review

### Mass-Radius Relations

- **Wolfgang, Rogers & Ford (2016)** — "Probabilistic Mass-Radius Relationship for Sub-Neptune-Sized Planets" — Bayesian M-R framework with intrinsic scatter
- **Chen & Kipping (2017)** — "Probabilistic Forecasting of the Masses and Radii of Other Worlds" — Extended M-R relations

### Photoevaporation & Atmospheric Evolution

- **Lopez & Fortney (2013)** — "The Role of Core Mass in Controlling Evaporation" — Evaporation thresholds
- **Owen & Wu (2017)** — "The Evaporation Valley in the Kepler Planets" — Theoretical explanation of radius valley

### Classification Systems

- **Durand-Manterola (2011)** — Three-class system based on pressure support physics (basis for our four-class system)

---

## Appendix: Mass-Radius Relations

The generator uses class-specific power laws with intrinsic scatter, plus a more sophisticated envelope structure model for sub-Neptunes.

### Simple Power Laws (no envelope)

### Rocky (M < 2 M⊕)

```
R = M^0.29 R⊕
σ_R ≈ 15%
```

Density increases with mass (self-compression of rock/iron).

### Transitional (2-5 M⊕)

```
R = 1.05 × M^0.40 R⊕
σ_R ≈ 20%
```

High scatter due to composition diversity (rocky vs. volatile-rich).

### Volatile (5-160 M⊕)

```
R = 0.895 × (M)^0.56 R⊕
σ_R ≈ 15%
```

Density decreases with mass (envelope accumulation dominates).

### Giant (> 160 M⊕)

```
R ≈ 11.54 × (M)^(-0.024) R⊕
σ_R ≈ 10%
```

Radius nearly constant (electron degeneracy). Jupiter and a 10× Jupiter have similar radii.

### Envelope Structure Model (Lopez & Fortney 2014)

For planets with H/He envelopes, we use a layered model:

```rust
// Core radius from rocky scaling
let r_core = 1.07 * core_mass_earth.powf(0.27);

// Envelope contribution
let envelope_scaling = 0.57;  // From structure models
let envelope_contribution = r_core * (M_env / M_core).powf(envelope_scaling);

// Total radius
let radius = r_core + envelope_contribution;
```

This produces more realistic sub-Neptune radii than simple power laws.

---

## Appendix: Occurrence Rate Calibration

Current occurrence rates (per star, P < 400 days, FGK stars):

| Size Bin | Radius Range | Occurrence |
|----------|--------------|------------|
| Sub-Earth | < 1 R⊕ | ~10% |
| Earth-sized | 1-1.5 R⊕ | ~20% |
| Super-Earth | 1.5-2.0 R⊕ | ~25% |
| Sub-Neptune | 2.0-4.0 R⊕ | ~30% |
| Neptune | 4-6 R⊕ | ~5% |
| Sub-Saturn | 6-10 R⊕ | ~2% |
| Jupiter | 10-15 R⊕ | ~1% |
| Hot Jupiter | >10 R⊕, P<10d | ~1% |

These rates are approximate and will be refined as the generator is validated against actual Kepler/TESS catalogs.
