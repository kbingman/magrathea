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
| **Terrestrial** | M < 2 M⊕ | Self-compression | ρ increases with M, R ∝ M^0.29 |
| **Transitional** | 2–5 M⊕ | Envelope retention threshold | Stellar flux determines rocky vs. gaseous |
| **Volatile** | 5–160 M⊕ | Envelope accumulation | ρ *decreases* with M, R ∝ M^0.56 |
| **Degenerate** | > 160 M⊕ | Electron degeneracy | R nearly constant, ρ ∝ M |

The **Transitional** class is particularly interesting - it straddles the Fulton gap (radius valley). Whether a 3 M⊕ planet becomes a rocky Super-Terran or a puffy Mini-Neptune depends on stellar irradiation and photoevaporation history.

#### Tier 2: PlanetType (Observable Expression)

How a planet of a given class manifests based on temperature, composition, and stellar environment:

**Terrestrial expressions**: SubEarth, Barren, Lava, Desert, Frozen, Oceanic, Terran, Eyeball, Carbon, Iron

**Transitional expressions**: SuperTerran, MiniNeptune, WaterWorld, Hycean, Chthonian

**Volatile expressions**: SubNeptune, WarmNeptune, IceGiant

**Degenerate expressions**: GasGiant, HotJupiter, PuffySaturn

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

The system is split into three focused crates:

```
crates/
├── planetary/              # Planet classification (no generation logic)
│   ├── Cargo.toml
│   └── src/
│       ├── lib.rs          # Re-exports
│       ├── planet_class.rs # PlanetClass enum (Rocky, Transitional, Volatile, Giant)
│       ├── planet_type.rs  # PlanetType enum (21 expressions)
│       ├── planet.rs       # Planet struct, M-R relations, HostStar
│       └── composition.rs  # Composition (iron/silicate/water/H-He)
│
├── celestial/            # Unified output types for all backends
│   ├── Cargo.toml
│   └── src/
│       ├── lib.rs          # Re-exports
│       ├── system.rs       # PlanetarySystem (Vec<StellarObject>, Vec<Planet>)
│       ├── metadata.rs     # SystemMetadata (UUID, GenerationMethod)
│       └── architecture.rs # SystemArchitecture enum
│
└── planetary-generator/    # Statistical generation (one backend)
    ├── Cargo.toml
    └── src/
        ├── lib.rs          # Re-exports generation functions
        ├── generation.rs   # Main pipeline: generate_planetary_system()
        └── sampling.rs     # Occurrence rates, period/mass sampling
```

### Dependencies

```toml
# planetary/Cargo.toml
[dependencies]
units = { path = "../units" }
serde = { version = "1.0", features = ["derive"] }
rand = "0.9"
rand_chacha = "0.9"

# celestial/Cargo.toml
[dependencies]
planetary = { path = "../planetary" }
stellar = { path = "../stellar" }
units = { path = "../units" }
uuid = { version = "1.11", features = ["v4", "v5", "serde"] }
serde = "1.0"

# planetary-generator/Cargo.toml
[dependencies]
planetary = { path = "../planetary" }
celestial = { path = "../celestial" }
stellar = { path = "../stellar" }
units = { path = "../units" }
uuid = { version = "1.11", features = ["v4", "v5"] }
rand = "0.9"
rand_chacha = "0.9"
```

### Core Types

```rust
// Primary classification (in planetary crate)
pub enum PlanetClass {
    Rocky,         // M < 2 M⊕ (formerly Terrestrial)
    Transitional,  // 2-10 M⊕
    Volatile,      // 10-160 M⊕
    Giant,         // > 160 M⊕ (formerly Degenerate)
}

// Observable expression (21+ variants, in planetary crate)
pub enum PlanetType {
    SubTerran,
    Terran,
    SuperTerran,
    MiniNeptune,
    // ... etc
}

// Bulk composition as mass fractions (in planetary crate)
pub struct Composition {
    pub iron_mass_fraction: f64,
    pub rock_mass_fraction: f64,
    pub water_mass_fraction: f64,
    pub gas_mass_fraction: f64,
    pub detail: Option<CompositionDetail>,
}

// Complete planet (in planetary crate)
pub struct Planet {
    pub mass: Mass,
    pub radius: Length,
    pub semi_major_axis: Length,
    pub eccentricity: f64,
    pub inclination: f64,
    pub class: PlanetClass,
    pub planet_type: PlanetType,
    pub composition: Composition,
    pub equilibrium_temp: f64,  // Kelvin
    pub surface_gravity: f64,   // m/s²
    pub escape_velocity: f64,   // m/s
}

// System architecture classification (in celestial crate)
pub enum SystemArchitecture {
    CompactMulti,    // TRAPPIST-1 style
    Mixed,           // Solar System style
    GiantDominated,  // Hot/cold Jupiter systems
    Sparse,          // Few or no planets
}

// System metadata (in celestial crate)
pub struct SystemMetadata {
    pub id: Uuid,
    pub generation_method: GenerationMethod,
    pub architecture: SystemArchitecture,
    pub name: Option<String>,
}

// Complete system (in celestial crate)
pub struct PlanetarySystem {
    pub stars: Vec<StellarObject>,  // At least one star
    pub planets: Vec<Planet>,       // Sorted by semi-major axis
    pub metadata: SystemMetadata,
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
use stellar_forge::solar_analog;
use forge::{generate_planetary_system, from_star};
use uuid::Uuid;

// From a StellarObject with specific UUID
let star = stellar::StellarObject::MainSequence(solar_analog());
let system = generate_planetary_system(star, Uuid::new_v4());

println!("System: {}", system.metadata.catalog_name());
println!("Architecture: {:?}", system.metadata.architecture);
println!("Planets: {}", system.planets.len());

for planet in &system.planets {
    println!(
        "  {:?}: {:.2} M⊕ at {:.2} AU ({:?})",
        planet.planet_type,
        planet.mass.to_earth_masses(),
        planet.semi_major_axis.to_au(),
        planet.class,
    );
}

// From a MainSequenceStar (convenience function, random UUID)
let star = stellar_forge::sample_main_sequence_star(&mut rng);
let system = from_star(&star);
```

### Reproducible Generation

```rust
use forge::generate_planetary_system_named;
use protodisk::solar_analog;
use stellar::StellarObject;

// Same name always produces same system
let star = StellarObject::MainSequence(solar_analog());
let system1 = generate_planetary_system_named(star.clone(), "my-system-42");
let system2 = generate_planetary_system_named(star, "my-system-42");

assert_eq!(system1.planets.len(), system2.planets.len());
assert_eq!(system1.metadata.id, system2.metadata.id);
```

### Population Studies

```rust
use rand::SeedableRng;
use rand_chacha::ChaChaRng;
use protodisk::sample_main_sequence_star;
use forge::from_star;

let mut rng = ChaChaRng::seed_from_u64(42);

// Generate 1000 systems for statistical analysis
let systems: Vec<_> = (0..1000)
    .map(|_| {
        let star = sample_main_sequence_star(&mut rng);
        from_star(&star)
    })
    .collect();

// Count architectures
let compact = systems.iter()
    .filter(|s| matches!(s.metadata.architecture, SystemArchitecture::CompactMulti))
    .count();
println!("Compact multi-planet systems: {}/1000", compact);
```

---

## Next Steps

### Observational Bias Philosophy

**Important**: Our goal is NOT to reproduce Kepler/TESS occurrence rates exactly. Those surveys have massive observational biases:

- **Transit geometry**: Only ~1% of planets at 1 AU transit their star
- **Size bias**: Small planets produce shallow, hard-to-detect transits
- **Period bias**: Need multiple transits, favoring short periods
- **Host star bias**: Active stars hide planetary signals

Kepler occurrence rates are **lower bounds**. The true population is likely 2-10× richer, especially for small and long-period planets. Our goal is to generate what's *actually there*, not what we can detect. Multi-planet systems are probably the norm. The Solar System is probably typical, not special.

### Phase 1: Core Implementation ✅ (Complete)

- [x] Basic crate structure (split into `planetary`, `celestial`, `planetary-generator`)
- [x] PlanetClass and PlanetType enums
- [x] Planet struct with M-R relations
- [x] Composition system with detailed breakdowns
- [x] Occurrence rate sampling (inner system via Kepler + outer system via RV/microlensing)
- [x] System architecture classification (CompactMulti, Mixed, GiantDominated, Sparse)
- [x] Stability filtering (Hill criterion)
- [x] **Stellar integration**: `from_star(&MainSequenceStar)`, `generate_planetary_system(StellarObject, Uuid)`
- [x] UUID-based identification and RNG seeding
- [x] TNO/Kuiper Belt object generation
- [x] **Validation tests** (sanity checks, not Kepler-matching)
- [ ] **Unit tests for edge cases** (very low/high metallicity, extreme masses)

### Phase 2: Enhanced Occurrence Rates ✅ (Complete)

- [x] **Boost planet frequency**: Increased to ~5.8 planets/system (was ~2.5)
- [x] **Ensure most stars have planets**: Sparse systems now ~7% (was ~45%)
- [x] **Rich outer systems**: Cold giants 25%, ice giants 60%
- [x] **Compact multi-planet systems**: CompactMulti dominant for M/K dwarfs

### Phase 3: Enhanced Physics

- [x] **Photoevaporation logic**: Mass-dependent envelope stripping (Owen & Wu 2017)
- [ ] **Resonant chain generation**: MMR period ratios for compact systems
- [ ] **Orbital element correlations**: e-i coupling, secular architecture
- [ ] **Tidal effects**: Circularization, spin-orbit coupling for close-in planets
- [ ] **Atmospheric escape rates**: Mass loss for hot Jupiters

### Phase 4: System Completeness

- [ ] **Moon systems**: Apply same classification to satellites (see plan doc)
- [ ] **Debris belts**: Asteroid belt / Kuiper belt analogs
- [ ] **Ring systems**: Probability and extent for giant planets
- [ ] **Binary star support**: S-type and P-type orbits ✅ (implemented)

### Phase 5: Validation & Tuning

- [x] **Population synthesis tests**: Sanity checks for plausibility
- [x] **Period ratio distribution**: Median ~1.27, MMR proximity present
- [x] **Radius distribution**: Fulton gap via mass-dependent photoevaporation
- [x] **Giant planet eccentricities**: Higher than small planets ✅
- [x] **Hot Jupiter lonely phenomenon**: 53× fewer companions ✅

### Phase 5: Integration ✅ (Complete)

- [x] **WASM bindings**: `planetary-wasm` crate for web visualization
- [x] **Stellar integration**: `from_star(&MainSequenceStar)` and `generate_planetary_system(StellarObject, Uuid)`
- [ ] **Emergent simulation handoff**: Use statistical generator for initial conditions

---

## Integration with Stellar Forge

The `planetary-generator` crate integrates directly with the stellar generation code via `StellarObject`:

```rust
use protodisk::{sample_main_sequence_star, solar_analog};
use stellar::StellarObject;
use forge::{generate_planetary_system, from_star, from_star_with_id};
use uuid::Uuid;

// Method 1: From MainSequenceStar (convenience)
let mut rng = ChaChaRng::seed_from_u64(42);
let star = sample_main_sequence_star(&mut rng);
let system = from_star(&star);  // Random UUID

// Method 2: With specific UUID for reproducibility
let star = sample_main_sequence_star(&mut rng);
let id = Uuid::new_v4();
let system = from_star_with_id(&star, id);

// Method 3: From any StellarObject (supports giants, white dwarfs, etc.)
let star = StellarObject::MainSequence(solar_analog());
let system = generate_planetary_system(star, Uuid::new_v4());

// The generation method is tracked in metadata
assert!(matches!(system.metadata.generation_method, GenerationMethod::Statistical));
```

### Multiple Backends

The `celestial` crate defines the unified output format, allowing multiple generation backends:

```rust
use celestial::{PlanetarySystem, GenerationMethod};

// Statistical generator (fast, occurrence-rate based)
let system = forge::from_star(&star);
assert_eq!(system.metadata.generation_method, GenerationMethod::Statistical);

// Future: Physics-based formation simulation
// let system = stellar_forge::simulate_formation(disk, star, id);
// assert_eq!(system.metadata.generation_method, GenerationMethod::StellarForge);

// Manual/imported systems
// let system = PlanetarySystem::new(stars, planets, metadata);
// metadata.generation_method = GenerationMethod::Manual;
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

The generator uses class-specific power laws with intrinsic scatter:

### Terrestrial (M < 2 M⊕)

```
R = M^0.29 R⊕
σ_R ≈ 15%
```

Density increases with mass (self-compression of rock/iron).

### Transitional (2-5 M⊕)

```
R = M^0.50 R⊕
σ_R ≈ 20%
```

High scatter due to composition diversity (rocky vs. volatile-rich).

### Volatile (5-160 M⊕)

```
R = 2.5 × (M/10)^0.56 R⊕
σ_R ≈ 15%
```

Density decreases with mass (envelope accumulation dominates).

### Degenerate (> 160 M⊕)

```
R ≈ 11.2 × (M/318)^(-0.04) R⊕
σ_R ≈ 10%
```

Radius nearly constant (electron degeneracy). Jupiter and a 10× Jupiter have similar radii.

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
