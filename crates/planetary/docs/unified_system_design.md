# Unified Planetary System Design

## Overview

This document describes the design for a unified `PlanetarySystem` structure that serves as the common output format for all system generation approaches in Magrathea: statistical sampling, stellar-forge formation simulation, and manual construction.

## Problem Statement

The current `PlanetarySystem` struct duplicates stellar properties as scalar fields:

```rust
// Current implementation (system.rs)
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

This creates several issues:

1. **Duplicated data**: The `stellar` crate already has rich `StellarObject` types with mass, luminosity, temperature, metallicity, plus additional properties (age, variability, color, spectral subtype)
2. **Lost information**: Converting a `MainSequenceStar` to `PlanetarySystem` discards useful data
3. **Single-star assumption**: No support for binary or multiple star systems
4. **No provenance**: No tracking of how/when the system was generated
5. **Stringly-typed**: `spectral_type: String` loses type safety

## Proposed Design

### Core Structure

```rust
use stellar::StellarObject;
use serde::{Deserialize, Serialize};
use uuid::Uuid;

/// A complete planetary system with one or more stellar hosts
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PlanetarySystem {
    /// Stellar host(s) - at least one required
    ///
    /// For single-star systems, this contains one element.
    /// For binary systems, contains two elements (primary first).
    /// For hierarchical systems, elements are ordered by hierarchy.
    pub stars: Vec<StellarObject>,

    /// Planets in the system, sorted by semi-major axis
    pub planets: Vec<Planet>,

    /// System metadata (generation info, architecture)
    pub metadata: SystemMetadata,
}

/// Metadata about system generation and classification
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SystemMetadata {
    /// Unique identifier for this system (also used as RNG seed source)
    ///
    /// UUIDs are JSON-safe (serialized as strings) and avoid JavaScript's
    /// Number.MAX_SAFE_INTEGER limitation that corrupts large u64 values.
    /// The seed for RNG is derived from the UUID's bytes.
    pub id: Uuid,

    /// How this system was generated
    pub generation_method: GenerationMethod,

    /// System architecture classification
    pub architecture: SystemArchitecture,

    /// Optional proper name for "hero" systems (e.g., "New Eden", "Cygnus Prime")
    ///
    /// Most systems use only the auto-generated catalog_name() (e.g., "KV-4729").
    /// This field is for notable systems that deserve memorable names.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
}

impl SystemMetadata {
    /// Derive a u64 seed from the UUID for RNG initialization
    ///
    /// Uses the first 8 bytes of the UUID, providing deterministic
    /// seed generation from any UUID.
    pub fn seed(&self) -> u64 {
        self.id.as_u64_pair().0
    }

    /// Generate a short catalog designation from the UUID
    ///
    /// Format: Two uppercase letters + 4 digits (e.g., "KV-4729", "AN-0821")
    /// Deterministic - same UUID always produces same designation.
    /// Provides ~6.76 million unique combinations (26² × 10000).
    pub fn catalog_name(&self) -> String {
        let bytes = self.id.as_bytes();
        let prefix1 = (bytes[0] % 26 + b'A') as char;
        let prefix2 = (bytes[1] % 26 + b'A') as char;
        let number = u16::from_le_bytes([bytes[2], bytes[3]]) % 10000;
        format!("{}{}-{:04}", prefix1, prefix2, number)
    }

    /// Returns the display name: proper name if set, otherwise catalog name
    pub fn display_name(&self) -> String {
        self.name.clone().unwrap_or_else(|| self.catalog_name())
    }

    /// Create metadata with a random UUID
    pub fn new_random(generation_method: GenerationMethod, architecture: SystemArchitecture) -> Self {
        Self {
            id: Uuid::new_v4(),
            generation_method,
            architecture,
            name: None,
        }
    }

    /// Create metadata with a deterministic UUID derived from a seed string
    ///
    /// Useful for reproducible generation from human-readable identifiers.
    /// The same seed_name always produces the same UUID (and thus the same RNG seed).
    /// Note: This does NOT set the display name - use `with_name()` for that.
    pub fn from_seed_name(seed_name: &str, generation_method: GenerationMethod, architecture: SystemArchitecture) -> Self {
        Self {
            id: Uuid::new_v5(&Uuid::NAMESPACE_OID, seed_name.as_bytes()),
            generation_method,
            architecture,
            name: None,
        }
    }

    /// Set a proper name for this system (builder pattern)
    pub fn with_name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }
}

/// How the system was generated
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GenerationMethod {
    /// Fast occurrence-rate sampling (current `planetary` generator)
    Statistical,

    /// Physics-based formation simulation (`stellar-forge`)
    StellarForge,

    /// User-defined or imported from external source
    Manual,

    /// Real observed system (e.g., from exoplanet archive)
    Observed,
}
```

### Derived Properties via Methods

Instead of storing computed values, provide methods that derive from the stellar hosts:

```rust
impl PlanetarySystem {
    /// Returns the primary (most massive) star
    pub fn primary_star(&self) -> &StellarObject {
        &self.stars[0]
    }

    /// Total luminosity of all stellar components (L☉)
    ///
    /// Used for habitable zone calculations in multi-star systems.
    pub fn total_luminosity(&self) -> f64 {
        self.stars.iter().map(|s| s.luminosity()).sum()
    }

    /// Effective stellar mass for orbital dynamics (M☉)
    ///
    /// For single stars, returns the star's mass.
    /// For close binaries, returns combined mass.
    pub fn effective_mass(&self) -> f64 {
        self.stars.iter().map(|s| s.mass().to_solar_masses()).sum()
    }

    /// Primary star's metallicity [Fe/H]
    pub fn metallicity(&self) -> f64 {
        self.primary_star().metallicity()
    }

    /// Primary star's spectral type as string (e.g., "G2")
    pub fn spectral_type(&self) -> String {
        self.primary_star().spectral_type_string()
    }

    /// Whether this is a multi-star system
    pub fn is_binary(&self) -> bool {
        self.stars.len() > 1
    }

    /// Habitable zone boundaries based on total luminosity
    pub fn habitable_zone(&self) -> HabitableZone {
        HabitableZone::from_luminosity(self.total_luminosity())
    }

    /// Snow line distance in AU
    pub fn snow_line(&self) -> f64 {
        2.7 * self.total_luminosity().sqrt()
    }
}
```

### StellarObject Trait Extensions

The `StellarObject` enum needs methods to extract common properties:

```rust
// In stellar crate
impl StellarObject {
    /// Luminosity in solar luminosities
    pub fn luminosity(&self) -> f64 {
        match self {
            Self::MainSequence(s) => s.luminosity,
            Self::Giant(s) => s.luminosity,
            Self::WhiteDwarf(s) => s.luminosity,
            Self::NeutronStar(_) => 0.0,  // Negligible optical
            Self::BlackHole(_) => 0.0,     // No luminosity (unless accreting)
        }
    }

    /// Mass as units::Mass
    pub fn mass(&self) -> Mass {
        match self {
            Self::MainSequence(s) => s.mass,
            Self::Giant(s) => s.mass,
            Self::WhiteDwarf(s) => s.mass,
            Self::NeutronStar(s) => s.mass,
            Self::BlackHole(s) => s.mass,
        }
    }

    /// Effective temperature in Kelvin
    pub fn temperature(&self) -> f64 {
        match self {
            Self::MainSequence(s) => s.temperature.to_kelvin(),
            Self::Giant(s) => s.temperature.to_kelvin(),
            Self::WhiteDwarf(s) => s.temperature.to_kelvin(),
            Self::NeutronStar(_) => 0.0,  // Surface temp not meaningful
            Self::BlackHole(_) => 0.0,
        }
    }

    /// Metallicity [Fe/H], returns 0.0 for remnants
    pub fn metallicity(&self) -> f64 {
        match self {
            Self::MainSequence(s) => s.metallicity,
            Self::Giant(_) => 0.0,  // TODO: Add metallicity to GiantStar
            _ => 0.0,
        }
    }

    /// Spectral type as display string
    pub fn spectral_type_string(&self) -> String {
        match self {
            Self::MainSequence(s) => format!("{}{}", s.spectral_type, s.subtype),
            Self::Giant(s) => format!("{}{}", s.spectral_type, s.subtype),
            Self::WhiteDwarf(s) => format!("{:?}", s.spectral_type),
            Self::NeutronStar(_) => "NS".to_string(),
            Self::BlackHole(_) => "BH".to_string(),
        }
    }
}
```

## Multi-Star System Support

### Binary System Types

The unified structure naturally supports various binary configurations:

```rust
/// Binary system orbital configuration
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum BinaryConfiguration {
    /// Planets orbit individual stars (e.g., Alpha Centauri)
    SType {
        /// Semi-major axis of binary orbit (AU)
        binary_separation: f64,
        /// Eccentricity of binary orbit
        binary_eccentricity: f64,
    },

    /// Planets orbit both stars (e.g., Kepler-16)
    PType {
        /// Semi-major axis of binary orbit (AU)
        binary_separation: f64,
        /// Eccentricity of binary orbit
        binary_eccentricity: f64,
    },

    /// Single star system (default)
    Single,
}
```

This could be added to `SystemMetadata` when binary support is implemented.

### Example: Binary System

```rust
let kepler16_like = PlanetarySystem {
    stars: vec![
        StellarObject::MainSequence(k_dwarf),
        StellarObject::MainSequence(m_dwarf),
    ],
    planets: vec![circumbinary_planet],
    metadata: SystemMetadata::from_name(
        "Kepler-16-analog",
        GenerationMethod::Statistical,
        SystemArchitecture::Sparse,
    ),
};
```

## Migration Strategy

### Phase 1: Add New Types (Non-Breaking)

1. Create `SystemMetadata` and `GenerationMethod` types
2. Add accessor methods to `StellarObject` in stellar crate
3. Keep existing `PlanetarySystem` unchanged

### Phase 2: Create Unified Type

1. Create new `UnifiedPlanetarySystem` alongside existing type
2. Add `From<PlanetarySystem>` conversion
3. Update generators to produce unified type
4. Add `into_legacy()` method for backwards compatibility

### Phase 3: Replace Original

1. Rename `UnifiedPlanetarySystem` to `PlanetarySystem`
2. Remove legacy struct
3. Update all consumers

### Backwards Compatibility Bridge

```rust
impl PlanetarySystem {
    /// Convert to legacy format for compatibility
    pub fn to_legacy(&self) -> LegacyPlanetarySystem {
        let primary = self.primary_star();
        LegacyPlanetarySystem {
            stellar_mass: primary.mass().to_solar_masses(),
            stellar_luminosity: primary.luminosity(),
            stellar_temperature: primary.temperature(),
            stellar_metallicity: primary.metallicity(),
            spectral_type: primary.spectral_type_string(),
            planets: self.planets.clone(),
            architecture: self.metadata.architecture,
        }
    }
}
```

## Generator Updates

### Statistical Generator

```rust
pub fn generate_planetary_system(
    star: StellarObject,
    id: Uuid,
) -> PlanetarySystem {
    // Derive seed from UUID for RNG
    let seed = id.as_u64_pair().0;
    let mut rng = ChaChaRng::seed_from_u64(seed);

    let architecture = SystemArchitecture::sample(
        &mut rng,
        &star.spectral_type_string(),
        star.metallicity(),
    );

    let planets = generate_planets(&mut rng, &star, architecture);

    PlanetarySystem {
        stars: vec![star],
        planets,
        metadata: SystemMetadata {
            id,
            generation_method: GenerationMethod::Statistical,
            architecture,
        },
    }
}

// Convenience: generate with random UUID
pub fn generate_planetary_system_random(star: StellarObject) -> PlanetarySystem {
    generate_planetary_system(star, Uuid::new_v4())
}

// Convenience: generate with deterministic name
pub fn generate_planetary_system_named(star: StellarObject, name: &str) -> PlanetarySystem {
    let id = Uuid::new_v5(&Uuid::NAMESPACE_OID, name.as_bytes());
    generate_planetary_system(star, id)
}
```

### Stellar-Forge Integration

```rust
// stellar-forge output
pub fn simulate_formation(
    disk: ProtoplanetaryDisk,
    star: StellarObject,
    id: Uuid,
) -> PlanetarySystem {
    let seed = id.as_u64_pair().0;
    let mut rng = ChaChaRng::seed_from_u64(seed);

    let planets = run_simulation(&mut rng, disk, &star);
    let architecture = classify_architecture(&planets);

    PlanetarySystem {
        stars: vec![star],
        planets,
        metadata: SystemMetadata {
            id,
            generation_method: GenerationMethod::StellarForge,
            architecture,
        },
    }
}
```

## Serialization Example

```json
{
  "stars": [
    {
      "type": "MainSequence",
      "mass": { "value": 1.0, "unit": "solar_masses" },
      "luminosity": 1.0,
      "temperature": { "value": 5778, "unit": "kelvin" },
      "spectralType": "G",
      "subtype": 2,
      "metallicity": 0.0,
      "age": { "value": 4.6, "unit": "gyr" }
    }
  ],
  "planets": [
    {
      "mass": { "value": 1.0, "unit": "earth_masses" },
      "radius": { "value": 1.0, "unit": "earth_radii" },
      "semiMajorAxis": { "value": 1.0, "unit": "au" },
      "eccentricity": 0.017,
      "class": "Rocky",
      "planetType": "Terran"
    }
  ],
  "metadata": {
    "id": "f47ac10b-58cc-4372-a567-0e02b2c3d479",
    "generationMethod": "Statistical",
    "architecture": "Mixed"
  }
}
```

## Design Decisions (Resolved)

1. **Where should the unified type live?**
   - **Decision**: New `star-system` crate (Option B) - provides clean separation
   - The `star-system` crate depends on both `planetary` and `stellar`
   - The `system-generator` crate depends on `star-system` for output types

2. **Should `stars` be `NonEmpty<StellarObject>`?**
   - **Decision**: Use `Vec` with runtime assertion in constructor
   - `PlanetarySystem::new()` panics if `stars` is empty
   - Simple API, clear error at construction time

3. **Binary orbit parameters**
   - **Decision**: Deferred - single-star focus initially
   - `BinaryConfiguration` enum defined but not yet integrated
   - Can be added to `SystemMetadata` when binary support is implemented

4. **Planet-star assignment in binaries**
   - **Decision**: Deferred - all planets orbit primary/barycenter initially
   - Future: May add `host_star_index: Option<usize>` on `Planet`

## Current Crate Structure

```
crates/
├── planetary/              # Planet classification and characterization
│   └── src/
│       ├── lib.rs          # Re-exports: Planet, PlanetClass, PlanetType, Composition
│       ├── planet.rs       # Planet struct, HostStar
│       ├── planet_class.rs # PlanetClass enum (Rocky, Transitional, Volatile, Giant)
│       ├── planet_type.rs  # PlanetType enum (21+ variants)
│       └── composition.rs  # Composition struct
│
├── star-system/            # Unified system output types
│   └── src/
│       ├── lib.rs          # Re-exports all types
│       ├── system.rs       # PlanetarySystem, HabitableZone
│       ├── metadata.rs     # SystemMetadata, GenerationMethod
│       └── architecture.rs # SystemArchitecture enum
│
├── planetary-generator/    # Statistical system generation
│   └── src/
│       ├── lib.rs          # Re-exports generation functions
│       ├── generation.rs   # generate_planetary_system(), from_star()
│       └── sampling.rs     # Occurrence rate sampling functions
│
├── stellar/                # Stellar objects and classification
│   └── src/
│       ├── lib.rs          # Re-exports
│       └── stellar_objects.rs  # StellarObject enum with accessor methods
│
├── magrathea-wasm/         # WASM bindings for stellar-forge (disks, formation)
│   └── src/
│       ├── lib.rs          # Shared utilities (to_js, from_js)
│       ├── stellar.rs      # Stellar generation bindings
│       └── disk.rs         # Disk generation bindings
│
└── planetary-wasm/         # WASM bindings for statistical generation
    └── src/
        ├── lib.rs          # Shared utilities
        ├── stellar.rs      # Star generation for systems
        └── system.rs       # System generation and queries
```

## Implementation Plan

### Phase 1: Stellar Crate Updates

**Goal**: Add accessor methods to `StellarObject` so planetary crate can use it

**Files to modify**:
- `crates/stellar/src/stellar_objects.rs`

**Tasks**:
1. Add `luminosity()` method to `StellarObject`
2. Add `mass()` method to `StellarObject`
3. Add `temperature()` method to `StellarObject`
4. Add `metallicity()` method to `StellarObject`
5. Add `spectral_type_string()` method to `StellarObject`
6. Add `age()` method to `StellarObject` (returns `Option<Time>`)
7. Consider adding `metallicity` field to `GiantStar` if missing

**Tests**: Unit tests for each accessor across all stellar object variants

---

### Phase 2: New Types in Planetary Crate

**Goal**: Create `SystemMetadata`, `GenerationMethod`, and builder infrastructure

**Files to create/modify**:
- `crates/planetary/src/metadata.rs` (new)
- `crates/planetary/src/lib.rs` (add export)
- `crates/planetary/Cargo.toml` (add uuid dependency)

**New types**:
```rust
pub struct SystemMetadata { ... }
pub enum GenerationMethod { ... }
```

**Tasks**:
1. Add `uuid` crate dependency with `v4`, `v5`, and `serde` features
2. Create `metadata.rs` with `SystemMetadata` and `GenerationMethod`
3. Derive `Serialize`, `Deserialize`, `Clone`, `Debug`, `PartialEq`
4. Implement `SystemMetadata::seed()` to derive u64 from UUID
5. Implement `SystemMetadata::new_random()` and `from_name()` constructors
6. Export from `lib.rs`

---

### Phase 3: Unified PlanetarySystem

**Goal**: Replace scalar stellar fields with `Vec<StellarObject>`

**Files to modify**:
- `crates/planetary/src/system.rs`
- `crates/planetary/Cargo.toml` (add stellar dependency)

**Tasks**:
1. Add `stellar` crate as dependency
2. Replace struct fields with:
   ```rust
   pub stars: Vec<StellarObject>,
   pub planets: Vec<Planet>,
   pub metadata: SystemMetadata,
   ```
3. Update `PlanetarySystem::new()` constructor
4. Add validation: panic or error if `stars` is empty
5. Add `primary_star()`, `total_luminosity()`, `effective_mass()` methods
6. Add `metallicity()`, `spectral_type()` convenience methods
7. Update `is_stable()` to use `effective_mass()`
8. Update `habitable_zone_planets()` to use `total_luminosity()`
9. Keep `SystemArchitecture` in `metadata`

---

### Phase 4: Update Generation Functions

**Goal**: Update `generation.rs` to produce unified systems

**Files to modify**:
- `crates/planetary/src/generation.rs`

**Tasks**:
1. Update `generate_planetary_system()` signature:
   ```rust
   pub fn generate_planetary_system(
       star: StellarObject,
       id: Uuid,
   ) -> PlanetarySystem
   ```
2. RNG is now created internally from UUID-derived seed
3. Add `generate_planetary_system_random()` convenience function
4. Add `generate_planetary_system_named()` for deterministic generation from strings
5. Update `from_star()` to use new signature
6. Update all internal functions to extract stellar params from `StellarObject`
7. Set `generation_method: GenerationMethod::Statistical` in output

---

### Phase 5: Update Examples and Tests

**Goal**: Ensure all examples compile and tests pass

**Files to modify**:
- `crates/planetary/examples/*.rs`
- `crates/planetary/src/*.rs` (test modules)

**Tasks**:
1. Update `generate_systems.rs` example
2. Update `generate_planets.rs` example
3. Update `generate_planets_imf.rs` example
4. Update `generate_solar_systems.rs` example
5. Fix any broken tests
6. Add new tests for multi-star system construction
7. Add tests for metadata serialization

---

### Phase 6: Stellar-Forge Integration

**Goal**: Ensure stellar-forge can target the unified output

**Files to review**:
- `crates/stellar-forge/src/*.rs`

**Tasks**:
1. Review stellar-forge output format
2. Add conversion or direct construction of `PlanetarySystem`
3. Set `generation_method: GenerationMethod::StellarForge`
4. Document integration pattern

---

### Checklist

- [x] **Phase 1**: Stellar accessor methods ✅ (completed)
  - [x] `StellarObject::luminosity()`
  - [x] `StellarObject::mass()`
  - [x] `StellarObject::temperature()`
  - [x] `StellarObject::metallicity()`
  - [x] `StellarObject::spectral_type_string()`
  - [x] Unit tests

- [x] **Phase 2**: New metadata types ✅ (completed in `star-system` crate)
  - [x] Add `uuid` dependency with `v4`, `v5`, `serde` features
  - [x] `SystemMetadata` struct with `id: Uuid`, `name: Option<String>`
  - [x] `GenerationMethod` enum
  - [x] `SystemMetadata::seed()` - derive u64 from UUID
  - [x] `SystemMetadata::catalog_name()` - generate "XX-0000" format
  - [x] `SystemMetadata::display_name()` - proper name or catalog name
  - [x] `SystemMetadata::new_random()` constructor
  - [x] `SystemMetadata::from_seed_name()` constructor
  - [x] `SystemMetadata::with_name()` builder method
  - [x] Export from lib.rs

- [x] **Phase 3**: Unified PlanetarySystem ✅ (completed in `star-system` crate)
  - [x] Add stellar dependency
  - [x] Replace fields with `stars: Vec<StellarObject>`
  - [x] Add `metadata: SystemMetadata`
  - [x] Implement accessor methods (`primary_star`, `total_luminosity`, `effective_mass`, etc.)
  - [x] Update existing methods (`is_stable`, `habitable_zone_planets`, etc.)

- [x] **Phase 4**: Update generators ✅ (completed in `planetary-generator` crate)
  - [x] Update signature: `generate_planetary_system(star: StellarObject, id: Uuid)`
  - [x] Add `generate_planetary_system_random()`
  - [x] Add `generate_planetary_system_named()`
  - [x] Add `from_star()` and `from_star_with_id()` convenience functions
  - [x] Update internal helpers

- [x] **Phase 5**: Examples and tests ✅ (completed)
  - [x] Fix all examples
  - [x] Fix all tests
  - [x] Add UUID serialization tests

- [ ] **Phase 6**: Stellar-forge integration (deferred)
  - [ ] Review/update stellar-forge
  - [ ] Document integration

- [x] **Phase 7**: WASM bindings ✅ (completed in `planetary-wasm` crate)
  - [x] Create new `planetary-wasm` crate (separate from `magrathea-wasm`)
  - [x] System generation functions (`generate_system_*`)
  - [x] Batch generation (`generate_systems_batch`, `generate_solar_systems_batch`)
  - [x] System query functions (`system_habitable_zone`, `system_snow_line`, etc.)
  - [x] Stellar generation for custom stars

---

## References

- Current `PlanetarySystem`: `crates/planetary/src/system.rs`
- `StellarObject` definition: `crates/stellar/src/stellar_objects.rs`
- Statistical generator: `crates/planetary/src/generation.rs`
- Stellar-forge plans: `crates/stellar-forge/docs/`
