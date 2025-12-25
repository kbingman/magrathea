# Composition Module Migration Plan

## Overview

Migrate MK3's detailed `PlanetaryComposition` to replace main's simple `Composition`.

**Current state (main):**
```rust
struct Composition {
    iron: f64,
    silicate: f64,
    water: f64,
    h_he_gas: f64,
}
```

**Target state (from MK3):**
```rust
struct PlanetaryComposition {
    core: CoreComposition,      // total, iron, sulfides
    mantle: MantleComposition,  // total, silicates, carbonates
    ices: IceComposition,       // total, water_ice, ammonia_ice, methane_ice, co2_ice
    atmosphere: AtmosphereComposition,  // total, H2, He, N2, O2, CO2, CH4, H2O, other
}
```

---

## Phase 1: Dependencies

### 1.1 Port `variation.rs`

MK3's composition uses `PlanetaryVariation` for stochastic variation.

**Source:** `magrathea-mk3/crates/planetary/src/variation.rs`

**Key types:**
- `PlanetaryVariation` - seeded RNG wrapper for deterministic variation
- `vary_composition()` - applies bounded random variation to fractions

**Adaptation needed:**
- Change `f32` to `f64` to match main's conventions
- Update RNG usage if needed

**Estimated size:** ~100-150 lines

---

## Phase 2: Type Mapping Strategy

### Problem

MK3's composition generation is based on `PlanetaryType`:
```rust
enum PlanetaryType {
    DwarfPlanet,
    Terrestrial,
    SuperTerrestrial,
    OceanWorld,
    SubNeptune,
    IceGiant,
    GasGiant,
}
```

Main uses a two-tier system:
```rust
enum PlanetClass {      // Mass regime
    Compact,            // < 2 M⊕
    Transitional,       // 2-5 M⊕
    Volatile,           // 5-160 M⊕
    Giant,              // > 160 M⊕
}

enum PlanetType {       // Observable expression (21+ variants)
    SubEarth, Barren, Lava, Desert, Frozen, Oceanic, Terran, ...
    SuperTerran, WaterWorld, Hycean, MiniNeptune, ...
    IceGiant, HotNeptune, WarmNeptune, ...
    GasGiant, HotJupiter, PuffySaturn, ...
}
```

### Solution: Map main's types to composition categories

Create a `CompositionCategory` that maps from main's types:

```rust
enum CompositionCategory {
    DwarfPlanet,
    Terrestrial,
    SuperTerrestrial,
    OceanWorld,
    SubNeptune,
    IceGiant,
    GasGiant,
}

impl CompositionCategory {
    fn from_planet(class: &PlanetClass, planet_type: &PlanetType) -> Self {
        match (class, planet_type) {
            // Compact class
            (PlanetClass::Compact, PlanetType::SubEarth) => Self::DwarfPlanet,
            (PlanetClass::Compact, PlanetType::DwarfPlanet) => Self::DwarfPlanet,
            (PlanetClass::Compact, PlanetType::KuiperBeltObject) => Self::DwarfPlanet,
            (PlanetClass::Compact, PlanetType::Oceanic) => Self::OceanWorld,
            (PlanetClass::Compact, _) => Self::Terrestrial,

            // Transitional class
            (PlanetClass::Transitional, PlanetType::SuperTerran) => Self::SuperTerrestrial,
            (PlanetClass::Transitional, PlanetType::WaterWorld) => Self::OceanWorld,
            (PlanetClass::Transitional, PlanetType::Hycean) => Self::OceanWorld,
            (PlanetClass::Transitional, PlanetType::MiniNeptune) => Self::SubNeptune,
            (PlanetClass::Transitional, _) => Self::SuperTerrestrial,

            // Volatile class
            (PlanetClass::Volatile, _) => Self::IceGiant,

            // Giant class
            (PlanetClass::Giant, _) => Self::GasGiant,
        }
    }
}
```

---

## Phase 3: Core Migration

### 3.1 Port layer structs

Port from MK3 with adaptations:

```rust
// All use f64 instead of f32

pub struct CoreComposition {
    pub total: f64,      // Fraction of total planet mass
    pub iron: f64,       // Fe/Ni (relative within core)
    pub sulfides: f64,   // FeS (relative within core)
}

pub struct MantleComposition {
    pub total: f64,
    pub silicates: f64,  // MgSiO3, Mg2SiO4
    pub carbonates: f64, // CaCO3, MgCO3
}

pub struct IceComposition {
    pub total: f64,
    pub water_ice: f64,
    pub ammonia_ice: f64,
    pub methane_ice: f64,
    pub co2_ice: f64,
}

pub struct AtmosphereComposition {
    pub total: f64,
    pub hydrogen: f64,
    pub helium: f64,
    pub nitrogen: f64,
    pub oxygen: f64,
    pub co2: f64,
    pub methane: f64,
    pub water_vapor: f64,
    pub other: f64,
}
```

### 3.2 Port composition logic

Key methods to port from MK3:

1. `from_bulk_fractions()` - Creates layered structure from bulk (iron, rock, water, gas)
2. `distribute_rock_relative()` - Splits rock into silicates/carbonates/sulfides
3. `distribute_water_relative()` - Splits ice into water/ammonia/methane/CO₂
4. `get_atmosphere_ratios()` - Maps AtmosphereType to gas species ratios
5. `apply_variation()` - Adds stochastic variation

### 3.3 Add backward-compatible API

For callers that just need bulk fractions:

```rust
impl Composition {
    /// Bulk fractions for simple use cases
    pub fn bulk(&self) -> BulkComposition {
        BulkComposition {
            iron: self.core.total * self.core.iron,
            rock: self.mantle.total,
            water: self.ices.total,
            gas: self.atmosphere.total,
        }
    }

    /// Legacy compatibility
    pub fn iron(&self) -> f64 { self.core.total }
    pub fn silicate(&self) -> f64 { self.mantle.total }
    pub fn water(&self) -> f64 { self.ices.total }
    pub fn h_he_gas(&self) -> f64 { self.atmosphere.total }
}
```

---

## Phase 4: Update Callers

### 4.1 Identify callers

Search for usages of `Composition`:
```bash
rg "Composition" --type rust crates/
```

Expected locations:
- `forge/src/generation.rs` - Planet generation
- `planetary/src/planet.rs` - Planet struct
- `planetary-wasm/` - WASM bindings

### 4.2 Update generation

In `forge`, composition generation needs to:
1. Determine `CompositionCategory` from class + type
2. Call new `Composition::new()` with appropriate params
3. Use `AtmosphereType` for gas ratios

### 4.3 Update WASM bindings

Ensure all new types have:
- `#[cfg_attr(feature = "tsify", derive(Tsify))]`
- Proper `#[serde(rename_all = "camelCase")]`

---

## Phase 5: Testing

### 5.1 Unit tests

Port MK3's composition tests with adaptations:
- Solar system examples (Earth, Mars, Jupiter, etc.)
- Layer fraction bounds (totals sum to ~1.0)
- Relative fractions within layers sum to 1.0
- Temperature-dependent ice distribution
- Atmosphere ratios match AtmosphereType

### 5.2 Integration tests

- Generate 1000 planets, verify all compositions valid
- Check no NaN or negative fractions
- Verify backward-compatible API works

---

## File Changes Summary

| File | Action |
|------|--------|
| `planetary/src/variation.rs` | Create (port from MK3) |
| `planetary/src/variation_test.rs` | Create |
| `planetary/src/composition.rs` | Replace (expand significantly) |
| `planetary/src/composition_test.rs` | Replace (expand significantly) |
| `planetary/src/lib.rs` | Add variation module, update exports |
| `forge/src/generation.rs` | Update composition generation |
| `planetary-wasm/src/*.rs` | Update if needed for new types |

---

## Estimated Effort

| Phase | Effort |
|-------|--------|
| Phase 1: Port variation.rs | Small (~1 hour) |
| Phase 2: Type mapping | Small (~30 min) |
| Phase 3: Core migration | Medium (~2-3 hours) |
| Phase 4: Update callers | Medium (~1-2 hours) |
| Phase 5: Testing | Medium (~1-2 hours) |
| **Total** | **~6-8 hours** |

---

## Risks & Mitigations

| Risk | Mitigation |
|------|------------|
| Breaking API changes | Add backward-compatible accessors |
| Composition fractions don't sum to 1.0 | Normalize in constructors, validate in tests |
| WASM binding issues | Test TypeScript output early |
| Performance regression from complexity | Profile if needed; composition is generated once per planet |

---

## Success Criteria

- [ ] All existing tests pass
- [ ] New composition tests pass (solar system examples)
- [ ] Bulk fractions sum to ~1.0 (within rounding)
- [ ] Relative fractions within each layer sum to 1.0
- [ ] Atmosphere composition matches AtmosphereType
- [ ] WASM bindings work correctly
- [ ] No clippy warnings
