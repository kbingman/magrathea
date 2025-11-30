# Rust Coding Guidelines for Magrathea

This document outlines the Rust coding patterns and conventions used in the Magrathea codebase.

---

## Core Principles

1. **Idiomatic Rust** - Prefer Rust patterns over imperative styles from other languages
2. **Functional over Imperative** - Use iterators, `map`, `filter`, etc. over for loops when practical
3. **Pattern Matching** - Leverage Rust's powerful `match` expressions
4. **Explicitness** - Prefer clear, explicit code over clever tricks
5. **Documentation** - All public APIs must have doc comments with examples

---

## Functional Programming Patterns

### ✅ Prefer Iterators over For Loops

**Good:**
```rust
// Functional style with iterator chains
let total: f32 = planets
    .iter()
    .filter(|p| p.mass > 1.0)
    .map(|p| p.mass)
    .sum();

// Collecting results
let earth_like: Vec<_> = planets
    .iter()
    .filter(|p| p.is_potentially_habitable())
    .collect();

// Early return with iterator methods
let first_gas_giant = planets
    .iter()
    .find(|p| matches!(p.planetary_type, PlanetaryType::GasGiant));
```

**Avoid:**
```rust
// Imperative style with for loops
let mut total = 0.0;
for planet in &planets {
    if planet.mass > 1.0 {
        total += planet.mass;
    }
}

// Manual collection
let mut earth_like = Vec::new();
for planet in &planets {
    if planet.is_potentially_habitable() {
        earth_like.push(planet);
    }
}
```

**When For Loops Are OK:**
- Complex state mutations that don't map well to iterators
- Performance-critical tight loops where iterators add overhead
- When iterator chains become deeply nested and hard to read

---

## Pattern Matching

### ✅ Use `match` for Multiple Conditions

**Good:**
```rust
fn classify_temperature(temp: f32) -> TemperatureClass {
    match temp {
        t if t < 150.0 => TemperatureClass::Cold,
        t if t < 400.0 => TemperatureClass::Temperate,
        t if t < 1000.0 => TemperatureClass::Warm,
        t if t < 2000.0 => TemperatureClass::Hot,
        _ => TemperatureClass::UltraHot,
    }
}

// Pattern matching on tuples
match (mass, temperature) {
    (m, _) if m < 0.1 => Self::None,
    (m, t) if m < 0.5 && t < 250.0 => Self::ThinCO2,
    (m, t) if (250.0..=350.0).contains(&t) => Self::NitrogenOxygen,
    _ => Self::None,
}

// Matching on enums
match planet_type {
    PlanetaryType::Terrestrial => classify_terrestrial(mass, temp),
    PlanetaryType::GasGiant => classify_gas_giant(temp),
    PlanetaryType::OceanWorld => classify_ocean_world(temp),
    _ => AtmosphereType::Unknown,
}
```

**Avoid:**
```rust
// Long if-else chains
if temp < 150.0 {
    TemperatureClass::Cold
} else if temp < 400.0 {
    TemperatureClass::Temperate
} else if temp < 1000.0 {
    TemperatureClass::Warm
} else if temp < 2000.0 {
    TemperatureClass::Hot
} else {
    TemperatureClass::UltraHot
}
```

### ✅ Use `matches!` for Simple Checks

**Good:**
```rust
if matches!(planet.planetary_type, PlanetaryType::GasGiant | PlanetaryType::IceGiant) {
    // Handle giants
}

if matches!(temp_class, TemperatureClass::Temperate) {
    // Handle temperate
}
```

**Avoid:**
```rust
if planet.planetary_type == PlanetaryType::GasGiant
    || planet.planetary_type == PlanetaryType::IceGiant {
    // Handle giants
}
```

---

## Option and Result Handling

### ✅ Use Combinators Over Explicit Matching

**Good:**
```rust
// Option combinators
let name = planet.ocean_subtype
    .map(|subtype| format!("{} {}", temp_class, subtype))
    .unwrap_or_else(|| format!("{} {}", temp_class, planet_type));

// Early return with ?
fn get_detail(&self) -> Result<&CompositionDetail, String> {
    self.composition.detail
        .as_ref()
        .ok_or_else(|| "No detail available".to_string())
}

// if let for single pattern
if let Some(detail) = &planet.composition.detail {
    println!("Silicates: {}", detail.rock_breakdown.silicate_fraction);
}
```

**Avoid:**
```rust
// Verbose explicit matching for simple cases
let name = match planet.ocean_subtype {
    Some(subtype) => format!("{} {}", temp_class, subtype),
    None => format!("{} {}", temp_class, planet_type),
};
```

---

## Type System Usage

### ✅ Use Enums for State

**Good:**
```rust
#[derive(Clone, Debug, PartialEq)]
pub enum DifferentiationState {
    Undifferentiated,
    PartiallyDifferentiated,
    FullyDifferentiated,
}

// Type-safe, impossible states ruled out
pub enum VolcanismLevel {
    None,
    Ancient,
    Active,
    Hyperactive,
    Cryovolcanic,
}
```

**Avoid:**
```rust
// Stringly-typed or bool flags
pub struct Planet {
    pub volcanism_active: bool,
    pub volcanism_type: String, // "active", "cryovolcanic", etc.
}
```

### ✅ Use Newtype Pattern for Units

**Example from codebase:**
```rust
use units::mass::Mass;
use units::length::Length;

let mass = Mass::from_earth_masses(1.0);
let radius = Length::from_earth_radii(1.0);
```

**Benefits:**
- Type safety prevents unit confusion
- Self-documenting API
- Compile-time guarantees

---

## Serialization

### ✅ Derive Serde Traits on Public Types

**Standard pattern:**
```rust
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PlanetaryComposition {
    pub iron_mass_fraction: f32,
    pub rock_mass_fraction: f32,
    // ...
}
```

**All public API types should:**
- Derive `Serialize, Deserialize` for JSON export
- Derive `Clone, Debug` for ergonomics
- Derive `PartialEq` for testing

---

## Documentation Standards

### ✅ Document All Public APIs

**Required format:**
```rust
/// Brief one-line summary
///
/// Longer description explaining the concept, behavior, and any important
/// details. Can be multiple paragraphs.
///
/// # Arguments
/// * `mass` - Planet mass in Earth masses
/// * `temperature` - Equilibrium temperature in Kelvin
///
/// # Returns
/// Classification of the planet's type
///
/// # Example
/// ```rust
/// use planetary::PlanetaryType;
///
/// let planet_type = PlanetaryType::classify(1.0, 1.0, 288.0, 0.2);
/// assert_eq!(planet_type, PlanetaryType::Terrestrial);
/// ```
///
/// # References
/// - Seager et al. (2007) - "Mass-Radius Relationships for Solid Exoplanets"
pub fn classify(mass: f32, radius: f32, temperature: f32, volatiles: f32) -> Self {
    // ...
}
```

**Always include:**
- Brief description
- Arguments documentation
- Return value description
- Example showing actual usage
- Scientific references where applicable

---

## Testing Patterns

### ✅ Organize Tests by Module

**Pattern:**
```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_specific_behavior() {
        // Arrange
        let input = 1.0;

        // Act
        let result = calculate(input);

        // Assert
        assert_eq!(result, expected);
    }

    #[test]
    fn test_edge_case() {
        // ...
    }
}
```

### ✅ Use Descriptive Test Names

**Good:**
```rust
#[test]
fn test_rock_breakdown_sums_correctly() { }

#[test]
fn test_hot_terrestrial_low_carbonates() { }

#[test]
fn test_ice_giant_has_multiple_ice_types() { }
```

**Avoid:**
```rust
#[test]
fn test_1() { }

#[test]
fn test_rocks() { }
```

---

## Error Handling

### ✅ Use Result for Fallible Operations

**Good:**
```rust
pub fn generate_system(seed: u64, index: u32) -> Result<ValidatedSystem, String> {
    let config = SystemConfig::new(seed)?;
    let bodies = generate_bodies(&config)?;
    validate_system(bodies)
}
```

### ✅ Use Custom Error Types for Libraries

**Pattern:**
```rust
#[derive(Debug)]
pub enum PlanetError {
    InvalidMass(f32),
    InvalidRadius(f32),
    ComputationFailed(String),
}

impl std::fmt::Display for PlanetError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidMass(m) => write!(f, "Invalid mass: {}", m),
            Self::InvalidRadius(r) => write!(f, "Invalid radius: {}", r),
            Self::ComputationFailed(msg) => write!(f, "Computation failed: {}", msg),
        }
    }
}

impl std::error::Error for PlanetError {}
```

---

## Performance Considerations

### ✅ Avoid Unnecessary Allocations

**Good:**
```rust
// Borrow instead of clone when possible
fn calculate(&self, composition: &PlanetaryComposition) -> f32 {
    composition.iron_mass_fraction * 2.0
}

// Use references in collections
let results: Vec<&Planet> = planets.iter()
    .filter(|p| p.mass > 1.0)
    .collect();
```

**Avoid:**
```rust
// Unnecessary clones
fn calculate(&self, composition: PlanetaryComposition) -> f32 {
    composition.iron_mass_fraction * 2.0
}
```

### ✅ Use const for Constants

**Good:**
```rust
const EARTH_MASS_KG: f64 = 5.972e24;
const GRAVITATIONAL_CONSTANT: f64 = 6.674e-11;
const BOLTZMANN_CONSTANT: f64 = 1.380649e-23;
```

---

## Code Organization

### ✅ Module Structure

```
crates/planetary/src/
├── lib.rs              // Public API, re-exports
├── planet.rs           // Main Planet struct
├── planetary_type.rs   // Classification enums
├── composition.rs      // Composition calculations
├── atmosphere.rs       // Atmosphere types
├── geology.rs          // Geological activity (new)
├── interior.rs         // Interior structure (new)
└── surface.rs          // Surface conditions (new)
```

**Guidelines:**
- One primary type per module
- Related functionality grouped together
- Re-export public API in `lib.rs`

### ✅ Import Organization

**Standard order:**
```rust
// 1. Standard library
use std::fmt;

// 2. External crates
use serde::{Deserialize, Serialize};

// 3. Internal crates (absolute paths)
use crate::planetary_type::PlanetaryType;
use crate::atmosphere::AtmosphereType;
```

---

## Common Anti-Patterns to Avoid

### ❌ Stringly-Typed APIs

**Avoid:**
```rust
pub fn classify(planet_type: &str) -> String { }
```

**Prefer:**
```rust
pub fn classify(planet_type: PlanetaryType) -> Classification { }
```

### ❌ Excessive `unwrap()` in Production Code

**Avoid:**
```rust
let value = some_option.unwrap(); // Panics if None!
```

**Prefer:**
```rust
let value = some_option.unwrap_or_default();
// or
let value = some_option.ok_or("Error message")?;
```

### ❌ Public Mutable State

**Avoid:**
```rust
pub struct BadPlanet {
    pub mass: f32, // Can be changed breaking invariants
}
```

**Prefer:**
```rust
pub struct GoodPlanet {
    mass: f32, // Private
}

impl GoodPlanet {
    pub fn mass(&self) -> f32 { self.mass }

    // Controlled mutation
    pub fn with_mass(mut self, mass: f32) -> Result<Self, PlanetError> {
        if mass <= 0.0 {
            return Err(PlanetError::InvalidMass(mass));
        }
        self.mass = mass;
        Ok(self)
    }
}
```

---

## Scientific Code Standards

### ✅ Always Include Units

**Good:**
```rust
/// Calculate surface gravity
///
/// # Arguments
/// * `mass` - Mass in Earth masses (M⊕)
/// * `radius` - Radius in Earth radii (R⊕)
///
/// # Returns
/// Surface gravity in m/s²
pub fn surface_gravity(mass: f32, radius: f32) -> f32 {
    9.81 * mass / (radius * radius)
}
```

### ✅ Cite References

**Good:**
```rust
/// Calculate magnetic field strength using dynamo theory
///
/// Based on scaling laws from Christensen (2010) and empirical
/// relationships for terrestrial planet magnetic fields.
///
/// # References
/// - Christensen (2010) - "Dynamo scaling laws and applications"
/// - Driscoll & Barnes (2015) - "Tidal effects on magnetic fields"
pub fn magnetic_field_strength(/* ... */) -> f32 {
    // ...
}
```

---

## Summary: Quick Reference

**When writing new code:**

1. ✅ Use iterators over for loops when practical
2. ✅ Use `match` for multiple conditions
3. ✅ Use builder pattern for optional expensive computations
4. ✅ Derive `Serialize, Deserialize, Clone, Debug, PartialEq` on public types
5. ✅ Document all public APIs with examples
6. ✅ Use enums for type-safe state
7. ✅ Include units in documentation
8. ✅ Cite scientific references
9. ✅ Use `Result` for fallible operations
10. ✅ Prefer borrowing over cloning

**When in doubt:**
- Look at existing code in the crate for patterns
- Prioritize clarity over cleverness
- Run `cargo clippy` and fix warnings
- Write tests demonstrating the API

---

## Before Committing

### ✅ Always Run `cargo fmt`

**The project uses automated formatting checks in CI/CD.** All code must be formatted with `rustfmt` before committing.

```bash
# Format all code in the workspace
cargo fmt

# Check formatting without making changes
cargo fmt -- --check
```

**Why this matters:**
- GitHub Actions will fail if code is not formatted
- Consistent formatting across the codebase
- Avoids unnecessary formatting diffs in PRs

**Recommended workflow:**
```bash
# Before committing
cargo fmt              # Format code
cargo clippy           # Check for warnings
cargo test             # Run tests
git add .
git commit -m "..."
```

**Pro tip:** Configure your editor to run `cargo fmt` on save:
- VS Code: Install `rust-analyzer` extension, enable format-on-save
- IntelliJ: Enable rustfmt in Preferences → Languages & Frameworks → Rust

---

## Git Commit Messages

Message:
  - Have a clear subject line summarizing the commit
  - Explain the project purpose (planetary formation generator)
  - List the key components being added
  - Use imperative mood ("Set up" not "Sets up")
  - Keep lines under 72 characters


**Last Updated:** November 29, 2025
