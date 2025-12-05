//! WASM bindings for planetary system generation.
//!
//! Provides functions to generate complete planetary systems using
//! occurrence-rate based statistical sampling.

use rand::SeedableRng;
use rand_chacha::ChaChaRng;
use uuid::Uuid;
use wasm_bindgen::prelude::*;

use planetary_generator::{
    from_star, from_star_with_id, generate_planetary_system, generate_planetary_system_named,
};
use star_system::PlanetarySystem;
use stellar::StellarObject;
use stellar_forge::{sample_main_sequence_star, solar_analog};

use crate::{from_js, to_js};

// =============================================================================
// System Generation - Convenience Functions
// =============================================================================

/// Generate a planetary system around a solar analog star.
///
/// Creates a Sun-like star and generates a statistically realistic
/// planetary system around it. Uses a random UUID for identification.
///
/// # Returns
/// A complete PlanetarySystem with star, planets, and metadata.
#[wasm_bindgen]
pub fn generate_system_for_solar_analog() -> Result<JsValue, JsError> {
    let star = StellarObject::MainSequence(solar_analog());
    let system = generate_planetary_system(star, Uuid::new_v4());
    to_js(&system)
}

/// Generate a planetary system around a random star.
///
/// Samples a star from the IMF distribution and generates a planetary
/// system around it. The seed controls both star and system generation.
///
/// # Arguments
/// * `seed` - Random seed for reproducible generation
#[wasm_bindgen]
pub fn generate_system_random(seed: u64) -> Result<JsValue, JsError> {
    let mut rng = ChaChaRng::seed_from_u64(seed);
    let star = sample_main_sequence_star(&mut rng);
    let system = from_star(&star);
    to_js(&system)
}

/// Generate a planetary system with a deterministic name-based seed.
///
/// The same name always produces the same system (with a solar analog star).
/// Useful for sharing specific systems or creating reproducible examples.
///
/// # Arguments
/// * `name` - A string used to derive the system's UUID and RNG seed
///
/// # Example (JavaScript)
/// ```javascript
/// const system1 = generate_system_from_seed("my-system");
/// const system2 = generate_system_from_seed("my-system");
/// // system1 and system2 are identical
/// ```
#[wasm_bindgen]
pub fn generate_system_from_seed(name: &str) -> Result<JsValue, JsError> {
    let star = StellarObject::MainSequence(solar_analog());
    let system = generate_planetary_system_named(star, name);
    to_js(&system)
}

/// Generate a planetary system from a UUID string.
///
/// Allows precise control over the system's identity and RNG seed.
/// Uses a solar analog star.
///
/// # Arguments
/// * `uuid_str` - A valid UUID string (e.g., "f47ac10b-58cc-4372-a567-0e02b2c3d479")
#[wasm_bindgen]
pub fn generate_system_from_uuid(uuid_str: &str) -> Result<JsValue, JsError> {
    let id =
        Uuid::parse_str(uuid_str).map_err(|e| JsError::new(&format!("Invalid UUID: {}", e)))?;
    let star = StellarObject::MainSequence(solar_analog());
    let system = generate_planetary_system(star, id);
    to_js(&system)
}

// =============================================================================
// System Generation - With Custom Stars
// =============================================================================

/// Generate a planetary system around a provided star.
///
/// Takes a MainSequenceStar (from stellar generation functions) and
/// generates a planetary system around it with a random UUID.
///
/// # Arguments
/// * `star` - A MainSequenceStar object from stellar generation functions
#[wasm_bindgen]
pub fn generate_system_for_star(star: JsValue) -> Result<JsValue, JsError> {
    let star: stellar_forge::MainSequenceStar = from_js(star)?;
    let system = from_star(&star);
    to_js(&system)
}

/// Generate a planetary system around a star with a specific UUID.
///
/// # Arguments
/// * `star` - A MainSequenceStar object
/// * `uuid_str` - UUID string for identification and RNG seeding
#[wasm_bindgen]
pub fn generate_system_for_star_with_uuid(
    star: JsValue,
    uuid_str: &str,
) -> Result<JsValue, JsError> {
    let star: stellar_forge::MainSequenceStar = from_js(star)?;
    let id =
        Uuid::parse_str(uuid_str).map_err(|e| JsError::new(&format!("Invalid UUID: {}", e)))?;
    let system = from_star_with_id(&star, id);
    to_js(&system)
}

/// Generate a planetary system around a star with a seed name.
///
/// # Arguments
/// * `star` - A MainSequenceStar object
/// * `seed_name` - String used to derive the system's UUID
#[wasm_bindgen]
pub fn generate_system_for_star_named(star: JsValue, seed_name: &str) -> Result<JsValue, JsError> {
    let star: stellar_forge::MainSequenceStar = from_js(star)?;
    let stellar_obj = StellarObject::MainSequence(star);
    let system = generate_planetary_system_named(stellar_obj, seed_name);
    to_js(&system)
}

// =============================================================================
// Batch Generation
// =============================================================================

/// Generate multiple planetary systems in batch.
///
/// Efficiently generates many systems for population studies or visualization.
/// Each system gets a unique star sampled from the IMF.
///
/// # Arguments
/// * `count` - Number of systems to generate
/// * `base_seed` - Base seed (each system uses base_seed + index)
///
/// # Returns
/// Array of PlanetarySystem objects
#[wasm_bindgen]
pub fn generate_systems_batch(count: u32, base_seed: u64) -> Result<JsValue, JsError> {
    let systems: Vec<PlanetarySystem> = (0..count)
        .map(|i| {
            let seed = base_seed.wrapping_add(i as u64);
            let mut rng = ChaChaRng::seed_from_u64(seed);
            let star = sample_main_sequence_star(&mut rng);
            from_star(&star)
        })
        .collect();

    to_js(&systems)
}

/// Generate multiple systems around solar analog stars.
///
/// All systems have Sun-like host stars but different planetary configurations.
///
/// # Arguments
/// * `count` - Number of systems to generate
/// * `base_seed` - Base seed for reproducibility
#[wasm_bindgen]
pub fn generate_solar_systems_batch(count: u32, base_seed: u64) -> Result<JsValue, JsError> {
    let star = solar_analog();

    let systems: Vec<PlanetarySystem> = (0..count)
        .map(|i| {
            let seed = base_seed.wrapping_add(i as u64);
            // Create deterministic UUID from seed
            let id = Uuid::new_v5(&Uuid::NAMESPACE_OID, &seed.to_le_bytes());
            from_star_with_id(&star, id)
        })
        .collect();

    to_js(&systems)
}

// =============================================================================
// System Queries
// =============================================================================

/// Get the habitable zone boundaries for a system.
///
/// # Arguments
/// * `system` - A PlanetarySystem object
///
/// # Returns
/// Object with innerEdge and outerEdge in AU
#[wasm_bindgen]
pub fn system_habitable_zone(system: JsValue) -> Result<JsValue, JsError> {
    let system: PlanetarySystem = from_js(system)?;
    let hz = system.habitable_zone();

    to_js(&HabitableZoneInfo {
        inner_edge_au: hz.inner_edge,
        outer_edge_au: hz.outer_edge,
    })
}

/// Get the snow line distance for a system.
///
/// The snow line is where water ice can condense, affecting planet composition.
///
/// # Arguments
/// * `system` - A PlanetarySystem object
///
/// # Returns
/// Snow line distance in AU
#[wasm_bindgen]
pub fn system_snow_line(system: JsValue) -> Result<f64, JsError> {
    let system: PlanetarySystem = from_js(system)?;
    Ok(system.snow_line())
}

/// Check if a system is dynamically stable.
///
/// Verifies that adjacent planet pairs have sufficient Hill sphere separation.
///
/// # Arguments
/// * `system` - A PlanetarySystem object
#[wasm_bindgen]
pub fn system_is_stable(system: JsValue) -> Result<bool, JsError> {
    let system: PlanetarySystem = from_js(system)?;
    Ok(system.is_stable())
}

/// Get planets in the habitable zone.
///
/// # Arguments
/// * `system` - A PlanetarySystem object
///
/// # Returns
/// Array of Planet objects within the habitable zone
#[wasm_bindgen]
pub fn system_habitable_planets(system: JsValue) -> Result<JsValue, JsError> {
    let system: PlanetarySystem = from_js(system)?;
    let hz_planets: Vec<_> = system
        .habitable_zone_planets()
        .into_iter()
        .cloned()
        .collect();
    to_js(&hz_planets)
}

#[derive(serde::Serialize)]
#[serde(rename_all = "camelCase")]
struct HabitableZoneInfo {
    inner_edge_au: f64,
    outer_edge_au: f64,
}
