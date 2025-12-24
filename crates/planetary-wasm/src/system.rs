//! WASM bindings for planetary system generation.
//!
//! Provides functions to generate complete planetary systems using
//! occurrence-rate based statistical sampling.

use rand::SeedableRng;
use rand_chacha::ChaChaRng;
use uuid::Uuid;
use wasm_bindgen::prelude::*;

use celestial::PlanetarySystem;
use forge::{
    from_star, from_star_with_id, generate_planetary_system, generate_planetary_system_named,
};
use planetary::Planet;
use protodisk::{MainSequenceStar, sample_main_sequence_star, solar_analog};
use stellar::StellarObject;

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
#[wasm_bindgen(js_name = "generateSystemForSolarAnalog")]
pub fn generate_system_for_solar_analog() -> PlanetarySystem {
    let star = StellarObject::MainSequence(solar_analog());

    generate_planetary_system(star, Uuid::new_v4())
}

/// Generate a planetary system around a random star.
///
/// Samples a star from the IMF distribution and generates a planetary
/// system around it. The seed controls both star and system generation.
///
/// # Arguments
/// * `seed` - Random seed for reproducible generation
#[wasm_bindgen(js_name = "generateSystemRandom")]
pub fn generate_system_random(seed: u64) -> PlanetarySystem {
    let mut rng = ChaChaRng::seed_from_u64(seed);
    let star = sample_main_sequence_star(&mut rng);

    from_star(&star)
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
/// const system1 = generateSystemFromSeed("my-system");
/// const system2 = generateSystemFromSeed("my-system");
/// // system1 and system2 are identical
/// ```
#[wasm_bindgen(js_name = "generateSystemFromSeed")]
pub fn generate_system_from_seed(name: &str) -> PlanetarySystem {
    let star = StellarObject::MainSequence(solar_analog());

    generate_planetary_system_named(star, name)
}

/// Generate a planetary system from a UUID string.
///
/// Allows precise control over the system's identity and RNG seed.
/// Uses a solar analog star.
///
/// # Arguments
/// * `uuid_str` - A valid UUID string (e.g., "f47ac10b-58cc-4372-a567-0e02b2c3d479")
#[wasm_bindgen(js_name = "generateSystemFromUuid")]
pub fn generate_system_from_uuid(
    #[wasm_bindgen(js_name = "uuidStr")] uuid_str: &str,
) -> Result<PlanetarySystem, JsError> {
    let id =
        Uuid::parse_str(uuid_str).map_err(|e| JsError::new(&format!("Invalid UUID: {}", e)))?;
    let star = StellarObject::MainSequence(solar_analog());

    Ok(generate_planetary_system(star, id))
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
#[wasm_bindgen(js_name = "generateSystemForStar")]
pub fn generate_system_for_star(star: MainSequenceStar) -> PlanetarySystem {
    from_star(&star)
}

/// Generate a planetary system around a star with a specific UUID.
///
/// # Arguments
/// * `star` - A MainSequenceStar object
/// * `uuid_str` - UUID string for identification and RNG seeding
#[wasm_bindgen(js_name = "generateSystemForStarWithUuid")]
pub fn generate_system_for_star_with_uuid(
    star: MainSequenceStar,
    #[wasm_bindgen(js_name = "uuidStr")] uuid_str: &str,
) -> Result<PlanetarySystem, JsError> {
    let id =
        Uuid::parse_str(uuid_str).map_err(|e| JsError::new(&format!("Invalid UUID: {}", e)))?;

    Ok(from_star_with_id(&star, id))
}

/// Generate a planetary system around a star with a seed name.
///
/// # Arguments
/// * `star` - A MainSequenceStar object
/// * `seed_name` - String used to derive the system's UUID
#[wasm_bindgen(js_name = "generateSystemForStarNamed")]
pub fn generate_system_for_star_named(
    star: MainSequenceStar,
    #[wasm_bindgen(js_name = "seedName")] seed_name: &str,
) -> PlanetarySystem {
    let stellar_obj = StellarObject::MainSequence(star);

    generate_planetary_system_named(stellar_obj, seed_name)
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
/// * `seed_name` - String seed for reproducible batch generation
/// * `count` - Number of systems to generate
///
/// # Returns
/// Array of PlanetarySystem objects
#[wasm_bindgen(js_name = "generateSystemsBatch")]
pub fn generate_systems_batch(
    #[wasm_bindgen(js_name = "seedName")] seed_name: &str,
    count: u32,
) -> Vec<PlanetarySystem> {
    let base_id = Uuid::new_v5(&Uuid::NAMESPACE_OID, seed_name.as_bytes());
    let base_seed = u64::from_le_bytes(base_id.as_bytes()[..8].try_into().unwrap());

    (0..count)
        .map(|i| {
            let seed = base_seed.wrapping_add(i as u64);
            let mut rng = ChaChaRng::seed_from_u64(seed);
            let star = sample_main_sequence_star(&mut rng);

            // Create deterministic UUID from seed for reproducibility
            let id = Uuid::new_v5(&Uuid::NAMESPACE_OID, &seed.to_le_bytes());
            from_star_with_id(&star, id)
        })
        .collect()
}

/// Generate multiple systems around solar analog stars.
///
/// All systems have Sun-like host stars but different planetary configurations.
///
/// # Arguments
/// * `seed_name` - String seed for reproducible batch generation
/// * `count` - Number of systems to generate
#[wasm_bindgen(js_name = "generateSolarSystemsBatch")]
pub fn generate_solar_systems_batch(
    #[wasm_bindgen(js_name = "seedName")] seed_name: &str,
    count: u32,
) -> Vec<PlanetarySystem> {
    let star = solar_analog();
    let base_id = Uuid::new_v5(&Uuid::NAMESPACE_OID, seed_name.as_bytes());
    let base_seed = u64::from_le_bytes(base_id.as_bytes()[..8].try_into().unwrap());

    (0..count)
        .map(|i| {
            let seed = base_seed.wrapping_add(i as u64);
            // Create deterministic UUID from seed
            let id = Uuid::new_v5(&Uuid::NAMESPACE_OID, &seed.to_le_bytes());
            from_star_with_id(&star, id)
        })
        .collect()
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
#[wasm_bindgen(js_name = "systemHabitableZone")]
pub fn system_habitable_zone(system: &PlanetarySystem) -> HabitableZoneInfo {
    let hz = system.habitable_zone();

    HabitableZoneInfo {
        inner_edge_au: hz.inner_edge,
        outer_edge_au: hz.outer_edge,
    }
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
#[wasm_bindgen(js_name = "systemSnowLine")]
pub fn system_snow_line(system: &PlanetarySystem) -> f64 {
    system.snow_line()
}

/// Check if a system is dynamically stable.
///
/// Verifies that adjacent planet pairs have sufficient Hill sphere separation.
///
/// # Arguments
/// * `system` - A PlanetarySystem object
#[wasm_bindgen(js_name = "systemIsStable")]
pub fn system_is_stable(system: &PlanetarySystem) -> bool {
    system.is_stable()
}

/// Get planets in the habitable zone.
///
/// # Arguments
/// * `system` - A PlanetarySystem object
///
/// # Returns
/// Array of Planet objects within the habitable zone
#[wasm_bindgen(js_name = "systemHabitablePlanets")]
pub fn system_habitable_planets(system: &PlanetarySystem) -> Vec<Planet> {
    system
        .habitable_zone_planets()
        .into_iter()
        .cloned()
        .collect()
}

#[derive(serde::Serialize, tsify_next::Tsify)]
#[serde(rename_all = "camelCase")]
#[tsify(into_wasm_abi)]
pub struct HabitableZoneInfo {
    pub inner_edge_au: f64,
    pub outer_edge_au: f64,
}
