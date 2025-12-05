//! WASM bindings for stellar object generation.
//!
//! Provides functions to create stars for use with system generation.

use rand::SeedableRng;
use rand_chacha::ChaChaRng;
use wasm_bindgen::prelude::*;

use stellar_forge::{
    sample_main_sequence_star as forge_sample_main_sequence, solar_analog as forge_solar_analog,
};

use crate::to_js;

/// Generate a solar analog star.
///
/// Returns a MainSequenceStar with properties matching our Sun.
#[wasm_bindgen]
pub fn solar_analog() -> Result<JsValue, JsError> {
    to_js(&forge_solar_analog())
}

/// Sample a random main sequence star from the IMF distribution.
///
/// Uses the Kroupa IMF to sample a realistic stellar mass, then generates
/// the corresponding star with random age and solar metallicity.
///
/// # Arguments
/// * `seed` - Random seed for reproducible generation
#[wasm_bindgen]
pub fn sample_main_sequence_star(seed: u64) -> Result<JsValue, JsError> {
    let mut rng = ChaChaRng::seed_from_u64(seed);
    to_js(&forge_sample_main_sequence(&mut rng))
}

/// Generate a main sequence star with specified mass.
///
/// Creates a star with the given mass, solar metallicity, and a middle-aged star.
///
/// # Arguments
/// * `mass_solar` - Stellar mass in solar masses (0.08 to ~100)
#[wasm_bindgen]
pub fn main_sequence_star_with_mass(mass_solar: f64) -> Result<JsValue, JsError> {
    // Use a reasonable default age (middle-aged for the mass)
    // More massive stars have shorter lifetimes, so scale accordingly
    let age_myr = if mass_solar > 2.0 {
        100.0 // Young for massive stars
    } else if mass_solar > 0.8 {
        4600.0 // Sun-like age
    } else {
        8000.0 // Older for low-mass stars (they live longer)
    };

    to_js(&stellar_forge::main_sequence_star(mass_solar, 0.0, age_myr))
}
