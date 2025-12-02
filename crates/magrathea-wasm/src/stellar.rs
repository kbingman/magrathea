//! WASM bindings for stellar object generation.

use rand::SeedableRng;
use rand_chacha::ChaChaRng;
use wasm_bindgen::prelude::*;

use stellar_forge::{
    giant_star as forge_giant_star, main_sequence_star as forge_main_sequence_star,
    neutron_star as forge_neutron_star, sample_main_sequence_star as forge_sample_main_sequence,
    sample_stellar_object as forge_sample_stellar_object, solar_analog as forge_solar_analog,
    white_dwarf as forge_white_dwarf,
};

use crate::to_js;

/// Generate a solar analog star.
///
/// Returns a MainSequenceStar with properties matching our Sun.
#[wasm_bindgen]
pub fn solar_analog() -> Result<JsValue, JsError> {
    to_js(&forge_solar_analog())
}

/// Generate a main sequence star with specified parameters.
///
/// # Arguments
/// * `mass_solar` - Stellar mass in solar masses (0.08 to 150)
/// * `metallicity` - Metallicity [Fe/H] in dex (0.0 = solar)
/// * `age_myr` - Stellar age in millions of years
#[wasm_bindgen]
pub fn main_sequence_star(
    mass_solar: f64,
    metallicity: f64,
    age_myr: f64,
) -> Result<JsValue, JsError> {
    to_js(&forge_main_sequence_star(mass_solar, metallicity, age_myr))
}

/// Sample a random main sequence star from the IMF distribution.
///
/// Uses the Kroupa IMF to sample a realistic stellar mass, then generates
/// the corresponding star.
///
/// # Arguments
/// * `seed` - Random seed for reproducible generation
#[wasm_bindgen]
pub fn sample_main_sequence_star(seed: u64) -> Result<JsValue, JsError> {
    let mut rng = ChaChaRng::seed_from_u64(seed);
    to_js(&forge_sample_main_sequence(&mut rng))
}

/// Generate a giant star from an initial mass.
///
/// # Arguments
/// * `seed` - Random seed for reproducible generation
/// * `initial_mass` - Initial stellar mass in solar masses
#[wasm_bindgen]
pub fn giant_star(seed: u64, initial_mass: f64) -> Result<JsValue, JsError> {
    let mut rng = ChaChaRng::seed_from_u64(seed);
    to_js(&forge_giant_star(&mut rng, initial_mass))
}

/// Generate a white dwarf with random properties.
///
/// # Arguments
/// * `seed` - Random seed for reproducible generation
#[wasm_bindgen]
pub fn white_dwarf(seed: u64) -> Result<JsValue, JsError> {
    let mut rng = ChaChaRng::seed_from_u64(seed);
    to_js(&forge_white_dwarf(&mut rng))
}

/// Generate a neutron star with random properties.
///
/// # Arguments
/// * `seed` - Random seed for reproducible generation
#[wasm_bindgen]
pub fn neutron_star(seed: u64) -> Result<JsValue, JsError> {
    let mut rng = ChaChaRng::seed_from_u64(seed);
    to_js(&forge_neutron_star(&mut rng))
}

/// Sample a random stellar object from realistic distributions.
///
/// Samples stellar type and mass from astrophysically motivated distributions.
///
/// # Arguments
/// * `seed` - Random seed for reproducible generation
#[wasm_bindgen]
pub fn sample_stellar_object(seed: u64) -> Result<JsValue, JsError> {
    let mut rng = ChaChaRng::seed_from_u64(seed);
    to_js(&forge_sample_stellar_object(&mut rng))
}
