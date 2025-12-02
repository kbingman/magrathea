//! WASM bindings for Magrathea stellar and disk generation.
//!
//! This crate provides JavaScript/TypeScript bindings for the core Magrathea
//! library using `wasm-bindgen` and `serde-wasm-bindgen` for seamless type
//! conversion.

use rand::SeedableRng;
use rand_chacha::ChaChaRng;
use wasm_bindgen::prelude::*;

use stellar_forge::{
    DiskModel, GasDisk, MainSequenceStar, giant_star as forge_giant_star,
    main_sequence_star as forge_main_sequence_star, neutron_star as forge_neutron_star,
    sample_main_sequence_star as forge_sample_main_sequence,
    sample_stellar_object as forge_sample_stellar_object, solar_analog as forge_solar_analog,
    white_dwarf as forge_white_dwarf,
};

fn to_js<T: serde::Serialize>(value: &T) -> Result<JsValue, JsError> {
    serde_wasm_bindgen::to_value(value).map_err(|e| JsError::new(&e.to_string()))
}

fn from_js<T: serde::de::DeserializeOwned>(value: JsValue) -> Result<T, JsError> {
    serde_wasm_bindgen::from_value(value).map_err(|e| JsError::new(&e.to_string()))
}

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

/// Create a Minimum Mass Solar Nebula disk.
///
/// Returns a GasDisk with MMSN parameters around a solar-mass star.
#[wasm_bindgen]
pub fn mmsn_disk() -> Result<JsValue, JsError> {
    to_js(&GasDisk::mmsn())
}

/// Create a gas disk scaled to a star's properties.
///
/// # Arguments
/// * `star` - A MainSequenceStar object (as returned by other functions)
#[wasm_bindgen]
pub fn gas_disk_for_star(star: JsValue) -> Result<JsValue, JsError> {
    let star: MainSequenceStar = from_js(star)?;
    to_js(&GasDisk::for_star(&star))
}

/// Query disk properties at a given radius.
///
/// # Arguments
/// * `disk` - A GasDisk object
/// * `radius_au` - Radius in AU
///
/// # Returns
/// Object with surfaceDensity, temperature, scaleHeight, etc.
#[wasm_bindgen]
pub fn disk_properties_at(disk: JsValue, radius_au: f64) -> Result<JsValue, JsError> {
    let disk: GasDisk = from_js(disk)?;
    let r = units::Length::from_au(radius_au);

    to_js(&DiskProperties {
        radius_au,
        surface_density_g_cm2: disk.surface_density(r).to_grams_per_cm2(),
        temperature_k: disk.temperature(r).to_kelvin(),
        scale_height_au: disk.scale_height(r).to_au(),
        aspect_ratio: disk.aspect_ratio(r),
        midplane_density_g_cm3: disk.midplane_density(r).to_grams_per_cm3(),
        pressure_dyn_cm2: disk.pressure(r).to_dyn_per_cm2(),
        sound_speed_cm_s: disk.sound_speed(r).to_cm_per_sec(),
        keplerian_velocity_cm_s: disk.keplerian_velocity(r).to_cm_per_sec(),
        orbital_period_years: disk.orbital_period(r).to_years(),
    })
}

/// Disk properties at a specific radius.
#[derive(serde::Serialize)]
#[serde(rename_all = "camelCase")]
struct DiskProperties {
    radius_au: f64,
    surface_density_g_cm2: f64,
    temperature_k: f64,
    scale_height_au: f64,
    aspect_ratio: f64,
    midplane_density_g_cm3: f64,
    pressure_dyn_cm2: f64,
    sound_speed_cm_s: f64,
    keplerian_velocity_cm_s: f64,
    orbital_period_years: f64,
}
