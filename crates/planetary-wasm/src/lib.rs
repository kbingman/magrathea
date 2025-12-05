//! WASM bindings for statistical planetary system generation.
//!
//! This crate provides JavaScript/TypeScript bindings for generating and
//! classifying planetary systems using occurrence-rate based statistical sampling.
//!
//! ## Quick Start (JavaScript)
//!
//! ```javascript
//! import init, {
//!     generate_system_for_solar_analog,
//!     generate_system_random,
//!     generate_system_from_seed,
//! } from 'planetary-wasm';
//!
//! await init();
//!
//! // Generate a system around a Sun-like star
//! const system = generate_system_for_solar_analog();
//! console.log(`System ${system.metadata.catalogName} has ${system.planets.length} planets`);
//!
//! // Reproducible generation from a seed string
//! const system2 = generate_system_from_seed("my-favorite-system");
//! ```

use wasm_bindgen::prelude::*;

mod stellar;
mod system;

// Type aliases for unit types (serialized as numbers via serde(transparent))
#[wasm_bindgen(typescript_custom_section)]
const TS_UNIT_TYPES: &'static str = r#"
/** Mass in solar masses (M☉) */
export type Mass = number;
/** Length in AU (astronomical units) */
export type Length = number;
/** Temperature in Kelvin */
export type Temperature = number;
/** Time in Myr (million years) */
export type Time = number;
/** Stellar radius in solar radii (R☉) */
export type StellarRadius = number;
/** UUID string */
export type Uuid = string;
"#;
