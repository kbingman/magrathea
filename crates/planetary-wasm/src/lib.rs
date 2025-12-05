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

/// Convert a Rust type to JsValue via serde.
fn to_js<T: serde::Serialize>(value: &T) -> Result<JsValue, JsError> {
    serde_wasm_bindgen::to_value(value).map_err(|e| JsError::new(&e.to_string()))
}

/// Convert a JsValue to a Rust type via serde.
fn from_js<T: serde::de::DeserializeOwned>(value: JsValue) -> Result<T, JsError> {
    serde_wasm_bindgen::from_value(value).map_err(|e| JsError::new(&e.to_string()))
}
