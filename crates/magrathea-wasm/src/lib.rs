//! WASM bindings for Magrathea stellar and disk generation.
//!
//! This crate provides JavaScript/TypeScript bindings for the core Magrathea
//! library using `wasm-bindgen` and `serde-wasm-bindgen` for seamless type
//! conversion.

use wasm_bindgen::prelude::*;

mod disk;
mod stellar;

fn to_js<T: serde::Serialize>(value: &T) -> Result<JsValue, JsError> {
    serde_wasm_bindgen::to_value(value).map_err(|e| JsError::new(&e.to_string()))
}

fn from_js<T: serde::de::DeserializeOwned>(value: JsValue) -> Result<T, JsError> {
    serde_wasm_bindgen::from_value(value).map_err(|e| JsError::new(&e.to_string()))
}
