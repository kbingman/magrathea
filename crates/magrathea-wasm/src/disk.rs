//! WASM bindings for protoplanetary disk models.

use std::cell::RefCell;
use std::collections::HashMap;

use wasm_bindgen::prelude::*;

use stellar_forge::{DiskMass, DiskModel, GasDisk, GridDisk, MainSequenceStar};

use crate::{from_js, to_js};

// Thread-local storage for grid disks (WASM is single-threaded)
thread_local! {
    static GRID_DISKS: RefCell<HashMap<u32, GridDisk>> = RefCell::new(HashMap::new());
    static NEXT_DISK_ID: RefCell<u32> = const { RefCell::new(0) };
}

// =============================================================================
// GasDisk bindings (power-law model)
// =============================================================================

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

// =============================================================================
// GridDisk bindings (mutable, viscous evolution)
// =============================================================================

/// Create a GridDisk from a star's gas disk.
///
/// Returns a disk ID that can be used with other grid_disk_* functions.
///
/// # Arguments
/// * `star` - A MainSequenceStar object
/// * `n_radii` - Number of radial grid points (100-500 recommended)
#[wasm_bindgen]
pub fn grid_disk_create(star: JsValue, n_radii: usize) -> Result<u32, JsError> {
    let star: MainSequenceStar = from_js(star)?;
    let gas_disk = GasDisk::for_star(&star);
    let grid_disk = GridDisk::from_gas_disk(&gas_disk, n_radii);

    let id = NEXT_DISK_ID.with(|next_id| {
        let mut id = next_id.borrow_mut();
        let current = *id;
        *id += 1;
        current
    });

    GRID_DISKS.with(|disks| {
        disks.borrow_mut().insert(id, grid_disk);
    });

    Ok(id)
}

/// Create a GridDisk from a GasDisk.
///
/// # Arguments
/// * `gas_disk` - A GasDisk object
/// * `n_radii` - Number of radial grid points
#[wasm_bindgen]
pub fn grid_disk_from_gas_disk(gas_disk: JsValue, n_radii: usize) -> Result<u32, JsError> {
    let gas_disk: GasDisk = from_js(gas_disk)?;
    let grid_disk = GridDisk::from_gas_disk(&gas_disk, n_radii);

    let id = NEXT_DISK_ID.with(|next_id| {
        let mut id = next_id.borrow_mut();
        let current = *id;
        *id += 1;
        current
    });

    GRID_DISKS.with(|disks| {
        disks.borrow_mut().insert(id, grid_disk);
    });

    Ok(id)
}

/// Evolve a GridDisk for the specified number of timesteps.
///
/// Uses the maximum stable timestep for each step.
///
/// # Arguments
/// * `disk_id` - The disk ID returned by grid_disk_create
/// * `n_steps` - Number of evolution steps
///
/// # Returns
/// The total elapsed time in years.
#[wasm_bindgen]
pub fn grid_disk_evolve(disk_id: u32, n_steps: usize) -> Result<f64, JsError> {
    GRID_DISKS.with(|disks| {
        let mut disks = disks.borrow_mut();
        let disk = disks
            .get_mut(&disk_id)
            .ok_or_else(|| JsError::new(&format!("Disk {} not found", disk_id)))?;

        let dt = disk.max_timestep();
        for _ in 0..n_steps {
            disk.evolve_viscous(dt);
        }

        Ok(units::Time::from_seconds(n_steps as f64 * dt).to_years())
    })
}

/// Get the current density profile of a GridDisk.
///
/// # Returns
/// Object with radii (AU), sigma (g/cm²), innerAu, outerAu, totalMassSolar
#[wasm_bindgen]
pub fn grid_disk_profile(disk_id: u32) -> Result<JsValue, JsError> {
    GRID_DISKS.with(|disks| {
        let disks = disks.borrow();
        let disk = disks
            .get(&disk_id)
            .ok_or_else(|| JsError::new(&format!("Disk {} not found", disk_id)))?;

        let radii_au: Vec<f64> = disk
            .radii()
            .iter()
            .map(|&r_cm| units::Length::from_cm(r_cm).to_au())
            .collect();

        to_js(&DiskProfile {
            radii_au,
            sigma_g_cm2: disk.sigma().to_vec(),
            inner_au: disk.inner_radius().to_au(),
            outer_au: disk.outer_radius().to_au(),
            total_mass_solar: disk.total_mass().to_solar_masses(),
        })
    })
}

/// Delete a GridDisk to free memory.
#[wasm_bindgen]
pub fn grid_disk_delete(disk_id: u32) {
    GRID_DISKS.with(|disks| {
        disks.borrow_mut().remove(&disk_id);
    });
}

/// Profile data for a GridDisk.
#[derive(serde::Serialize)]
#[serde(rename_all = "camelCase")]
struct DiskProfile {
    radii_au: Vec<f64>,
    sigma_g_cm2: Vec<f64>,
    inner_au: f64,
    outer_au: f64,
    total_mass_solar: f64,
}
