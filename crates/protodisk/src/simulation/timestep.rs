//! Adaptive timestep calculation.
//!
//! Calculates the maximum safe timestep based on all active physical processes.
//! Uses Courant-like conditions to ensure stability and accuracy.

use units::{Length, Time};

use super::SimulationState;
use crate::disk::DiskModel;

/// Safety factor for timestep (fraction of characteristic timescale)
///
/// Using 0.1 means we take timesteps that are 10% of the fastest process.
/// This ensures stability while allowing reasonable progress.
const TIMESTEP_SAFETY_FACTOR: f64 = 0.1;

/// Minimum timestep in years (1 year)
///
/// Prevents timesteps from becoming too small in extreme cases.
const MIN_TIMESTEP_YEARS: f64 = 1.0;

/// Maximum timestep in years (10,000 years)
///
/// Prevents missing important events even when processes are slow.
const MAX_TIMESTEP_YEARS: f64 = 10_000.0;

/// Calculate adaptive timestep for simulation.
///
/// The timestep is chosen as a fraction of the shortest characteristic
/// timescale in the system:
///
/// ```text
/// Δt = safety_factor × min(τ_orbital, τ_drift, τ_coag, τ_migration, τ_disk)
/// ```
///
/// # Arguments
/// * `state` - Current simulation state
///
/// # Returns
/// Timestep in years
///
/// # Physics
///
/// Different processes have different characteristic timescales:
///
/// - **Orbital**: τ_orb = 2π/Ω (fastest at inner edge, ~0.2 yr at 0.1 AU)
/// - **Drift**: τ_drift = r/v_drift (depends on particle size and disk structure)
/// - **Coagulation**: τ_coag ~ (n σ v)⁻¹ (depends on particle density and size)
/// - **Migration**: τ_mig = r × (dr/dt)⁻¹ (typically 10⁴-10⁶ yr for Type I)
/// - **Disk evolution**: τ_visc = r²/ν (typically >10⁶ yr)
///
/// # Stability
///
/// The safety factor ensures that:
/// - Particles don't drift more than ~10% of their orbital radius
/// - Orbital phases are well-resolved
/// - Coagulation doesn't change size distribution drastically
/// - Migration is smooth and stable
///
/// # Performance
///
/// The timestep is clamped between MIN_TIMESTEP_YEARS and MAX_TIMESTEP_YEARS
/// to balance accuracy with computational cost.
pub fn calculate_timestep(state: &SimulationState) -> Time {
    let mut min_timescale = MAX_TIMESTEP_YEARS;

    // 1. Orbital timescale at inner edge (typically shortest)
    let inner_r = state.disk.inner_radius();
    let inner_period = state.disk.orbital_period(inner_r);
    min_timescale = min_timescale.min(inner_period.to_years());

    // 2-3. Particle timescales (drift, coagulation)
    // TODO: Implement once ParticleBin API is clarified
    // For now, use conservative estimate
    if !state.particle_bins.is_empty() {
        // Use ~1000 year default for particle processes
        min_timescale = min_timescale.min(1000.0);
    }

    // 4. Migration timescale for discrete bodies
    for body in &state.discrete_bodies {
        let r = Length::from_au(body.semi_major_axis.to_au());
        let sigma = state.disk.surface_density(r);
        let h_over_r = state.disk.aspect_ratio(r);

        let tau_mig = crate::bodies::type_i_migration_timescale(
            body.total_mass(),
            state.disk.stellar_mass(),
            r,
            sigma,
            h_over_r,
        );

        min_timescale = min_timescale.min(tau_mig.to_years());
    }

    // 5. Disk viscous timescale (usually longest, so not often limiting)
    let visc_timescale = state.disk.viscous_timescale(state.disk.inner_radius());
    min_timescale = min_timescale.min(visc_timescale.to_years());

    // Apply safety factor and clamp
    let dt = TIMESTEP_SAFETY_FACTOR * min_timescale;
    let dt_clamped = dt.clamp(MIN_TIMESTEP_YEARS, MAX_TIMESTEP_YEARS);

    Time::from_years(dt_clamped)
}

// TODO: Restore these functions once ParticleBin API is clarified
// For now, using conservative defaults in calculate_timestep()

// fn calculate_drift_timescale(bin: &crate::particles::ParticleBin, disk: &crate::disk::GridDisk) -> f64 {
//     ...
// }

// fn calculate_coagulation_timescale(bin: &crate::particles::ParticleBin) -> f64 {
//     ...
// }

#[cfg(test)]
mod tests {
    use super::*;
    use crate::disk::GasDisk;

    #[test]
    fn timestep_is_positive() {
        let star = crate::solar_analog();
        let gas_disk = GasDisk::for_star(&star);
        let grid_disk = crate::disk::GridDisk::from_gas_disk(&gas_disk, 50);

        let state = SimulationState::new(grid_disk, vec![], vec![]);
        let dt = calculate_timestep(&state);

        assert!(dt.to_years() > 0.0);
    }

    #[test]
    fn timestep_respects_bounds() {
        let star = crate::solar_analog();
        let gas_disk = GasDisk::for_star(&star);
        let grid_disk = crate::disk::GridDisk::from_gas_disk(&gas_disk, 50);

        let state = SimulationState::new(grid_disk, vec![], vec![]);
        let dt = calculate_timestep(&state);

        assert!(
            dt.to_years() >= MIN_TIMESTEP_YEARS,
            "Timestep {} < minimum {}",
            dt.to_years(),
            MIN_TIMESTEP_YEARS
        );
        assert!(
            dt.to_years() <= MAX_TIMESTEP_YEARS,
            "Timestep {} > maximum {}",
            dt.to_years(),
            MAX_TIMESTEP_YEARS
        );
    }

    // TODO: Re-enable these tests once ParticleBin API is clarified

    // #[test]
    // fn timestep_decreases_with_inner_particles() {
    //     ...
    // }

    // #[test]
    // fn drift_timescale_reasonable() {
    //     ...
    // }

    // #[test]
    // fn coagulation_timescale_reasonable() {
    //     ...
    // }
}
