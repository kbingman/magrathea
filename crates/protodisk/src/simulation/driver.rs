//! Main simulation driver.
//!
//! Orchestrates all physical processes in the correct order each timestep.

use units::{Length, Time};

use super::{SimulationState, calculate_timestep};
use crate::bodies::{
    OrbitalElements, migration_distance, orbital_elements_to_cartesian,
    supply_limited_accretion_rate, type_i_migration_timescale,
};
use crate::disk::DiskModel;

/// Execute one simulation timestep.
///
/// Updates all components of the simulation state in the correct order:
/// 1. Calculate adaptive timestep
/// 2. Evolve gas disk (viscous + photoevaporation)
/// 3. Evolve particles (drift, settling, coagulation)
/// 4. Check for gravitational instability → planetesimal formation
/// 5. Discrete body accretion from gas and particles
/// 6. Envelope evolution (including runaway accretion)
/// 7. Planetary migration
/// 8. Collision detection and mergers
///
/// # Arguments
/// * `state` - Current simulation state (will be modified in place)
///
/// # Returns
/// Timestep that was taken
pub fn step(state: &mut SimulationState) -> Time {
    // 1. Calculate adaptive timestep
    let dt = calculate_timestep(state);
    let dt_seconds = dt.to_seconds();

    // 2. Evolve disk (viscous spreading + photoevaporation)
    state.disk.evolve_viscous(dt_seconds);
    state.disk.apply_photoevaporation(dt_seconds);

    // 3. Evolve particles
    // TODO: Implement particle evolution

    // 4. Check for gravitational instability → planetesimal formation
    for bin in &mut state.particle_bins {
        if let Some(_event) = bin.attempt_planetesimal_formation(&state.disk, &mut state.rng) {
            // TODO: Convert planetesimals to discrete bodies
            // For now, planetesimals are removed from the bin but not tracked
            // This depletes the particle population as expected
        }
    }

    // 5. Discrete body accretion
    // TODO: Implement accretion

    // 6. Envelope evolution (hydrostatic + runaway accretion)
    const OPACITY: f64 = 0.01; // Rosseland mean opacity in cm²/g
    for body in &mut state.discrete_bodies {
        // Calculate core accretion rate from disk viscosity
        let core_accretion_rate = supply_limited_accretion_rate(&state.disk, body.semi_major_axis);

        // Evolve envelope (captures gas, transitions to runaway, etc.)
        body.evolve_envelope(&state.disk, core_accretion_rate, dt, OPACITY);
    }

    // 7. Migration (Type I for embedded planets)
    let stellar_mass = state.disk.stellar_mass();
    for body in &mut state.discrete_bodies {
        let r = body.semi_major_axis;
        let sigma = state.disk.surface_density(r);
        let h_over_r = state.disk.aspect_ratio(r);

        // Calculate migration timescale and distance
        let tau_mig =
            type_i_migration_timescale(body.total_mass(), stellar_mass, r, sigma, h_over_r);
        let da = migration_distance(r, tau_mig, dt);

        // Update semi-major axis (da is negative for inward migration)
        let new_a = Length::from_au((r.to_au() + da).max(state.disk.inner_radius().to_au()));

        // Create orbital elements with new semi-major axis
        let new_elements = OrbitalElements {
            semi_major_axis: new_a,
            eccentricity: body.eccentricity,
            inclination: body.inclination,
            longitude_ascending_node: 0.0,
            argument_of_periapsis: 0.0,
            mean_anomaly: 0.0, // Keep at same orbital phase (approximation)
        };

        // Convert to Cartesian and update body
        let (new_pos, new_vel) = orbital_elements_to_cartesian(&new_elements, stellar_mass);
        body.position = new_pos;
        body.velocity = new_vel;
        body.semi_major_axis = new_a;
    }

    // 8. Collision detection
    // TODO: Implement collision detection

    // Update simulation time
    state.time = state.time + dt;

    dt
}

/// Run simulation to completion.
///
/// Executes timesteps until either:
/// - Simulation time exceeds `max_time`
/// - Gas disk is fully dispersed
///
/// # Arguments
/// * `state` - Initial simulation state
/// * `max_time` - Maximum simulation time
///
/// # Returns
/// Final simulation state
///
/// # Example
/// ```rust
/// use protodisk::*;
/// use units::Time;
///
/// let star = solar_analog();
/// let gas_disk = GasDisk::for_star(&star);
/// let grid_disk = GridDisk::from_gas_disk(&gas_disk, 50);
/// let state = SimulationState::new(grid_disk, vec![], vec![]);
///
/// let max_time = Time::from_years(1_000_000.0);
/// let final_state = run_simulation(state, max_time);
///
/// println!("Final time: {} years", final_state.time.to_years());
/// println!("Bodies formed: {}", final_state.body_count());
/// ```
pub fn run_simulation(mut state: SimulationState, max_time: Time) -> SimulationState {
    while !state.is_finished(max_time) {
        step(&mut state);
    }
    state
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::disk::GasDisk;
    use crate::solar_analog;

    #[test]
    fn step_advances_time() {
        let star = solar_analog();
        let gas_disk = GasDisk::for_star(&star);
        let grid_disk = crate::disk::GridDisk::from_gas_disk(&gas_disk, 50);

        let mut state = SimulationState::new(grid_disk, vec![], vec![]);
        let initial_time = state.time;

        let dt = step(&mut state);

        assert!(dt.to_years() > 0.0);
        assert!(state.time.to_years() > initial_time.to_years());
    }

    #[test]
    fn run_simulation_completes() {
        let star = solar_analog();
        let gas_disk = GasDisk::for_star(&star);
        let grid_disk = crate::disk::GridDisk::from_gas_disk(&gas_disk, 50);

        let state = SimulationState::new(grid_disk, vec![], vec![]);
        let max_time = Time::from_years(100_000.0);

        let final_state = run_simulation(state, max_time);

        assert!(final_state.time.to_years() > 0.0);
        assert!(final_state.is_finished(max_time));
    }

    #[test]
    fn simulation_respects_max_time() {
        let star = solar_analog();
        let gas_disk = GasDisk::for_star(&star);
        let grid_disk = crate::disk::GridDisk::from_gas_disk(&gas_disk, 50);

        let state = SimulationState::new(grid_disk, vec![], vec![]);
        let max_time = Time::from_years(10_000.0);

        let final_state = run_simulation(state, max_time);

        assert!(final_state.time <= max_time);
    }

    #[test]
    fn disk_evolves_during_step() {
        use crate::disk::DiskMass;

        let star = solar_analog();
        let gas_disk = GasDisk::for_star(&star);
        let grid_disk = crate::disk::GridDisk::from_gas_disk(&gas_disk, 50);

        let mut state = SimulationState::new(grid_disk, vec![], vec![]);
        let initial_mass = state.disk.total_mass();

        // Take one timestep
        step(&mut state);

        let final_mass = state.disk.total_mass();

        // Disk mass should decrease due to photoevaporation and accretion
        assert!(
            final_mass < initial_mass,
            "Disk mass should decrease: initial={} final={}",
            initial_mass.to_solar_masses(),
            final_mass.to_solar_masses()
        );
    }

    #[test]
    fn bodies_migrate_inward() {
        use crate::bodies::DiscreteBody;
        use nalgebra::{Point2, Vector2};
        use planetary::composition::Composition;
        use units::Mass;

        const G: f64 = 39.478417; // AU³ M☉⁻¹ year⁻²

        let star = solar_analog();
        let stellar_mass = star.mass;
        let gas_disk = GasDisk::for_star(&star);
        let grid_disk = crate::disk::GridDisk::from_gas_disk(&gas_disk, 50);

        // Create a super-Earth at 5 AU in a circular orbit
        let a = 5.0; // AU
        let v_circ = (G * stellar_mass.to_solar_masses() / a).sqrt();

        let body = DiscreteBody::new(
            Mass::from_earth_masses(10.0), // core
            Mass::zero(),                  // no envelope
            Point2::new(a, 0.0),           // position at apoapse
            Vector2::new(0.0, v_circ),     // circular velocity
            stellar_mass,
            Composition::earth_like(),
        );

        let mut state = SimulationState::new(grid_disk, vec![], vec![body]);
        let initial_a = state.discrete_bodies[0].semi_major_axis;

        // Take multiple timesteps to see migration
        for _ in 0..100 {
            step(&mut state);
        }

        let final_a = state.discrete_bodies[0].semi_major_axis;

        // Body should migrate inward (Type I migration)
        assert!(
            final_a < initial_a,
            "Body should migrate inward: initial={} AU, final={} AU",
            initial_a.to_au(),
            final_a.to_au()
        );
    }

    #[test]
    fn envelope_evolution_works() {
        use crate::bodies::{DiscreteBody, EnvelopeState};
        use nalgebra::{Point2, Vector2};
        use planetary::composition::Composition;
        use units::Mass;

        const G: f64 = 39.478417; // AU³ M☉⁻¹ year⁻²

        let star = solar_analog();
        let stellar_mass = star.mass;
        let gas_disk = GasDisk::for_star(&star);
        let grid_disk = crate::disk::GridDisk::from_gas_disk(&gas_disk, 50);

        // Create a 1 M⊕ embryo at 5 AU (sufficient for envelope capture)
        let a = 5.0; // AU
        let v_circ = (G * stellar_mass.to_solar_masses() / a).sqrt();

        let body = DiscreteBody::new(
            Mass::from_earth_masses(1.0), // core
            Mass::zero(),                 // no initial envelope
            Point2::new(a, 0.0),
            Vector2::new(0.0, v_circ),
            stellar_mass,
            Composition::earth_like(),
        );

        let mut state = SimulationState::new(grid_disk, vec![], vec![body]);

        // Initially should have no envelope
        assert!(matches!(
            state.discrete_bodies[0].envelope_state,
            EnvelopeState::None
        ));

        // Run simulation for some time
        for _ in 0..10 {
            step(&mut state);
        }

        // Should have captured an envelope (hydrostatic or runaway)
        assert!(
            !matches!(state.discrete_bodies[0].envelope_state, EnvelopeState::None),
            "Body should have captured envelope"
        );

        // Envelope mass should be positive
        assert!(
            state.discrete_bodies[0].envelope_mass.to_earth_masses() > 0.0,
            "Envelope mass should be positive"
        );
    }

    #[test]
    fn mass_approximately_conserved() {
        let star = solar_analog();
        let gas_disk = GasDisk::for_star(&star);
        let grid_disk = crate::disk::GridDisk::from_gas_disk(&gas_disk, 50);

        let mut state = SimulationState::new(grid_disk, vec![], vec![]);
        let initial_mass = state.total_mass();

        println!("Initial total mass: {} M☉", initial_mass.to_solar_masses());

        // Run for 10 timesteps
        for i in 0..10 {
            step(&mut state);
            let current_mass = state.total_mass();
            let mass_lost = initial_mass - current_mass;

            println!(
                "Step {}: time = {:.1} yr, mass = {} M☉, lost = {} M☉",
                i + 1,
                state.time.to_years(),
                current_mass.to_solar_masses(),
                mass_lost.to_solar_masses()
            );
        }

        let final_mass = state.total_mass();
        let mass_lost = initial_mass - final_mass;

        // Mass should only decrease due to photoevaporation (and accretion onto star)
        // Loss should be relatively small over short timescales
        assert!(
            final_mass <= initial_mass,
            "Mass cannot increase: initial={}, final={}",
            initial_mass.to_solar_masses(),
            final_mass.to_solar_masses()
        );

        // Check that mass loss is reasonable (< 1% over this short time)
        let fractional_loss = mass_lost.to_solar_masses() / initial_mass.to_solar_masses();
        assert!(
            fractional_loss < 0.01,
            "Mass loss too large: {:.4}% in {} years",
            fractional_loss * 100.0,
            state.time.to_years()
        );

        println!(
            "✓ Mass conservation validated: {:.6}% lost over {:.1} years",
            fractional_loss * 100.0,
            state.time.to_years()
        );
    }
}
