//! Main simulation driver.
//!
//! Orchestrates all physical processes in the correct order each timestep.

use rand::Rng;
use units::{Length, Mass, SurfaceDensity, Time};

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

    // 3. Evolve particles (settling + coagulation)
    for bin in &mut state.particle_bins {
        // 3a. Settling: relax scale height toward equilibrium
        bin.relax_scale_height(&state.disk, dt);

        // 3b. Coagulation: evolve size distribution via collisions
        // Create kernel for this bin
        let kernel = crate::particles::Coagulation::new(bin, &state.disk);

        // Get current distribution parameters
        let surface_density = bin.surface_density().to_grams_per_cm2();
        let scale_height = bin.scale_height().to_cm();

        // Evolve via sticking (simplified - no fragmentation for now)
        let new_dist =
            kernel.evolve_sticking(bin.size_distribution(), surface_density, scale_height, dt);

        // Update the bin's size distribution
        *bin.size_distribution_mut() = new_dist;
    }

    // 4. Check for gravitational instability → planetesimal formation
    let stellar_mass = state.disk.stellar_mass();
    let mut new_planetesimals = Vec::new();

    for bin in &mut state.particle_bins {
        if let Some(event) = bin.attempt_planetesimal_formation(&state.disk, &mut state.rng) {
            // Convert planetesimal formation event to discrete bodies
            // Create a small number of representative bodies (max 5) to keep N manageable
            let num_formed = event.number_formed();
            let num_to_create = num_formed.clamp(1.0, 5.0) as usize;
            let mass_per_body = event.total_mass / (num_to_create as f64);

            for _ in 0..num_to_create {
                // Place in circular orbit at formation location
                // Add small random offsets to semi-major axis (±1% of width)
                let a_base = event.location.to_au();
                let da = 0.01 * state.rng.random::<f64>() - 0.005; // ±0.5%
                let a = a_base * (1.0 + da);

                // Small random eccentricity (0-0.05)
                let e = 0.05 * state.rng.random::<f64>();

                // Small random inclination (0-2 degrees)
                let i = (2.0 * std::f64::consts::PI / 180.0) * state.rng.random::<f64>();

                // Random mean anomaly (0-2π)
                let mean_anomaly = 2.0 * std::f64::consts::PI * state.rng.random::<f64>();

                // Create orbital elements
                let elements = OrbitalElements {
                    semi_major_axis: Length::from_au(a),
                    eccentricity: e,
                    inclination: i,
                    longitude_ascending_node: 2.0
                        * std::f64::consts::PI
                        * state.rng.random::<f64>(),
                    argument_of_periapsis: 2.0 * std::f64::consts::PI * state.rng.random::<f64>(),
                    mean_anomaly,
                };

                // Convert to Cartesian coordinates
                let (pos, vel) = orbital_elements_to_cartesian(&elements, stellar_mass);

                // Create planetesimal (rocky composition, no envelope)
                let planetesimal = crate::bodies::DiscreteBody::new(
                    mass_per_body,
                    Mass::zero(), // No envelope
                    pos,
                    vel,
                    stellar_mass,
                    planetary::composition::Composition::earth_like(), // Rocky
                );

                new_planetesimals.push(planetesimal);
            }
        }
    }

    // Add newly formed planetesimals to the simulation
    state.discrete_bodies.extend(new_planetesimals);

    // 5. Discrete body accretion from particle bins

    // Process each body-bin interaction
    for body in &mut state.discrete_bodies {
        let body_r = body.semi_major_axis;
        let feeding_zone_width = body.feeding_zone_width(stellar_mass);

        // Find particle bins that overlap with this body's feeding zone
        for bin in &mut state.particle_bins {
            let bin_r = bin.radial_center();
            let bin_width = bin.radial_width();

            // Check if bin overlaps with feeding zone (within 2.5 Hill radii)
            let distance = (body_r.to_au() - bin_r.to_au()).abs();
            let overlap_threshold = (feeding_zone_width.to_au() + bin_width.to_au()) / 2.0;

            if distance < overlap_threshold {
                // Calculate accretion rate from this bin
                let particle_sigma = bin.surface_density();
                let particle_v_disp = bin.velocity_dispersion();

                let dm_dt = body.accretion_rate_from_particles(
                    particle_sigma,
                    particle_v_disp,
                    stellar_mass,
                );

                // Calculate mass to accrete this timestep
                let mass_wanted = dm_dt.to_solar_masses_per_year() * dt.to_years();
                let mass_wanted_mass = Mass::from_solar_masses(mass_wanted);

                // Limit to available mass in bin (take at most 90% to avoid complete depletion)
                let mass_available = bin.total_mass();
                let mass_accreted = if mass_wanted_mass > mass_available {
                    mass_available * 0.9
                } else {
                    mass_wanted_mass
                };

                // Update body mass and radius
                body.core_mass = body.core_mass + mass_accreted;
                body.physical_radius = crate::bodies::DiscreteBody::estimate_radius(
                    body.total_mass(),
                    &body.composition,
                );

                // Remove mass from particle bin
                let sigma_to_remove = mass_accreted.to_grams() / bin.area();
                let new_sigma = SurfaceDensity::from_grams_per_cm2(
                    (particle_sigma.to_grams_per_cm2() - sigma_to_remove).max(0.0),
                );
                bin.set_surface_density(new_sigma);
            }
        }
    }

    // 6. Envelope evolution (hydrostatic + runaway accretion)
    const OPACITY: f64 = 0.01; // Rosseland mean opacity in cm²/g
    for body in &mut state.discrete_bodies {
        // Calculate core accretion rate from disk viscosity
        let core_accretion_rate = supply_limited_accretion_rate(&state.disk, body.semi_major_axis);

        // Evolve envelope (captures gas, transitions to runaway, etc.)
        body.evolve_envelope(&state.disk, core_accretion_rate, dt, OPACITY);
    }

    // 7. Migration (Type I for embedded planets)
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

    // 8. Collision detection and mergers
    // Check all pairs for collisions (physical contact)
    let mut merged_indices = Vec::new();
    let mut new_bodies = Vec::new();

    for i in 0..state.discrete_bodies.len() {
        if merged_indices.contains(&i) {
            continue; // Already merged
        }

        for j in (i + 1)..state.discrete_bodies.len() {
            if merged_indices.contains(&j) {
                continue; // Already merged
            }

            let body_i = &state.discrete_bodies[i];
            let body_j = &state.discrete_bodies[j];

            // Calculate distance between bodies
            let dx = body_i.position.x - body_j.position.x;
            let dy = body_i.position.y - body_j.position.y;
            let distance_au = (dx * dx + dy * dy).sqrt();
            let distance = Length::from_au(distance_au);

            // Check for physical collision (radii overlap)
            let collision_distance = body_i.physical_radius + body_j.physical_radius;

            if distance < collision_distance {
                // Collision! Merge the two bodies
                let mi = body_i.total_mass();
                let mj = body_j.total_mass();
                let total_mass = mi + mj;

                // Conserve momentum to get new velocity
                let px = mi.to_solar_masses() * body_i.velocity.x
                    + mj.to_solar_masses() * body_j.velocity.x;
                let py = mi.to_solar_masses() * body_i.velocity.y
                    + mj.to_solar_masses() * body_j.velocity.y;

                let new_velocity = nalgebra::Vector2::new(
                    px / total_mass.to_solar_masses(),
                    py / total_mass.to_solar_masses(),
                );

                // Center of mass position
                let new_position = nalgebra::Point2::new(
                    (mi.to_solar_masses() * body_i.position.x
                        + mj.to_solar_masses() * body_j.position.x)
                        / total_mass.to_solar_masses(),
                    (mi.to_solar_masses() * body_i.position.y
                        + mj.to_solar_masses() * body_j.position.y)
                        / total_mass.to_solar_masses(),
                );

                // Add masses
                let new_core_mass = body_i.core_mass + body_j.core_mass;
                let new_envelope_mass = body_i.envelope_mass + body_j.envelope_mass;

                // Mass-weighted composition
                let mi_solar = mi.to_solar_masses();
                let mj_solar = mj.to_solar_masses();
                let total_solar = total_mass.to_solar_masses();

                let new_composition = planetary::composition::Composition::new(
                    (body_i.composition.iron * mi_solar + body_j.composition.iron * mj_solar)
                        / total_solar,
                    (body_i.composition.silicate * mi_solar
                        + body_j.composition.silicate * mj_solar)
                        / total_solar,
                    (body_i.composition.water * mi_solar + body_j.composition.water * mj_solar)
                        / total_solar,
                    (body_i.composition.h_he_gas * mi_solar
                        + body_j.composition.h_he_gas * mj_solar)
                        / total_solar,
                );

                // Create merged body
                let merged = crate::bodies::DiscreteBody::new(
                    new_core_mass,
                    new_envelope_mass,
                    new_position,
                    new_velocity,
                    stellar_mass,
                    new_composition,
                );

                new_bodies.push(merged);
                merged_indices.push(i);
                merged_indices.push(j);

                break; // Body i is merged, move to next i
            }
        }
    }

    // Remove merged bodies (in reverse order to preserve indices)
    merged_indices.sort_unstable();
    merged_indices.reverse();
    for &idx in &merged_indices {
        state.discrete_bodies.remove(idx);
    }

    // Add new merged bodies
    state.discrete_bodies.extend(new_bodies);

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
