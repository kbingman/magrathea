//! Simulation state container.
//!
//! Holds all the evolving components of the planet formation simulation.

use units::{Mass, Time};

use crate::bodies::DiscreteBody;
use crate::disk::{DiskMass, GridDisk};
use crate::particles::ParticleBin;

/// Complete state of the planet formation simulation.
///
/// This structure contains all the evolving populations and tracks
/// the current simulation time.
#[derive(Debug, Clone)]
pub struct SimulationState {
    /// Current simulation time
    pub time: Time,

    /// Gas disk (evolves via viscosity and photoevaporation)
    pub disk: GridDisk,

    /// Particle populations in radial bins
    ///
    /// Each bin contains a size distribution of particles from dust to
    /// km-sized planetesimals. Particles evolve via drift, settling,
    /// and coagulation.
    pub particle_bins: Vec<ParticleBin>,

    /// Individually tracked bodies
    ///
    /// Bodies that have grown large enough to be tracked individually
    /// (typically > 0.01 M⊕). These interact gravitationally and accrete
    /// from particle bins.
    pub discrete_bodies: Vec<DiscreteBody>,
}

impl SimulationState {
    /// Create a new simulation state.
    ///
    /// # Arguments
    /// * `disk` - Initial gas disk configuration
    /// * `particle_bins` - Initial particle populations
    /// * `discrete_bodies` - Initial discrete bodies (often empty)
    ///
    /// # Returns
    /// New simulation state at t=0
    pub fn new(
        disk: GridDisk,
        particle_bins: Vec<ParticleBin>,
        discrete_bodies: Vec<DiscreteBody>,
    ) -> Self {
        Self {
            time: Time::zero(),
            disk,
            particle_bins,
            discrete_bodies,
        }
    }

    /// Calculate total mass in the system.
    ///
    /// Sums mass in:
    /// - Gas disk
    /// - Particle bins
    /// - Discrete bodies
    ///
    /// # Returns
    /// Total mass in solar masses
    ///
    /// # Conservation Check
    ///
    /// This should remain constant (within numerical error) throughout
    /// the simulation, except for photoevaporative mass loss.
    pub fn total_mass(&self) -> Mass {
        let disk_mass = self.disk.total_mass();

        let particle_mass = self
            .particle_bins
            .iter()
            .map(|bin| bin.total_mass())
            .fold(Mass::zero(), |acc, m| acc + m);

        let body_mass = self
            .discrete_bodies
            .iter()
            .map(|body| body.total_mass())
            .fold(Mass::zero(), |acc, m| acc + m);

        disk_mass + particle_mass + body_mass
    }

    /// Count total number of discrete bodies.
    pub fn body_count(&self) -> usize {
        self.discrete_bodies.len()
    }

    /// Count total number of particle bins.
    pub fn particle_bin_count(&self) -> usize {
        self.particle_bins.len()
    }

    /// Get the most massive discrete body.
    ///
    /// # Returns
    /// Reference to the most massive body, or None if no bodies exist
    pub fn most_massive_body(&self) -> Option<&DiscreteBody> {
        self.discrete_bodies
            .iter()
            .max_by(|a, b| a.total_mass().partial_cmp(&b.total_mass()).unwrap())
    }

    /// Check if simulation has ended.
    ///
    /// Simulation ends when:
    /// - Time exceeds maximum (typically 10 Myr)
    /// - Gas disk is fully dispersed
    ///
    /// # Arguments
    /// * `max_time` - Maximum simulation time
    ///
    /// # Returns
    /// true if simulation should stop
    pub fn is_finished(&self, max_time: Time) -> bool {
        // Time limit reached
        if self.time >= max_time {
            return true;
        }

        // Disk fully dispersed (less than 0.1% of Jupiter mass in gas)
        let min_mass = Mass::from_jupiter_masses(0.001);
        if self.disk.total_mass() < min_mass {
            return true;
        }

        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::disk::GasDisk;
    use crate::solar_analog;

    #[test]
    fn new_state_starts_at_t_zero() {
        let star = solar_analog();
        let gas_disk = GasDisk::for_star(&star);
        let grid_disk = GridDisk::from_gas_disk(&gas_disk, 50);

        let state = SimulationState::new(grid_disk, vec![], vec![]);

        assert_eq!(state.time.to_years(), 0.0);
    }

    #[test]
    fn total_mass_includes_all_components() {
        let star = solar_analog();
        let gas_disk = GasDisk::for_star(&star);
        let grid_disk = GridDisk::from_gas_disk(&gas_disk, 50);

        let disk_mass = grid_disk.total_mass();

        let state = SimulationState::new(grid_disk, vec![], vec![]);

        // With no particles or bodies, total mass should equal disk mass
        let total = state.total_mass();
        let ratio = total.to_solar_masses() / disk_mass.to_solar_masses();

        assert!((ratio - 1.0).abs() < 1e-10, "Mass mismatch: {}", ratio);
    }

    #[test]
    fn empty_state_has_zero_bodies() {
        let star = solar_analog();
        let gas_disk = GasDisk::for_star(&star);
        let grid_disk = GridDisk::from_gas_disk(&gas_disk, 50);

        let state = SimulationState::new(grid_disk, vec![], vec![]);

        assert_eq!(state.body_count(), 0);
        assert!(state.most_massive_body().is_none());
    }

    #[test]
    fn is_finished_when_time_exceeds_max() {
        let star = solar_analog();
        let gas_disk = GasDisk::for_star(&star);
        let grid_disk = GridDisk::from_gas_disk(&gas_disk, 50);

        let mut state = SimulationState::new(grid_disk, vec![], vec![]);
        let max_time = Time::from_years(1_000_000.0);

        assert!(!state.is_finished(max_time));

        state.time = Time::from_years(2_000_000.0);
        assert!(state.is_finished(max_time));
    }

    // TODO: Fix this test - needs GridDisk mutation API
    // #[test]
    // fn is_finished_when_disk_dispersed() {
    //     let star = solar_analog();
    //     let gas_disk = GasDisk::for_star(&star);
    //     let mut grid_disk = GridDisk::from_gas_disk(&gas_disk, 50);
    //
    //     // Artificially disperse disk by setting all surface densities to zero
    //     for i in 0..grid_disk.n_radii() {
    //         grid_disk.set_sigma(i, 0.0);
    //     }
    //
    //     let state = SimulationState::new(grid_disk, vec![], vec![]);
    //     let max_time = Time::from_years(10_000_000.0);
    //
    //     assert!(state.is_finished(max_time), "Should be finished when disk dispersed");
    // }
}
