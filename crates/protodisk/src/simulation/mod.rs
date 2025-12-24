//! Emergent planet formation simulation driver.
//!
//! This module orchestrates all the physics modules (disk evolution, particle
//! dynamics, planetesimal formation, discrete body growth, envelope accretion,
//! and migration) into a coherent time-stepping simulation.
//!
//! # Architecture
//!
//! The simulation maintains three populations:
//! 1. **Gas disk** - evolves via viscous spreading and photoevaporation
//! 2. **Particle bins** - statistical representation of small bodies (dust to km-sized)
//! 3. **Discrete bodies** - individually tracked protoplanets and planets
//!
//! Particles transition to discrete bodies when they reach the tracking threshold
//! (typically ~0.01-0.1 M⊕) or form via gravitational instability.
//!
//! # Time-Stepping
//!
//! The simulation uses adaptive timesteps based on the fastest physical process:
//! - Orbital period (shortest at inner edge)
//! - Particle drift timescale
//! - Coagulation timescale
//! - Migration timescale
//! - Disk viscous timescale
//!
//! # Update Sequence
//!
//! Each timestep proceeds in this order:
//! 1. Disk evolution (viscous spreading, photoevaporation)
//! 2. Particle evolution (drift, settling, coagulation)
//! 3. Gravitational instability check → planetesimal formation
//! 4. Discrete body accretion from particle bins
//! 5. Gas envelope evolution
//! 6. Planetary migration
//! 7. Collision detection and mergers
//!
//! # Validation
//!
//! Conservation laws are checked each timestep:
//! - Total mass (disk + particles + bodies)
//! - Angular momentum
//!
//! # References
//! - Ida & Lin (2004) - "Toward a Deterministic Model of Planetary Formation"
//! - Chambers (2006) - "A semi-analytic model for oligarchic growth"
//! - Bitsch et al. (2015) - "The growth of planets by pebble accretion"

mod driver;
mod state;
mod timestep;

#[cfg(test)]
mod integration_test;

pub use driver::{run_simulation, step};
pub use state::SimulationState;
pub use timestep::calculate_timestep;
