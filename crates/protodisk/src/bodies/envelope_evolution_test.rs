//! Tests for gas envelope evolution and runaway accretion.

use super::*;
use crate::disk::GasDisk;
use nalgebra::{Point2, Vector2};
use planetary::composition::Composition;
use units::{Mass, MassRate, Time};

#[test]
fn envelope_starts_none_for_small_core() {
    let body = DiscreteBody::new(
        Mass::from_earth_masses(0.05), // Too small for envelope
        Mass::zero(),
        Point2::new(5.0, 0.0),
        Vector2::new(0.0, 2.8),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    assert!(matches!(body.envelope_state, EnvelopeState::None));
}

#[test]
fn envelope_captured_when_core_sufficient() {
    let disk = GasDisk::mmsn();
    let mut body = DiscreteBody::new(
        Mass::from_earth_masses(1.0), // Sufficient for envelope
        Mass::zero(),
        Point2::new(5.0, 0.0),
        Vector2::new(0.0, 2.8),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    assert!(matches!(body.envelope_state, EnvelopeState::None));

    // Evolve with typical conditions
    let dt = Time::from_years(1000.0);
    let core_accretion_rate = MassRate::from_earth_masses_per_myr(0.1);
    let opacity = 1.0;

    body.evolve_envelope(&disk, core_accretion_rate, dt, opacity);

    // Should now have small envelope
    assert!(matches!(
        body.envelope_state,
        EnvelopeState::Hydrostatic { .. }
    ));
    assert!(body.envelope_mass > Mass::zero());
}

#[test]
fn hydrostatic_envelope_grows_slowly() {
    let disk = GasDisk::mmsn();
    let initial_envelope = Mass::from_earth_masses(0.1);
    let mut body = DiscreteBody::new(
        Mass::from_earth_masses(3.0), // Below critical mass
        initial_envelope,
        Point2::new(5.0, 0.0),
        Vector2::new(0.0, 2.8),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    body.envelope_state = EnvelopeState::Hydrostatic {
        envelope_mass: initial_envelope,
    };

    // Evolve for 10,000 years
    let dt = Time::from_years(10_000.0);
    let core_accretion_rate = MassRate::from_earth_masses_per_myr(0.1);
    let opacity = 1.0;

    body.evolve_envelope(&disk, core_accretion_rate, dt, opacity);

    // Envelope should have grown
    assert!(body.envelope_mass > initial_envelope);

    // But growth should be slow (still hydrostatic)
    assert!(matches!(
        body.envelope_state,
        EnvelopeState::Hydrostatic { .. }
    ));

    // Growth should be modest over 10 kyr
    let growth_fraction = (body.envelope_mass - initial_envelope).to_earth_masses()
        / initial_envelope.to_earth_masses();
    assert!(growth_fraction < 1.0); // Less than doubled
}

#[test]
fn runaway_triggered_above_critical_mass() {
    let disk = GasDisk::mmsn();
    let mut body = DiscreteBody::new(
        Mass::from_earth_masses(15.0), // Above typical critical mass (~10 M⊕)
        Mass::from_earth_masses(1.0),
        Point2::new(5.0, 0.0),
        Vector2::new(0.0, 2.8),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    body.envelope_state = EnvelopeState::Hydrostatic {
        envelope_mass: Mass::from_earth_masses(1.0),
    };

    // Evolve with conditions that give low critical mass
    let dt = Time::from_years(1000.0);
    let core_accretion_rate = MassRate::from_earth_masses_per_myr(0.01); // Low rate → low M_crit
    let opacity = 0.1; // Low opacity → low M_crit

    body.evolve_envelope(&disk, core_accretion_rate, dt, opacity);

    // Should have transitioned to runaway
    assert!(matches!(body.envelope_state, EnvelopeState::Runaway { .. }));
}

#[test]
fn runaway_accretion_is_rapid() {
    let disk = GasDisk::mmsn();
    let initial_envelope = Mass::from_earth_masses(5.0);
    let mut body = DiscreteBody::new(
        Mass::from_earth_masses(15.0),
        initial_envelope,
        Point2::new(5.0, 0.0),
        Vector2::new(0.0, 2.8),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    body.envelope_state = EnvelopeState::Runaway {
        envelope_mass: initial_envelope,
        accretion_rate: MassRate::from_solar_masses_per_year(1e-6),
    };

    // Evolve for 10,000 years
    let dt = Time::from_years(10_000.0);
    let core_accretion_rate = MassRate::from_earth_masses_per_myr(0.1);
    let opacity = 1.0;

    body.evolve_envelope(&disk, core_accretion_rate, dt, opacity);

    // Envelope should have grown significantly
    let growth_factor = body.envelope_mass.to_earth_masses() / initial_envelope.to_earth_masses();
    assert!(growth_factor > 2.0); // Should at least double
}

#[test]
fn gap_opening_stops_runaway() {
    let disk = GasDisk::mmsn();
    // Create a massive body that will open gap
    let mut body = DiscreteBody::new(
        Mass::from_earth_masses(100.0), // Very massive core
        Mass::from_jupiter_masses(0.5), // Substantial envelope
        Point2::new(5.0, 0.0),
        Vector2::new(0.0, 2.8),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    body.envelope_state = EnvelopeState::Runaway {
        envelope_mass: Mass::from_jupiter_masses(0.5),
        accretion_rate: MassRate::from_solar_masses_per_year(1e-6),
    };

    // Evolve
    let dt = Time::from_years(1000.0);
    let core_accretion_rate = MassRate::from_earth_masses_per_myr(0.1);
    let opacity = 1.0;

    body.evolve_envelope(&disk, core_accretion_rate, dt, opacity);

    // Should have stopped (Final state)
    assert!(matches!(body.envelope_state, EnvelopeState::Final { .. }));
}

#[test]
fn gap_opening_criterion_for_jupiter() {
    let disk = GasDisk::mmsn();
    let jupiter = DiscreteBody::new(
        Mass::from_jupiter_masses(0.3), // Core
        Mass::from_jupiter_masses(0.7), // Envelope
        Point2::new(5.2, 0.0),
        Vector2::new(0.0, 2.755),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    let stellar_mass = Mass::from_solar_masses(1.0);

    // Jupiter-mass planet should open gap
    assert!(jupiter.opens_gap(&disk, stellar_mass));
}

#[test]
fn sub_neptune_does_not_open_gap() {
    let disk = GasDisk::mmsn();
    let sub_neptune = DiscreteBody::new(
        Mass::from_earth_masses(5.0), // Core
        Mass::from_earth_masses(5.0), // Envelope (10 M⊕ total)
        Point2::new(5.0, 0.0),
        Vector2::new(0.0, 2.8),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    let stellar_mass = Mass::from_solar_masses(1.0);

    // Sub-Neptune should NOT open gap
    assert!(!sub_neptune.opens_gap(&disk, stellar_mass));
}

#[test]
fn gap_criterion_depends_on_mass_ratio() {
    let disk = GasDisk::mmsn();
    let stellar_mass = Mass::from_solar_masses(1.0);

    // Small planet
    let small = DiscreteBody::new(
        Mass::from_earth_masses(1.0),
        Mass::zero(),
        Point2::new(5.0, 0.0),
        Vector2::new(0.0, 2.8),
        stellar_mass,
        Composition::earth_like(),
    );

    // Large planet
    let large = DiscreteBody::new(
        Mass::from_jupiter_masses(1.0),
        Mass::zero(),
        Point2::new(5.0, 0.0),
        Vector2::new(0.0, 2.8),
        stellar_mass,
        Composition::earth_like(),
    );

    assert!(!small.opens_gap(&disk, stellar_mass));
    assert!(large.opens_gap(&disk, stellar_mass));
}

#[test]
fn local_gas_mass_scales_with_hill_radius() {
    let disk = GasDisk::mmsn();

    // Small body (small Hill radius)
    let small = DiscreteBody::new(
        Mass::from_earth_masses(1.0),
        Mass::zero(),
        Point2::new(5.0, 0.0),
        Vector2::new(0.0, 2.8),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    // Large body (large Hill radius)
    let large = DiscreteBody::new(
        Mass::from_earth_masses(100.0),
        Mass::zero(),
        Point2::new(5.0, 0.0),
        Vector2::new(0.0, 2.8),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    let m_gas_small = small.local_gas_mass(&disk);
    let m_gas_large = large.local_gas_mass(&disk);

    // Larger body has more gas available
    assert!(m_gas_large > m_gas_small);

    // Local gas mass ∝ R_H² ∝ M^(2/3)
    let mass_ratio = large.total_mass().to_solar_masses() / small.total_mass().to_solar_masses();
    let gas_ratio = m_gas_large.to_solar_masses() / m_gas_small.to_solar_masses();
    let expected_ratio = mass_ratio.powf(2.0 / 3.0);

    assert!((gas_ratio / expected_ratio - 1.0).abs() < 0.1);
}

#[test]
fn runaway_stops_when_local_gas_depleted() {
    let disk = GasDisk::mmsn();

    // Create body with envelope nearly equal to local gas mass
    let initial_envelope = Mass::from_earth_masses(10.0);
    let mut body = DiscreteBody::new(
        Mass::from_earth_masses(15.0),
        initial_envelope,
        Point2::new(5.0, 0.0),
        Vector2::new(0.0, 2.8),
        Mass::from_solar_masses(1.0),
        Composition::earth_like(),
    );

    let local_gas = body.local_gas_mass(&disk);

    // Set envelope to just over half the local gas
    let large_envelope = local_gas * 0.6;
    body.envelope_mass = large_envelope;
    body.envelope_state = EnvelopeState::Runaway {
        envelope_mass: large_envelope,
        accretion_rate: MassRate::from_solar_masses_per_year(1e-6),
    };

    // Evolve
    let dt = Time::from_years(1000.0);
    let core_accretion_rate = MassRate::from_earth_masses_per_myr(0.1);
    let opacity = 1.0;

    body.evolve_envelope(&disk, core_accretion_rate, dt, opacity);

    // Should have stopped (consumed most local gas)
    assert!(matches!(body.envelope_state, EnvelopeState::Final { .. }));
}
