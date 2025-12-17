//! Binary star orbital parameter types.

use serde::{Deserialize, Serialize};
use units::{Length, Time, Velocity};

/// Type of planetary orbit in a binary star system
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BinaryOrbitType {
    /// Planet orbits the primary star (circumstellar around primary)
    STypePrimary,

    /// Planet orbits the secondary star (circumstellar around secondary)
    STypeSecondary,

    /// Planet orbits both stars (circumbinary)
    PType,
}

/// Orbital parameters for a binary star system
///
/// These parameters define the Keplerian orbit of the two stars around
/// their common center of mass (barycenter).
///
/// # Coordinate System
/// The orbit is defined in a reference frame centered on the barycenter.
/// The orbital plane is the xy-plane, with the z-axis perpendicular.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct OrbitalParameters {
    /// Semi-major axis of the relative orbit (distance between stars at mean separation)
    pub semi_major_axis: Length,

    /// Orbital eccentricity (0 = circular, 0-1 = elliptical)
    pub eccentricity: f64,

    /// Inclination angle relative to the reference plane (0-180 degrees)
    /// 0° = orbit in xy-plane viewed face-on
    /// 90° = orbit viewed edge-on
    pub inclination: f64,

    /// Longitude of ascending node (0-360 degrees)
    /// Defines where orbit crosses reference plane
    pub longitude_of_ascending_node: f64,

    /// Argument of periapsis (0-360 degrees)
    /// Angle from ascending node to periapsis
    pub argument_of_periapsis: f64,

    /// Mean anomaly at epoch (0-360 degrees)
    /// Position of stars in orbit at reference time
    pub mean_anomaly: f64,

    /// Orbital period
    pub period: Time,
}

impl OrbitalParameters {
    /// Calculate the periapsis distance (closest approach)
    pub fn periapsis(&self) -> Length {
        Length::from_au(self.semi_major_axis.to_au() * (1.0 - self.eccentricity))
    }

    /// Calculate the apoapsis distance (farthest separation)
    pub fn apoapsis(&self) -> Length {
        Length::from_au(self.semi_major_axis.to_au() * (1.0 + self.eccentricity))
    }

    /// Calculate the mean orbital velocity
    ///
    /// Uses the approximation v_mean ≈ 2πa/P
    pub fn mean_orbital_velocity(&self) -> Velocity {
        let a_au = self.semi_major_axis.to_au();
        let period_years = self.period.to_years();

        Velocity::from_au_per_year(2.0 * std::f64::consts::PI * a_au / period_years)
    }
}

/// Binary star configuration for a planetary system
///
/// This type is stored in `SystemMetadata` to track binary orbital
/// parameters and planetary orbit type.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct BinaryConfiguration {
    /// Type of planetary orbits in this system
    pub orbit_type: BinaryOrbitType,

    /// Binary star orbital parameters
    pub orbital_params: OrbitalParameters,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_periapsis_apoapsis() {
        let orbit = OrbitalParameters {
            semi_major_axis: Length::from_au(1.0),
            eccentricity: 0.3,
            inclination: 0.0,
            longitude_of_ascending_node: 0.0,
            argument_of_periapsis: 0.0,
            mean_anomaly: 0.0,
            period: Time::from_years(1.0),
        };

        let periapsis = orbit.periapsis().to_au();
        let apoapsis = orbit.apoapsis().to_au();

        assert!((periapsis - 0.7).abs() < 1e-6);
        assert!((apoapsis - 1.3).abs() < 1e-6);
    }

    #[test]
    fn test_circular_orbit() {
        let orbit = OrbitalParameters {
            semi_major_axis: Length::from_au(5.0),
            eccentricity: 0.0,
            inclination: 0.0,
            longitude_of_ascending_node: 0.0,
            argument_of_periapsis: 0.0,
            mean_anomaly: 0.0,
            period: Time::from_years(10.0),
        };

        let periapsis = orbit.periapsis().to_au();
        let apoapsis = orbit.apoapsis().to_au();

        // Circular orbit: periapsis == apoapsis == semi_major_axis
        assert!((periapsis - 5.0).abs() < 1e-6);
        assert!((apoapsis - 5.0).abs() < 1e-6);
    }
}
