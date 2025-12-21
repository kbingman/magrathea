//! Orbital element representations and conversions.
//!
//! Provides conversion between Cartesian state (position, velocity) and
//! Keplerian orbital elements (a, e, i, Ω, ω, M).
//!
//! # Coordinate Systems
//!
//! - **Cartesian**: Position (x, y) in AU, velocity (vx, vy) in AU/year
//! - **Orbital Elements**: Semi-major axis, eccentricity, inclination, etc.
//!
//! # Usage
//!
//! Cartesian coordinates are used for numerical integration (no singularities).
//! Orbital elements are computed when needed for physics calculations.

use nalgebra::{Point2, Vector2};
use units::{Length, Mass};

/// Gravitational constant in AU³ M☉⁻¹ year⁻²
/// G = 4π² ≈ 39.478417
const G: f64 = 39.478417;

/// Keplerian orbital elements.
///
/// Describes an orbit around a central body. For 2D simulations,
/// inclination and related angles are zero or undefined.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrbitalElements {
    /// Semi-major axis (AU)
    pub semi_major_axis: Length,

    /// Eccentricity (dimensionless, 0 ≤ e < 1 for bound orbits)
    pub eccentricity: f64,

    /// Inclination (radians, always 0 for 2D)
    pub inclination: f64,

    /// Longitude of ascending node (radians, undefined in 2D)
    pub longitude_ascending_node: f64,

    /// Argument of periapsis (radians)
    pub argument_of_periapsis: f64,

    /// Mean anomaly (radians)
    pub mean_anomaly: f64,
}

impl OrbitalElements {
    /// Orbital period for this orbit.
    ///
    /// Uses Kepler's third law: T = 2π√(a³/GM)
    ///
    /// # Arguments
    /// * `stellar_mass` - Mass of central star
    ///
    /// # Returns
    /// Orbital period in years
    pub fn period(&self, stellar_mass: Mass) -> f64 {
        let a_au = self.semi_major_axis.to_au();
        let m_solar = stellar_mass.to_solar_masses();
        std::f64::consts::TAU * (a_au.powi(3) / (G * m_solar)).sqrt()
    }

    /// Orbital angular frequency (mean motion).
    ///
    /// n = 2π/T = √(GM/a³)
    ///
    /// # Arguments
    /// * `stellar_mass` - Mass of central star
    ///
    /// # Returns
    /// Angular frequency in radians/year
    pub fn mean_motion(&self, stellar_mass: Mass) -> f64 {
        let a_au = self.semi_major_axis.to_au();
        let m_solar = stellar_mass.to_solar_masses();
        (G * m_solar / a_au.powi(3)).sqrt()
    }
}

/// Convert Keplerian orbital elements to Cartesian state.
///
/// Converts orbital elements back to position and velocity vectors.
/// Uses iterative solution of Kepler's equation for eccentric anomaly.
///
/// # Arguments
/// * `elements` - Orbital elements
/// * `stellar_mass` - Mass of central star (solar masses)
///
/// # Returns
/// Tuple of (position, velocity) in AU and AU/year
///
/// # Examples
/// ```
/// use nalgebra::{Point2, Vector2};
/// use units::{Length, Mass};
/// use protodisk::bodies::{OrbitalElements, orbital_elements_to_cartesian};
///
/// let elements = OrbitalElements {
///     semi_major_axis: Length::from_au(1.0),
///     eccentricity: 0.2,
///     inclination: 0.0,
///     longitude_ascending_node: 0.0,
///     argument_of_periapsis: 0.0,
///     mean_anomaly: 0.0,
/// };
///
/// let (pos, vel) = orbital_elements_to_cartesian(&elements, Mass::from_solar_masses(1.0));
/// ```
pub fn orbital_elements_to_cartesian(
    elements: &OrbitalElements,
    stellar_mass: Mass,
) -> (Point2<f64>, Vector2<f64>) {
    let mu = G * stellar_mass.to_solar_masses();
    let a = elements.semi_major_axis.to_au();
    let e = elements.eccentricity;
    let omega = elements.argument_of_periapsis;
    let mean_anomaly = elements.mean_anomaly;

    // Solve Kepler's equation: M = E - e sin(E) for eccentric anomaly E
    let ecc_anomaly = solve_keplers_equation(mean_anomaly, e);

    // Convert eccentric anomaly to true anomaly
    let true_anomaly = if e < 0.999 {
        // Elliptical orbit
        let cos_nu = (ecc_anomaly.cos() - e) / (1.0 - e * ecc_anomaly.cos());
        let sin_nu = ((1.0 - e * e).sqrt() * ecc_anomaly.sin()) / (1.0 - e * ecc_anomaly.cos());
        sin_nu.atan2(cos_nu)
    } else {
        // Near-parabolic, use mean anomaly as approximation
        mean_anomaly
    };

    // Radius in orbit
    let r = a * (1.0 - e * ecc_anomaly.cos());

    // Position in orbital plane (periapsis at x-axis)
    let x_orb = r * true_anomaly.cos();
    let y_orb = r * true_anomaly.sin();

    // Velocity in orbital plane
    let v_factor = (mu / (a * (1.0 - e * e))).sqrt();
    let vx_orb = -v_factor * true_anomaly.sin();
    let vy_orb = v_factor * (e + true_anomaly.cos());

    // Rotate by argument of periapsis (2D rotation)
    let cos_omega = omega.cos();
    let sin_omega = omega.sin();

    let x = cos_omega * x_orb - sin_omega * y_orb;
    let y = sin_omega * x_orb + cos_omega * y_orb;

    let vx = cos_omega * vx_orb - sin_omega * vy_orb;
    let vy = sin_omega * vx_orb + cos_omega * vy_orb;

    (Point2::new(x, y), Vector2::new(vx, vy))
}

/// Solve Kepler's equation using Newton-Raphson iteration.
///
/// Solves M = E - e sin(E) for E (eccentric anomaly).
///
/// # Arguments
/// * `mean_anomaly` - Mean anomaly M (radians)
/// * `eccentricity` - Eccentricity e
///
/// # Returns
/// Eccentric anomaly E (radians)
pub(crate) fn solve_keplers_equation(mean_anomaly: f64, eccentricity: f64) -> f64 {
    // Initial guess: E₀ = M + e sin(M)
    let mut ecc_anomaly = mean_anomaly + eccentricity * mean_anomaly.sin();

    // Newton-Raphson: E_{n+1} = E_n - f(E_n)/f'(E_n)
    // where f(E) = E - e sin(E) - M
    // and f'(E) = 1 - e cos(E)
    for _ in 0..10 {
        // Usually converges in 3-4 iterations
        let f = ecc_anomaly - eccentricity * ecc_anomaly.sin() - mean_anomaly;
        let f_prime = 1.0 - eccentricity * ecc_anomaly.cos();

        let delta = f / f_prime;
        ecc_anomaly -= delta;

        if delta.abs() < 1e-10 {
            break;
        }
    }

    ecc_anomaly
}

/// Convert Cartesian state to Keplerian orbital elements.
///
/// Computes orbital elements from position and velocity vectors in 2D.
/// Uses the standard algorithm from orbital mechanics textbooks.
///
/// # Arguments
/// * `position` - Position vector (AU)
/// * `velocity` - Velocity vector (AU/year)
/// * `stellar_mass` - Mass of central star (solar masses)
///
/// # Returns
/// Orbital elements for the orbit
///
/// # Examples
/// ```
/// use nalgebra::{Point2, Vector2};
/// use units::Mass;
/// use protodisk::bodies::cartesian_to_orbital_elements;
///
/// // Circular orbit at 1 AU
/// let pos = Point2::new(1.0, 0.0);
/// let vel = Vector2::new(0.0, 6.28); // ~2π AU/year
/// let elements = cartesian_to_orbital_elements(pos, vel, Mass::from_solar_masses(1.0));
///
/// assert!((elements.semi_major_axis.to_au() - 1.0).abs() < 0.01);
/// assert!(elements.eccentricity < 0.01);
/// ```
pub fn cartesian_to_orbital_elements(
    position: Point2<f64>,
    velocity: Vector2<f64>,
    stellar_mass: Mass,
) -> OrbitalElements {
    let mu = G * stellar_mass.to_solar_masses(); // GM

    let r_vec = position.coords;
    let v_vec = velocity;

    let r = r_vec.magnitude();
    let v = v_vec.magnitude();

    // Specific orbital energy: ε = v²/2 - μ/r
    let specific_energy = v * v / 2.0 - mu / r;

    // Semi-major axis: a = -μ/(2ε)
    let a = if specific_energy.abs() < 1e-10 {
        // Nearly parabolic, use large value
        1e6
    } else {
        -mu / (2.0 * specific_energy)
    };

    // Specific angular momentum (2D, z-component only)
    // h = r × v = x*vy - y*vx
    let h = r_vec.x * v_vec.y - r_vec.y * v_vec.x;

    // Eccentricity vector: e_vec = (v × h)/μ - r/|r|
    // In 2D: e_vec = (v_vec × h_z) / μ - r_vec / r
    // Cross product in 2D: (vx, vy) × h_z = (vy*h, -vx*h)
    let e_vec = Vector2::new(
        v_vec.y * h / mu - r_vec.x / r,
        -v_vec.x * h / mu - r_vec.y / r,
    );

    let e = e_vec.magnitude();

    // Argument of periapsis: angle of eccentricity vector from x-axis
    let omega = if e > 1e-8 {
        e_vec.y.atan2(e_vec.x)
    } else {
        0.0 // Undefined for circular orbits
    };

    // True anomaly: angle of position vector from periapsis
    let true_anomaly = if e > 1e-8 {
        let cos_nu = e_vec.dot(&r_vec) / (e * r);
        // Use 2D cross product for sin: e × r = e_x*r_y - e_y*r_x
        let sin_nu = (e_vec.x * r_vec.y - e_vec.y * r_vec.x) / (e * r);
        sin_nu.atan2(cos_nu)
    } else {
        // For circular orbits, use angle from x-axis
        r_vec.y.atan2(r_vec.x)
    };

    // Eccentric anomaly (from true anomaly)
    let ecc_anomaly = if e < 0.999 {
        // Elliptical orbit
        let cos_e = (e + true_anomaly.cos()) / (1.0 + e * true_anomaly.cos());
        let sin_e = ((1.0 - e * e).sqrt() * true_anomaly.sin()) / (1.0 + e * true_anomaly.cos());
        sin_e.atan2(cos_e)
    } else {
        true_anomaly // Nearly parabolic
    };

    // Mean anomaly: M = E - e sin(E)
    let mean_anomaly = ecc_anomaly - e * ecc_anomaly.sin();

    OrbitalElements {
        semi_major_axis: Length::from_au(a),
        eccentricity: e.clamp(0.0, 0.9999), // Clamp to valid range
        inclination: 0.0,                   // 2D simulation
        longitude_ascending_node: 0.0,      // Undefined in 2D
        argument_of_periapsis: omega,
        mean_anomaly,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn circular_orbit_at_1_au() {
        // Circular orbit: v = sqrt(GM/r) = sqrt(39.478417) ≈ 6.283
        let pos = Point2::new(1.0, 0.0);
        let vel = Vector2::new(0.0, (G * 1.0 / 1.0).sqrt());

        let elements = cartesian_to_orbital_elements(pos, vel, Mass::from_solar_masses(1.0));

        assert!((elements.semi_major_axis.to_au() - 1.0).abs() < 0.01);
        assert!(elements.eccentricity < 0.01);
    }

    #[test]
    fn eccentric_orbit() {
        // Start at periapsis (r = a(1-e)) with velocity for e=0.3 orbit
        let a = 1.5;
        let e = 0.3;
        let r_peri = a * (1.0 - e);

        // At periapsis: v = sqrt(GM(1+e)/(a(1-e)))
        let mu = G * 1.0;
        let v_peri = (mu * (1.0 + e) / (a * (1.0 - e))).sqrt();

        let pos = Point2::new(r_peri, 0.0);
        let vel = Vector2::new(0.0, v_peri);

        let elements = cartesian_to_orbital_elements(pos, vel, Mass::from_solar_masses(1.0));

        assert!((elements.semi_major_axis.to_au() - a).abs() < 0.01);
        assert!((elements.eccentricity - e).abs() < 0.01);
    }

    #[test]
    fn period_matches_keplers_third_law() {
        let elements = OrbitalElements {
            semi_major_axis: Length::from_au(1.0),
            eccentricity: 0.0,
            inclination: 0.0,
            longitude_ascending_node: 0.0,
            argument_of_periapsis: 0.0,
            mean_anomaly: 0.0,
        };

        let period = elements.period(Mass::from_solar_masses(1.0));

        // T = 2π√(a³/GM) = 2π√(1/39.478417) ≈ 1.0 year
        assert!((period - 1.0).abs() < 0.01);
    }

    #[test]
    fn mean_motion_is_two_pi_over_period() {
        let elements = OrbitalElements {
            semi_major_axis: Length::from_au(1.0),
            eccentricity: 0.0,
            inclination: 0.0,
            longitude_ascending_node: 0.0,
            argument_of_periapsis: 0.0,
            mean_anomaly: 0.0,
        };

        let n = elements.mean_motion(Mass::from_solar_masses(1.0));
        let period = elements.period(Mass::from_solar_masses(1.0));

        let expected = std::f64::consts::TAU / period;

        assert!((n - expected).abs() < 0.01);
    }
}
