//! Direct N-body gravity (O(N²) implementation)

use crate::forces::{ForceModel, G};
use crate::state::SystemState;
use nalgebra::Vector2;

/// Direct O(N²) gravitational force computation
///
/// Computes gravitational acceleration by summing forces from all other bodies
/// and the central star. Simple and accurate, but scales poorly for large N.
///
/// Best for:
/// - Small systems (N < 100)
/// - Testing and validation
/// - Benchmarking against tree methods
///
/// # Examples
///
/// ```
/// use nbody::forces::{DirectGravity, ForceModel};
/// use nbody::state::SystemState;
/// use stellar::generation::main_sequence_star;
/// use nalgebra::{Point2, Vector2};
/// use units::Length;
///
/// let star = main_sequence_star(1.0, 0.0, 4_600.0);
/// let mut system = SystemState::new(star);
/// system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
///
/// let gravity = DirectGravity::new();
/// let accel = gravity.acceleration(0, &system);
///
/// // Should point toward star (negative x direction)
/// assert!(accel.x < 0.0);
/// ```
pub struct DirectGravity {
    /// Optional softening length to prevent singularities (AU)
    pub softening: f64,
}

impl DirectGravity {
    /// Creates a new direct gravity force with no softening
    pub fn new() -> Self {
        Self { softening: 0.0 }
    }

    /// Creates a new direct gravity force with specified softening length
    ///
    /// # Arguments
    ///
    /// * `softening` - Softening length in AU
    ///
    /// # Examples
    ///
    /// ```
    /// use nbody::forces::DirectGravity;
    ///
    /// // Use 0.01 AU softening to prevent close encounter singularities
    /// let gravity = DirectGravity::with_softening(0.01);
    /// ```
    pub fn with_softening(softening: f64) -> Self {
        Self { softening }
    }
}

impl Default for DirectGravity {
    fn default() -> Self {
        Self::new()
    }
}

impl ForceModel for DirectGravity {
    fn acceleration(&self, idx: usize, state: &SystemState) -> Vector2<f64> {
        let body = &state.bodies[idx];
        let eps2 = self.softening * self.softening;

        // Acceleration from star (at origin)
        let r_star = body.position.coords; // Vector from origin to body
        let r2_star = r_star.magnitude_squared() + eps2;
        let r_star_mag = r2_star.sqrt();
        let star_mass = state.star.mass.to_solar_masses();
        let a_star = -r_star * (G * star_mass / (r2_star * r_star_mag));

        // Acceleration from other bodies
        let a_bodies = state
            .bodies
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != idx)
            .map(|(_, other)| {
                let dr = other.position - body.position; // This is already a Vector2
                let r2 = dr.magnitude_squared() + eps2;
                let r = r2.sqrt();
                dr * (G * other.mass / (r2 * r))
            })
            .fold(Vector2::zeros(), |acc, a| acc + a);

        a_star + a_bodies
    }

    fn potential_energy(&self, state: &SystemState) -> f64 {
        let eps2 = self.softening * self.softening;
        let star_mass = state.star.mass.to_solar_masses();

        // Star-body potential
        let star_potential: f64 = state
            .bodies
            .iter()
            .map(|b| {
                let r = (b.position.coords.magnitude_squared() + eps2).sqrt();
                -G * star_mass * b.mass / r
            })
            .sum();

        // Body-body potential (each pair counted once)
        let body_potential: f64 = state
            .bodies
            .iter()
            .enumerate()
            .flat_map(|(i, a)| {
                state.bodies[i + 1..].iter().map(move |b| {
                    let dr = a.position - b.position; // Vector2
                    let r = (dr.magnitude_squared() + eps2).sqrt();
                    -G * a.mass * b.mass / r
                })
            })
            .sum();

        star_potential + body_potential
    }
}

/// Compute Hill radius for a body orbiting a star
///
/// The Hill radius is the approximate radius of gravitational influence
/// of a smaller body in orbit around a larger body.
///
/// # Arguments
///
/// * `mass` - Body mass in solar masses
/// * `orbital_radius` - Orbital radius in AU
/// * `star_mass` - Star mass in solar masses
///
/// # Returns
///
/// Hill radius in AU
///
/// # Examples
///
/// ```
/// use nbody::forces::gravity::hill_radius;
///
/// // Earth's Hill radius at 1 AU
/// let r_hill = hill_radius(3.0e-6, 1.0, 1.0);
/// assert!((r_hill - 0.01).abs() < 0.001); // About 0.01 AU
/// ```
pub fn hill_radius(mass: f64, orbital_radius: f64, star_mass: f64) -> f64 {
    orbital_radius * (mass / (3.0 * star_mass)).powf(1.0 / 3.0)
}

/// Compute mutual Hill radius of two bodies
///
/// Uses the average of the two orbital radii and sum of masses.
///
/// # Arguments
///
/// * `body1_mass` - First body mass in solar masses
/// * `body1_radius` - First body orbital radius in AU
/// * `body2_mass` - Second body mass in solar masses
/// * `body2_radius` - Second body orbital radius in AU
/// * `star_mass` - Star mass in solar masses
///
/// # Returns
///
/// Mutual Hill radius in AU
///
/// # Examples
///
/// ```
/// use nbody::forces::gravity::mutual_hill_radius;
///
/// // Two Earth-mass planets at 1 and 1.5 AU
/// let r_hill = mutual_hill_radius(3.0e-6, 1.0, 3.0e-6, 1.5, 1.0);
/// ```
pub fn mutual_hill_radius(
    body1_mass: f64,
    body1_radius: f64,
    body2_mass: f64,
    body2_radius: f64,
    star_mass: f64,
) -> f64 {
    let a_avg = (body1_radius + body2_radius) / 2.0;
    let m_sum = body1_mass + body2_mass;

    a_avg * (m_sum / (3.0 * star_mass)).powf(1.0 / 3.0)
}
