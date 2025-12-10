use nalgebra::{Point2, Vector2};
use units::{Length, Mass};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BodyId(pub u32);

#[derive(Debug, Clone, Copy)]
pub struct Body {
    pub id: BodyId,
    pub mass: f64,              // Solar masses
    pub radius: f64,            // AU (physical radius for collisions)
    pub position: Point2<f64>,  // AU (heliocentric Cartesian, 2D)
    pub velocity: Vector2<f64>, // AU/year
}

impl Body {
    /// Creates a new body with mass specified in solar masses (for doctest examples)
    ///
    /// This is a convenience constructor for doctests and examples.
    /// Position and velocity are in AU and AU/year respectively.
    pub fn new_solar_masses(mass_solar: f64, position: [f64; 2], velocity: [f64; 2]) -> Self {
        Body {
            id: BodyId(0),
            mass: mass_solar,
            radius: Length::from_solar_radii(1.0).to_au(), // 1 solar radius in AU
            position: Point2::new(position[0], position[1]),
            velocity: Vector2::new(velocity[0], velocity[1]),
        }
    }

    /// Creates a new body with mass specified in Earth masses (for doctest examples)
    ///
    /// This is a convenience constructor for doctests and examples.
    /// Position and velocity are in AU and AU/year respectively.
    pub fn new_earth_masses(mass_earth: f64, position: [f64; 2], velocity: [f64; 2]) -> Self {
        let mass = Mass::from_earth_masses(mass_earth);

        Body {
            id: BodyId(0),
            mass: mass.to_solar_masses(),
            radius: Length::from_earth_radii(1.0).to_au(), // 1 Earth radius in AU
            position: Point2::new(position[0], position[1]),
            velocity: Vector2::new(velocity[0], velocity[1]),
        }
    }

    pub fn momentum(&self) -> Vector2<f64> {
        self.velocity * self.mass
    }

    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * self.velocity.magnitude_squared()
    }

    pub fn distance_to(&self, other: &Body) -> f64 {
        (self.position - other.position).magnitude()
    }

    pub fn orbital_radius(&self) -> f64 {
        self.position.coords.magnitude()
    }

    pub fn orbital_velocity(&self) -> f64 {
        self.velocity.magnitude()
    }

    /// Angular momentum scalar (r × v, not multiplied by mass)
    /// In 2D, this returns the z-component of the angular momentum vector
    pub fn specific_angular_momentum(&self) -> f64 {
        self.position.x * self.velocity.y - self.position.y * self.velocity.x
    }
}

// Implement Massive trait for Body to work with Barnes-Hut tree
impl crate::arena_bhtree::Massive for Body {
    fn position(&self) -> Point2<f64> {
        self.position
    }

    fn mass(&self) -> Mass {
        Mass::from_solar_masses(self.mass)
    }

    fn mass_solar(&self) -> f64 {
        self.mass
    }
}
