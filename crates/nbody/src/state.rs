use crate::body::{Body, BodyId};
use nalgebra::{Point2, Vector2};
use stellar::stellar_objects::MainSequenceStar;

/// Complete state of an N-body system at a given time
#[derive(Debug, Clone)]
pub struct SystemState {
    /// Current simulation time in years
    pub time: f64,
    /// Central star
    pub star: MainSequenceStar,
    /// Collection of orbiting bodies
    pub bodies: Vec<Body>,
    /// Next available body ID
    next_id: u32,
}

impl SystemState {
    /// Creates a new system with only a central star
    ///
    /// # Arguments
    ///
    /// * `star` - The central star
    ///
    /// # Examples
    ///
    /// ```
    /// use nbody::state::SystemState;
    /// use stellar::generation::main_sequence_star;
    ///
    /// // Generate a Sun-like star
    /// let star = main_sequence_star(1.0, 0.0, 4_600.0);
    /// let system = SystemState::new(star);
    ///
    /// assert_eq!(system.body_count(), 0);
    /// assert_eq!(system.time, 0.0);
    /// ```
    pub fn new(star: MainSequenceStar) -> Self {
        Self {
            time: 0.0,
            star,
            bodies: Vec::new(),
            next_id: 0,
        }
    }

    /// Adds a new body to the system and returns its ID
    ///
    /// # Arguments
    ///
    /// * `mass` - Body mass in solar masses
    /// * `radius` - Physical radius in AU
    /// * `position` - Position in AU (heliocentric)
    /// * `velocity` - Velocity in AU/year
    ///
    /// # Returns
    ///
    /// BodyId that can be used to reference this body
    ///
    /// # Examples
    ///
    /// ```
    /// use nbody::state::SystemState;
    /// use nalgebra::{Point2, Vector2};
    /// use units::{Mass, Length};
    /// use stellar::generation::main_sequence_star;
    ///
    /// let star = main_sequence_star(1.0, 0.0, 4_600.0);
    /// let mut system = SystemState::new(star);
    ///
    /// let earth_id = system.add_body(
    ///     Mass::from_earth_masses(1.0).to_solar_masses(),
    ///     Length::from_earth_radii(1.0).to_au(),
    ///     Point2::new(1.0, 0.0),
    ///     Vector2::new(0.0, 6.28),
    /// );
    ///
    /// assert_eq!(system.body_count(), 1);
    /// ```
    pub fn add_body(
        &mut self,
        mass: f64,
        radius: f64,
        position: Point2<f64>,
        velocity: Vector2<f64>,
    ) -> BodyId {
        let id = BodyId(self.next_id);
        self.next_id += 1;
        self.bodies.push(Body {
            id,
            mass,
            radius,
            position,
            velocity,
        });
        id
    }

    /// Removes a body from the system
    ///
    /// # Arguments
    ///
    /// * `id` - The BodyId to remove
    ///
    /// # Returns
    ///
    /// The removed Body if found, None otherwise
    ///
    /// # Examples
    ///
    /// ```
    /// use nbody::state::SystemState;
    /// use nalgebra::{Point2, Vector2};
    /// use stellar::generation::main_sequence_star;
    ///
    /// let star = main_sequence_star(1.0, 0.0, 4_600.0);
    /// let mut system = SystemState::new(star);
    ///
    /// let id = system.add_body(1.0, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
    /// assert_eq!(system.body_count(), 1);
    ///
    /// let removed = system.remove_body(id);
    /// assert!(removed.is_some());
    /// assert_eq!(system.body_count(), 0);
    /// ```
    pub fn remove_body(&mut self, id: BodyId) -> Option<Body> {
        self.bodies
            .iter()
            .position(|b| b.id == id)
            .map(|idx| self.bodies.remove(idx))
    }

    /// Gets a reference to a body by ID
    ///
    /// # Arguments
    ///
    /// * `id` - The BodyId to look up
    ///
    /// # Returns
    ///
    /// Reference to the Body if found, None otherwise
    pub fn get_body(&self, id: BodyId) -> Option<&Body> {
        self.bodies.iter().find(|b| b.id == id)
    }

    /// Gets a mutable reference to a body by ID
    ///
    /// # Arguments
    ///
    /// * `id` - The BodyId to look up
    ///
    /// # Returns
    ///
    /// Mutable reference to the Body if found, None otherwise
    pub fn get_body_mut(&mut self, id: BodyId) -> Option<&mut Body> {
        self.bodies.iter_mut().find(|b| b.id == id)
    }

    /// Returns the number of bodies in the system
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }

    /// Returns the total mass of all bodies (excluding the star)
    ///
    /// # Examples
    ///
    /// ```
    /// use nbody::state::SystemState;
    /// use nalgebra::{Point2, Vector2};
    /// use stellar::generation::main_sequence_star;
    ///
    /// let star = main_sequence_star(1.0, 0.0, 4_600.0);
    /// let mut system = SystemState::new(star);
    ///
    /// system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
    /// system.add_body(0.002, 0.01, Point2::new(2.0, 0.0), Vector2::new(0.0, 4.44));
    ///
    /// assert_eq!(system.total_planet_mass(), 0.003);
    /// ```
    pub fn total_planet_mass(&self) -> f64 {
        self.bodies.iter().map(|b| b.mass).sum()
    }

    /// Returns the total momentum of all bodies
    ///
    /// This should be approximately zero for an isolated system
    /// (useful for checking numerical drift)
    pub fn total_momentum(&self) -> Vector2<f64> {
        self.bodies
            .iter()
            .map(|b| b.momentum())
            .fold(Vector2::zeros(), |acc, p| acc + p)
    }

    /// Returns the total angular momentum of all bodies
    ///
    /// This should be conserved in an isolated system
    pub fn total_angular_momentum(&self) -> f64 {
        self.bodies
            .iter()
            .map(|b| b.specific_angular_momentum() * b.mass)
            .sum()
    }
}
