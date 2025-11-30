use serde::{Deserialize, Serialize};
use std::ops::{Add, Div, Mul, Sub};

use crate::mass::{Mass, EARTH_MASS_G, SOLAR_MASS_G};

/// A physical mass rate (mass per time) quantity using f64 precision.
///
/// The `MassRate` struct represents mass flow rates with solar masses per year as the base unit.
/// This is particularly useful for modeling accretion processes, stellar mass loss, and
/// disk evolution in planet formation simulations.
///
/// # Examples
///
/// ```rust
/// use units::mass_rate::MassRate;
/// use units::time::Time;
///
/// // Create mass rates using different units
/// let stellar_wind = MassRate::from_solar_masses_per_year(1e-14);  // Typical solar wind
/// let disk_accretion = MassRate::from_earth_masses_per_myr(60.0);  // Disk pebble flux
///
/// // Integrate over time to get total mass
/// let time_span = Time::from_years(1_000_000.0);
/// let accreted_mass = disk_accretion.integrate(time_span);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Deserialize, Serialize)]
pub struct MassRate(f64); // Base unit: Solar Masses per year

impl MassRate {
    /// Creates a new `MassRate` from a value in solar masses per year.
    ///
    /// This is the most direct constructor since M☉/year is the base unit.
    ///
    /// # Arguments
    ///
    /// * `value` - The mass rate in solar masses per year
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::mass_rate::MassRate;
    ///
    /// let stellar_wind = MassRate::from_solar_masses_per_year(1e-14);  // Solar wind loss rate
    /// let agn_accretion = MassRate::from_solar_masses_per_year(1.0);   // AGN accretion rate
    /// ```
    pub fn from_solar_masses_per_year(value: f64) -> Self {
        Self(value)
    }

    /// Creates a new `MassRate` from a value in Earth masses per million years.
    ///
    /// This unit is commonly used in planet formation simulations where
    /// timescales are millions of years and masses are planetary.
    ///
    /// # Arguments
    ///
    /// * `value` - The mass rate in Earth masses per million years
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::mass_rate::MassRate;
    ///
    /// let pebble_flux = MassRate::from_earth_masses_per_myr(60.0);     // Typical pebble accretion
    /// let gas_accretion = MassRate::from_earth_masses_per_myr(100.0);  // Gas giant formation
    /// ```
    pub fn from_earth_masses_per_myr(value: f64) -> Self {
        // Convert: (M⊕ / Myr) -> (M☉ / year)
        // = (M⊕ / Myr) * (M⊕/M☉) * (Myr/year)
        // = value * (EARTH_MASS_G / SOLAR_MASS_G) / 1e6
        Self(value * EARTH_MASS_G / SOLAR_MASS_G / 1e6)
    }

    /// Creates a new `MassRate` from a value in grams per year.
    ///
    /// # Arguments
    ///
    /// * `value` - The mass rate in grams per year
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::mass_rate::MassRate;
    ///
    /// let small_flux = MassRate::from_grams_per_year(1e20);  // Small object accretion
    /// ```
    pub fn from_grams_per_year(value: f64) -> Self {
        Self(value / SOLAR_MASS_G)
    }

    /// Returns the mass rate value in solar masses per year.
    ///
    /// # Returns
    ///
    /// The mass rate in M☉/year
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::mass_rate::MassRate;
    ///
    /// let rate = MassRate::from_earth_masses_per_myr(60.0);
    /// let solar_rate = rate.to_solar_masses_per_year();
    /// ```
    pub fn to_solar_masses_per_year(&self) -> f64 {
        self.0
    }

    /// Converts the mass rate to Earth masses per million years.
    ///
    /// # Returns
    ///
    /// The mass rate in M⊕/Myr
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::mass_rate::MassRate;
    ///
    /// let rate = MassRate::from_solar_masses_per_year(1e-10);
    /// let earth_rate = rate.to_earth_masses_per_myr();
    /// ```
    pub fn to_earth_masses_per_myr(&self) -> f64 {
        // Convert: (M☉ / year) -> (M⊕ / Myr)
        // = value * (M☉/M⊕) * (year/Myr)
        // = value * (SOLAR_MASS_G / EARTH_MASS_G) * 1e6
        self.0 * SOLAR_MASS_G / EARTH_MASS_G * 1e6
    }

    /// Converts the mass rate to grams per year.
    ///
    /// # Returns
    ///
    /// The mass rate in g/year
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::mass_rate::MassRate;
    ///
    /// let rate = MassRate::from_solar_masses_per_year(1e-14);
    /// let gram_rate = rate.to_grams_per_year();
    /// ```
    pub fn to_grams_per_year(&self) -> f64 {
        self.0 * SOLAR_MASS_G
    }

    /// Integrates the mass rate over a time period to get total mass.
    ///
    /// This is a convenience method for calculating total accreted/lost mass.
    ///
    /// # Arguments
    ///
    /// * `duration` - The time period
    ///
    /// # Returns
    ///
    /// The integrated mass
    ///
    /// # Examples
    ///
    /// ```rust
    /// use units::mass_rate::MassRate;
    /// use units::time::Time;
    ///
    /// let flux = MassRate::from_earth_masses_per_myr(60.0);
    /// let duration = Time::from_years(1_000_000.0);  // 1 Myr
    /// let total_mass = flux.integrate(duration);
    /// // total_mass should be approximately 60 Earth masses
    /// ```
    pub fn integrate(&self, duration: crate::time::Time) -> Mass {
        let years = duration.to_years();
        Mass::from_solar_masses(self.0 * years)
    }
}

impl Add for MassRate {
    type Output = MassRate;

    fn add(self, rhs: MassRate) -> MassRate {
        MassRate(self.0 + rhs.0)
    }
}

impl Sub for MassRate {
    type Output = MassRate;

    fn sub(self, rhs: MassRate) -> MassRate {
        MassRate(self.0 - rhs.0)
    }
}

impl Mul<f64> for MassRate {
    type Output = MassRate;

    fn mul(self, rhs: f64) -> MassRate {
        MassRate(self.0 * rhs)
    }
}

impl Div<f64> for MassRate {
    type Output = MassRate;

    fn div(self, rhs: f64) -> MassRate {
        MassRate(self.0 / rhs)
    }
}

/// Allow f64 * MassRate (commutative multiplication)
impl Mul<MassRate> for f64 {
    type Output = MassRate;

    fn mul(self, rhs: MassRate) -> MassRate {
        rhs * self
    }
}
