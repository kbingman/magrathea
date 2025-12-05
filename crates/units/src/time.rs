use serde::{Deserialize, Serialize};
use std::ops::{Add, Div, Mul, Sub};

pub(crate) const DAYS_PER_YEAR: f64 = 365.25;
pub(crate) const HOURS_PER_YEAR: f64 = 8_766.0; // 365.25 * 24
pub const SECONDS_PER_YEAR: f64 = 31_557_600.0; // 365.25 days per year

/// Million years in regular years
const MYR_TO_YEARS: f64 = 1_000_000.0;

/// A physical time quantity using f64 precision.
///
/// The `Time` struct represents time with years as the base unit,
/// which is natural for stellar and planetary evolution timescales.
///
/// # Examples
///
/// ```rust
/// use units::Time;
///
/// // Create times in different units
/// let orbital_period = Time::from_years(1.0);
/// let rotation = Time::from_days(1.0);
/// let evolution = Time::from_myr(100.0);
///
/// // Convert between units
/// let days = orbital_period.to_days();
/// let myr = evolution.to_myr();
/// ```
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Time(f64); // Base unit: Years

impl Time {
    /// Creates a zero time value
    pub fn zero() -> Self {
        Self(0.0)
    }

    /// Creates a new `Time` from a value in years.
    pub fn from_years(value: f64) -> Self {
        Self(value)
    }

    /// Creates a time from a value in million years (Myr)
    pub fn from_myr(value: f64) -> Self {
        Self(value * MYR_TO_YEARS)
    }

    /// Creates a new `Time` from a value in days.
    pub fn from_days(value: f64) -> Self {
        Self(value / DAYS_PER_YEAR)
    }

    /// Creates a new `Time` from a value in hours.
    pub fn from_hours(value: f64) -> Self {
        Self(value / HOURS_PER_YEAR)
    }

    /// Creates a new `Time` from a value in seconds.
    pub fn from_seconds(value: f64) -> Self {
        Self(value / SECONDS_PER_YEAR)
    }

    /// Returns the time in years.
    pub fn to_years(&self) -> f64 {
        self.0
    }

    /// Returns the time in million years
    pub fn to_myr(&self) -> f64 {
        self.0 / MYR_TO_YEARS
    }

    /// Converts the time to days.
    pub fn to_days(&self) -> f64 {
        self.0 * DAYS_PER_YEAR
    }

    /// Converts the time to hours.
    pub fn to_hours(&self) -> f64 {
        self.0 * HOURS_PER_YEAR
    }

    /// Converts the time to seconds.
    pub fn to_seconds(&self) -> f64 {
        self.0 * SECONDS_PER_YEAR
    }
}

impl Add for Time {
    type Output = Time;

    fn add(self, rhs: Time) -> Time {
        Time(self.0 + rhs.0)
    }
}

impl Sub for Time {
    type Output = Time;

    fn sub(self, rhs: Time) -> Time {
        Time(self.0 - rhs.0)
    }
}

impl Mul<f64> for Time {
    type Output = Time;

    fn mul(self, rhs: f64) -> Time {
        Time(self.0 * rhs)
    }
}

impl Div<f64> for Time {
    type Output = Time;

    fn div(self, rhs: f64) -> Time {
        Time(self.0 / rhs)
    }
}

/// Division of Time by Time returns a dimensionless ratio
impl Div for Time {
    type Output = f64;

    fn div(self, rhs: Time) -> f64 {
        self.0 / rhs.0
    }
}

impl Div<&Time> for Time {
    type Output = f64;

    fn div(self, rhs: &Time) -> f64 {
        self.0 / rhs.0
    }
}

impl Div<&Time> for &Time {
    type Output = f64;

    fn div(self, rhs: &Time) -> f64 {
        self.0 / rhs.0
    }
}

/// Allow f64 * Time (commutative multiplication)
impl Mul<Time> for f64 {
    type Output = Time;

    fn mul(self, rhs: Time) -> Time {
        rhs * self
    }
}
