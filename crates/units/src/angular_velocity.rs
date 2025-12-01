/// Angular velocity in rad/s
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, serde::Serialize, serde::Deserialize)]
pub struct AngularVelocity(pub f64);

impl AngularVelocity {
    pub fn from_rad_per_sec(value: f64) -> Self {
        Self(value)
    }

    pub fn to_rad_per_sec(&self) -> f64 {
        self.0
    }
}
