//! System metadata types for generation provenance and identification.

use serde::{Deserialize, Serialize};
use uuid::Uuid;

#[cfg(feature = "tsify")]
use tsify_next::Tsify;

use crate::architecture::SystemArchitecture;

/// How the planetary system was generated
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "tsify", derive(Tsify))]
#[cfg_attr(feature = "tsify", tsify(into_wasm_abi, from_wasm_abi))]
pub enum GenerationMethod {
    /// Fast occurrence-rate sampling (statistical generator)
    Statistical,

    /// Physics-based formation simulation (stellar-forge)
    StellarForge,

    /// User-defined or imported from external source
    Manual,

    /// Real observed system (e.g., from exoplanet archive)
    Observed,
}

/// Metadata about system generation and identification
///
/// Every planetary system has metadata that tracks:
/// - A unique identifier (UUID) that also serves as the RNG seed source
/// - How the system was generated
/// - The system's architectural classification
/// - A catalog designation (e.g., "KV-4729") derived from the UUID
/// - An optional proper name for notable "hero" systems
/// - Binary star configuration (if applicable)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
#[cfg_attr(feature = "tsify", derive(Tsify))]
#[cfg_attr(feature = "tsify", tsify(into_wasm_abi, from_wasm_abi))]
pub struct SystemMetadata {
    /// Unique identifier for this system (also used as RNG seed source)
    ///
    /// UUIDs are JSON-safe (serialized as strings) and avoid JavaScript's
    /// `Number.MAX_SAFE_INTEGER` limitation that corrupts large u64 values.
    pub id: Uuid,

    /// How this system was generated
    pub generation_method: GenerationMethod,

    /// System architecture classification
    pub architecture: SystemArchitecture,

    /// Short catalog designation (e.g., "KV-4729", "AN-0821")
    ///
    /// Format: Two uppercase letters + 4 digits, derived deterministically from the UUID.
    /// Provides ~6.76 million unique combinations (26² × 10000).
    pub catalog_name: String,

    /// Optional proper name for "hero" systems (e.g., "New Eden", "Cygnus Prime")
    ///
    /// Most systems use only the auto-generated `catalog_name` (e.g., "KV-4729").
    /// This field is for notable systems that deserve memorable names.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,

    /// Binary star configuration
    ///
    /// Present when the system contains a binary star pair (stars.len() == 2).
    /// Defines the orbital parameters of the binary and the type of planetary
    /// orbits (S-type around one star, or P-type circumbinary).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub binary_config: Option<crate::binary::BinaryConfiguration>,
}

/// Generate a catalog designation from a UUID
///
/// Format: Two uppercase letters + 4 digits (e.g., "KV-4729", "AN-0821")
/// Deterministic - same UUID always produces same designation.
/// Provides ~6.76 million unique combinations (26² × 10000).
fn catalog_name_from_uuid(id: &Uuid) -> String {
    let bytes = id.as_bytes();
    let prefix1 = (bytes[0] % 26 + b'A') as char;
    let prefix2 = (bytes[1] % 26 + b'A') as char;
    let number = u16::from_le_bytes([bytes[2], bytes[3]]) % 10000;
    format!("{}{}-{:04}", prefix1, prefix2, number)
}

impl SystemMetadata {
    /// Derive a u64 seed from the UUID for RNG initialization
    ///
    /// Uses the first 8 bytes of the UUID, providing deterministic
    /// seed generation from any UUID.
    ///
    /// # Example
    /// ```
    /// use star_system::{SystemMetadata, GenerationMethod, SystemArchitecture};
    ///
    /// let meta = SystemMetadata::new_random(
    ///     GenerationMethod::Statistical,
    ///     SystemArchitecture::Mixed,
    /// );
    /// let seed = meta.seed();
    /// // seed is a deterministic u64 derived from the UUID
    /// ```
    pub fn seed(&self) -> u64 {
        self.id.as_u64_pair().0
    }

    /// Returns the display name: proper name if set, otherwise catalog name
    ///
    /// # Example
    /// ```
    /// use star_system::{SystemMetadata, GenerationMethod, SystemArchitecture};
    ///
    /// // Without proper name - uses catalog name
    /// let meta = SystemMetadata::new_random(
    ///     GenerationMethod::Statistical,
    ///     SystemArchitecture::Mixed,
    /// );
    /// assert_eq!(meta.display_name(), meta.catalog_name);
    ///
    /// // With proper name
    /// let meta = meta.with_name("Cygnus Prime");
    /// assert_eq!(meta.display_name(), "Cygnus Prime");
    /// ```
    pub fn display_name(&self) -> String {
        self.name
            .clone()
            .unwrap_or_else(|| self.catalog_name.clone())
    }

    /// Create metadata with a random UUID
    ///
    /// # Example
    /// ```
    /// use star_system::{SystemMetadata, GenerationMethod, SystemArchitecture};
    ///
    /// let meta = SystemMetadata::new_random(
    ///     GenerationMethod::Statistical,
    ///     SystemArchitecture::CompactMulti,
    /// );
    /// ```
    pub fn new_random(
        generation_method: GenerationMethod,
        architecture: SystemArchitecture,
    ) -> Self {
        let id = Uuid::new_v4();
        Self {
            catalog_name: catalog_name_from_uuid(&id),
            id,
            generation_method,
            architecture,
            name: None,
            binary_config: None,
        }
    }

    /// Create metadata with a specific UUID
    ///
    /// Useful when you have an existing UUID or want to restore a system
    /// from a previously saved state.
    pub fn with_id(
        id: Uuid,
        generation_method: GenerationMethod,
        architecture: SystemArchitecture,
    ) -> Self {
        Self {
            catalog_name: catalog_name_from_uuid(&id),
            id,
            generation_method,
            architecture,
            name: None,
            binary_config: None,
        }
    }

    /// Create metadata with a deterministic UUID derived from a seed string
    ///
    /// Useful for reproducible generation from human-readable identifiers.
    /// The same seed_name always produces the same UUID (and thus the same RNG seed).
    ///
    /// Note: This does NOT set the display name - use `with_name()` for that.
    ///
    /// # Example
    /// ```
    /// use star_system::{SystemMetadata, GenerationMethod, SystemArchitecture};
    ///
    /// let meta1 = SystemMetadata::from_seed_name(
    ///     "test-system-42",
    ///     GenerationMethod::Statistical,
    ///     SystemArchitecture::Mixed,
    /// );
    /// let meta2 = SystemMetadata::from_seed_name(
    ///     "test-system-42",
    ///     GenerationMethod::Statistical,
    ///     SystemArchitecture::Mixed,
    /// );
    /// // Same seed name produces same UUID
    /// assert_eq!(meta1.id, meta2.id);
    /// assert_eq!(meta1.seed(), meta2.seed());
    /// ```
    pub fn from_seed_name(
        seed_name: &str,
        generation_method: GenerationMethod,
        architecture: SystemArchitecture,
    ) -> Self {
        let id = Uuid::new_v5(&Uuid::NAMESPACE_OID, seed_name.as_bytes());
        Self {
            catalog_name: catalog_name_from_uuid(&id),
            id,
            generation_method,
            architecture,
            name: None,
            binary_config: None,
        }
    }

    /// Set a proper name for this system (builder pattern)
    ///
    /// # Example
    /// ```
    /// use star_system::{SystemMetadata, GenerationMethod, SystemArchitecture};
    ///
    /// let meta = SystemMetadata::new_random(
    ///     GenerationMethod::Manual,
    ///     SystemArchitecture::Mixed,
    /// ).with_name("Sol");
    ///
    /// assert_eq!(meta.name, Some("Sol".to_string()));
    /// assert_eq!(meta.display_name(), "Sol");
    /// ```
    pub fn with_name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    /// Set binary configuration for this system (builder pattern)
    ///
    /// # Example
    /// ```
    /// use star_system::{
    ///     SystemMetadata,
    ///     GenerationMethod,
    ///     SystemArchitecture,
    ///     BinaryConfiguration,
    ///     BinaryOrbitType,
    ///     OrbitalParameters
    /// };
    /// use units::{Length, Time};
    ///
    /// let orbit = OrbitalParameters {
    ///     semi_major_axis: Length::from_au(10.0),
    ///     eccentricity: 0.1,
    ///     inclination: 0.0,
    ///     longitude_of_ascending_node: 0.0,
    ///     argument_of_periapsis: 0.0,
    ///     mean_anomaly: 0.0,
    ///     period: Time::from_years(10.0),
    /// };
    ///
    /// let binary_config = BinaryConfiguration {
    ///     orbit_type: BinaryOrbitType::STypePrimary,
    ///     orbital_params: orbit,
    /// };
    ///
    /// let meta = SystemMetadata::new_random(
    ///     GenerationMethod::Statistical,
    ///     SystemArchitecture::Mixed,
    /// ).with_binary_config(binary_config);
    ///
    /// assert!(meta.binary_config.is_some());
    /// ```
    pub fn with_binary_config(mut self, config: crate::binary::BinaryConfiguration) -> Self {
        self.binary_config = Some(config);
        self
    }
}
