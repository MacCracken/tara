//! Hertzsprung-Russell diagram classification and stellar taxonomy.

use serde::{Deserialize, Serialize};

/// Placeholder for HR diagram position.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[non_exhaustive]
pub struct HrPosition {
    /// Effective temperature in Kelvin.
    pub temperature_k: f64,
    /// Absolute magnitude.
    pub absolute_magnitude: f64,
}
