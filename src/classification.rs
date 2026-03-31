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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hr_position_serde_roundtrip() {
        let pos = HrPosition {
            temperature_k: 5778.0,
            absolute_magnitude: 4.83,
        };
        let json = serde_json::to_string(&pos).unwrap();
        let back: HrPosition = serde_json::from_str(&json).unwrap();
        assert!((back.temperature_k - 5778.0).abs() < f64::EPSILON);
        assert!((back.absolute_magnitude - 4.83).abs() < f64::EPSILON);
    }
}
