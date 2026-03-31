//! Error types for the Tara stellar astrophysics engine.

use serde::{Deserialize, Serialize};

/// Errors that can occur in Tara operations.
#[derive(Debug, thiserror::Error, Serialize, Deserialize)]
#[non_exhaustive]
pub enum TaraError {
    /// An invalid parameter was provided.
    #[error("invalid parameter: {0}")]
    InvalidParameter(String),

    /// A mathematical computation failed.
    #[error("math error: {0}")]
    MathError(String),

    /// A spectral analysis operation failed.
    #[error("spectral error: {0}")]
    SpectralError(String),

    /// A model computation failed.
    #[error("model error: {0}")]
    ModelError(String),

    /// An I/O error occurred (stored as description for serde compatibility).
    #[error("I/O error: {0}")]
    Io(String),
}

impl From<std::io::Error> for TaraError {
    fn from(err: std::io::Error) -> Self {
        Self::Io(err.to_string())
    }
}

/// A convenience `Result` type for Tara operations.
pub type Result<T> = std::result::Result<T, TaraError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tara_error_serde_roundtrip() {
        let errors = [
            TaraError::InvalidParameter("bad mass".into()),
            TaraError::MathError("overflow".into()),
            TaraError::SpectralError("no lines".into()),
            TaraError::ModelError("diverged".into()),
            TaraError::Io("file not found".into()),
        ];
        for err in &errors {
            let json = serde_json::to_string(err).unwrap();
            let back: TaraError = serde_json::from_str(&json).unwrap();
            assert_eq!(err.to_string(), back.to_string());
        }
    }

    #[test]
    fn io_error_converts() {
        let io_err = std::io::Error::new(std::io::ErrorKind::NotFound, "missing");
        let tara_err: TaraError = io_err.into();
        assert!(matches!(tara_err, TaraError::Io(_)));
    }
}
