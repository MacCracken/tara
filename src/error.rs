//! Error types for the Tara stellar astrophysics engine.

/// Errors that can occur in Tara operations.
#[derive(Debug, thiserror::Error)]
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

    /// An I/O error occurred.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}

/// A convenience `Result` type for Tara operations.
pub type Result<T> = std::result::Result<T, TaraError>;
