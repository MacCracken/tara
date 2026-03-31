//! Optional logging initialization using `tracing-subscriber`.
//!
//! Requires the `logging` feature.

/// Try to initialise the global tracing subscriber.
///
/// # Errors
///
/// Returns an error if a global subscriber is already set.
pub fn try_init() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    tracing_subscriber::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
        .try_init()
}
