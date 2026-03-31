//! Optional logging initialization using `tracing-subscriber`.
//!
//! Requires the `logging` feature. Configures structured logging via the
//! `RUST_LOG` environment variable (e.g. `RUST_LOG=tara=debug`).
//!
//! # Example
//!
//! ```ignore
//! tara::logging::try_init().expect("failed to initialize logging");
//! // Now all tara tracing events are emitted to stderr.
//! ```

/// Try to initialise the global tracing subscriber with env-filter support.
///
/// Reads the `RUST_LOG` environment variable for filtering directives.
/// Defaults to `warn` level if `RUST_LOG` is not set.
///
/// # Errors
///
/// Returns an error if a global subscriber is already set.
pub fn try_init() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("warn")),
        )
        .with_target(true)
        .with_thread_ids(false)
        .with_file(false)
        .with_line_number(false)
        .try_init()
}

/// Try to initialise a compact logging subscriber suited for development.
///
/// Includes source file and line numbers. Defaults to `debug` if `RUST_LOG`
/// is not set.
///
/// # Errors
///
/// Returns an error if a global subscriber is already set.
pub fn try_init_dev() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("debug")),
        )
        .with_target(true)
        .with_thread_ids(false)
        .with_file(true)
        .with_line_number(true)
        .compact()
        .try_init()
}
