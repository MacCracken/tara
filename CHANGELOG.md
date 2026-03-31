# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Added

- Initial scaffold — core `Star` type, `SpectralClass` enum, `TaraError`, module stubs
- Cross-crate bridge module — falak (gravitational parameter, absolute magnitude), prakash (Wien's law, limb darkening), tanmatra (lifetime scaling, CNO fraction)
- Soorat integration module — HR diagram, evolution track, spectral profile, star field visualization structs (feature-gated `soorat-compat`)
- Logging module with tracing-subscriber (feature-gated `logging`)
- Feature gates: `optics`, `thermal`, `nuclear` for optional dependency bridges
- HR diagram classification stub
