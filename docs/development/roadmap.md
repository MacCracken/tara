# Tara Roadmap

## Review Polish (pre-publish)

- [ ] Add `#[non_exhaustive]` to `ZamsProperties` and `MsProperties` structs in `sse.rs`
- [ ] Add `#[inline]` to `sse::coeff` — hot-path Horner evaluation used in all SSE calculations
- [ ] Add SSE benchmarks to `benches/benchmarks.rs` — ZAMS, TMS, lifetime, and MS interpolation functions have zero benchmark coverage
- [ ] Add convergence check to Newton iteration in `soorat.rs` (`temperature_to_bv`) — currently 20 fixed iterations with no early exit on convergence
- [ ] Document hardcoded phase multipliers (1.1, 1.3, 1.4, 1.5, 1.6) in `evolutionary_phase` — will drift if `main_sequence_lifetime` changes
- [ ] Clarify `class_temp_range` return order — returns (max, min) which is inverted relative to naming convention
- [ ] Add SSE integration tests to `tests/integration.rs` — solar ZAMS/TMS validation, metallicity consistency

## Backlog

- SSE hook terms: α_L, β_L, η_L, α_R, β_R, γ_R coefficients for full MS interpolation with hook feature
- SSE post-main-sequence phases (Hurley et al. 2000): subgiant, RGB, HB, AGB
- SSE low-mass R_TMS piecewise formula (a52–a56 branch for M < a62)

## Future

- Binary/multiple star systems
- Stellar population synthesis
- Isochrone generation
- Interstellar extinction / reddening corrections
