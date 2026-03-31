# Tara Roadmap

## 1.1.0 — Complete Single-Star Evolution

### SSE Post-MS Refinement

- Full `lgbtf` implementation using GB parameters for precise RGB luminosity (currently simplified)
- Core Helium Burning (HB) luminosity and radius: L_HeI, L_ZAHB, R_min prescriptions
- Early-AGB and TP-AGB luminosity/radius (thermal pulses, third dredge-up)
- SSE naked helium star evolution (kw=7,8,9)
- Detailed core mass growth on RGB/AGB via core mass-luminosity relation

### Stellar Winds & Mass Loss

- Reimers (1975) RGB wind prescription
- Vassiliadis & Wood (1993) / Bloecker (1995) AGB mass loss
- Vink et al. (2001) OB star line-driven winds
- Nugis & Lamers (2000) Wolf-Rayet winds
- de Jager et al. (1988) general fallback prescription

### Remnant Formation

- Fryer et al. (2012) rapid/delayed core-collapse remnant mass prescriptions
- Pair-instability and pulsational pair-instability supernova mass gaps
- Cummings et al. (2018) updated WD initial-final mass relation
- Metallicity-dependent remnant type boundaries

## 1.2.0 — Binary Evolution & Population Tools

### Binary Stellar Evolution (BSE)

- Hurley, Tout & Pols (2002) BSE fitting formulae
- Roche lobe overflow (RLOF) and mass transfer stability criteria
- Mass transfer rates (thermal / nuclear / dynamical timescale)
- Orbital evolution under mass loss and mass transfer
- Tidal circularization and synchronization (Zahn 1977)
- Magnetic braking (Rappaport et al. 1983)

### Common Envelope Evolution

- Alpha-lambda formalism (Webbink 1984, de Kool 1990)
- Binding energy parameter λ (tabulated or from structure models)
- Energy-balance with enthalpy option (Ivanova et al. 2013)

### Supernova Natal Kicks

- Maxwellian kick distribution (Hobbs et al. 2005, σ = 265 km/s)
- Reduced kicks for electron-capture supernovae
- Fallback-modulated kicks for black holes (Fryer et al. 2012)

### Population Synthesis Building Blocks

- Initial mass function sampling: Kroupa (2001), Chabrier (2003), Salpeter (1955)
- Isochrone generation from evolutionary tracks
- Interstellar extinction: Cardelli, Clayton & Mathis (1989), Fitzpatrick (1999)

## 1.3.0 — Observational Interface

### Compact Object Properties

- White dwarf mass-radius relation (Hamada & Salpeter 1961 / Chandrasekhar theory)
- White dwarf cooling curves (Mestel theory, Fontaine et al. 2001)
- Neutron star mass-radius from TOV equation with selectable EOS
- Black hole Kerr parameter from angular momentum conservation

### Synthetic Photometry

- Filter response curves: Johnson-Cousins UBVRI, 2MASS JHK, SDSS ugriz, Gaia G/BP/RP
- SED convolution with filter response for synthetic magnitudes
- Photometric system conversions

### Chemical Evolution

- AGB yields (Karakas & Lattanzio 2014)
- Massive star yields (Nomoto et al. 2013, Limongi & Chieffi 2018)
- Type Ia yields (Iwamoto et al. 1999)

## Future

- Stellar rotation and rotational mixing (Heger et al. 2000)
- Pre-main-sequence evolution (Hayashi tracks, Kelvin-Helmholtz contraction)
- Pulsation and variability (Cepheid/RR Lyrae instability strip, period-luminosity)
- Gravitational wave progenitor modeling (chirp mass, Peters 1964 inspiral)
- Circumstellar envelopes and AGB dust formation
- Opacity tables: OPAL (Iglesias & Rogers 1996), Ferguson et al. (2005)
- Equation of state: electron degeneracy, Coulomb corrections, NS EOS
