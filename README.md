# SiteResponseToolKit - SRTK
The Site Characterisation and Seismic Response Analysis Toolkit for Python

Current features:

- Site database and site building tools
- Parsing site model of arbitrary format (standard is csv) using generic I/O ASCII library
- Compute travel-time average velocity for variable depth (default is Vs30)
- Compute site class (presently only EC8, no special classes yet)
- Compute Quarter-Wavelength average parameters (velocity and density) and amplification
- Compute Kappa0 for arbitrary depth from Qs profile (default is whole profile)
- Compute SH-wave Transfer Function (elastic/anelastic) for arbitrary angle of incidence
- Compute resonance frequencies (from ShTf) and corresponding amplitudes

To do:
- Linear equivalent soil response
- Methods to adjust for reference Vs and Kappa
- Response spectral amplification using RVT
- Waveform convolution and basic signal processing methods
- Soil profile randomisation
- Implement Xml database file
