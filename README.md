<img alt="OQ-SRTK - The Seismic Site Response Toolkit" class="right" style="width: 60%" src="https://raw.githubusercontent.com/klunk386/SeismicSiteTool/master/Logo/OQ-SRTK-Logo.png" />

[![AGPLv3](https://www.gnu.org/graphics/agplv3-88x31.png)](https://www.gnu.org/licenses/agpl.html)

# OQ-SRTK - OpenQuake Site Response Toolkit

A Site Characterisation and Seismic Response Analysis Toolkit for Python

Current features:

  * Site database and site building tools
  * Parsing site model of arbitrary format (standard is csv) using generic I/O ASCII library
  * Compute travel-time average velocity for variable depth (default is Vs30)
  * Compute site class (presently only EC8, no special classes yet)
  * Compute Quarter-Wavelength average parameters (velocity and density) and amplification
  * Compute Kappa0 for arbitrary depth from Qs profile (default is whole profile)
  * Compute SH-wave Transfer Function (elastic/anelastic) for arbitrary angle of incidence
  * Compute resonance frequencies and corresponding amplitudes
  * Basic signal processing

To do:

  * Linear equivalent soil response
  * Methods to adjust for reference Vs and Kappa
  * Response spectral amplification using RVT
  * Waveform convolution and basic signal processing methods
  * Soil profile randomisation
  * Implement Xml database file

### Dependencies

OQ-SRTK requires the following dependencies:

  * [NumPy/Scipy](http://www.scipy.org/)
  * [Matplotlib](http://matplotlib.org/)

### License

Copyright (c) 2017 GEM Foundation

OQ-SRTK is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

You should have received a copy of the GNU Affero General Public License with this download. If not, see <http://www.gnu.org/licenses/>

### Disclaimer

The software provided herein is released as a prototype implementation on behalf of scientists and engineers working within the GEM Foundation (Global Earthquake Model).

It is distributed for the purpose of open collaboration and in the hope that it will be useful to the scientific, engineering, disaster risk and software design communities.

The software is NOT distributed as part of GEM’s OpenQuake suite (http://www.globalquakemodel.org/openquake) and must be considered as a separate entity. The software provided herein is designed and implemented by scientific staff. It is not developed to the design standards, nor subject to same level of critical review by professional software developers, as GEM’s OpenQuake software suite.

Feedback and contribution to the software is welcome, and can be directed to the hazard scientific staff of the GEM Model Facility (hazard@globalquakemodel.org).

The Site Response Toolkit is therefore distributed WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

The GEM Foundation, and the authors of the software, assume no liability for use of the software.
