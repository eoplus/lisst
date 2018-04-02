
# lisst (in development...)

An R package to read, manipulate and visualize data from the Laser In-Situ Scattering and Transmissometry (LISST) instruments, Sequoia Scientific, Inc.

### Provided functionality:
- Absolute particle VSF retrieval from binary file;
- Particle number concentraton from ppm Volume;
- PSD model fitting;
- Automatic data quality control;
- Subseting and ploting utilities.

### Instrument models suported:
- LISST-100(X)
- LISST-200X

### Note:
As the curent verion, is not yet possible to invert directly the raw data from LISST into particle size distribution (PSD) since the courtesy proprietary code provided by the vendor for the inversion is a source code in MATLAB P-code format, that can only be read by MATLAB (MathWorks, Inc).
