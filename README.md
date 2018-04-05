# lisst (in development...)

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/lisst)](https://cran.r-project.org/package=lisst)

An R package to read, manipulate and visualize data from the Laser In-Situ Scattering and Transmissometry (LISST) instruments, Sequoia Scientific, Inc.

```
# Install from Github:
devtools::install_github("AlexCast/lisst")
```

### Provided functionality:
- Absolute particle VSF retrieval from binary file;
- Particle number concentraton from ppm Volume;
- PSD model fitting; (in development...)
- Automatic data quality control; (in development...)
- Subseting and ploting utilities;
- Automatic units track and conversion with the units package;
- Automatic uncertainty track and propagation with the errors package. (in development...)

### Instrument models suported:
- LISST-100(X)
- LISST-200X (in development...)


### Note:
As the curent verion, is not yet possible to invert directly the raw binary data from LISST into particle size distribution (PSD) since the courtesy proprietary code provided by the vendor for the inversion is a source code in MATLAB P-code format, that can only be read by MATLAB (MathWorks, Inc).
