# lisst (in development...)

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/lisst)](https://cran.r-project.org/package=lisst)

An R package to read, manipulate and visualize data from the Laser In-Situ Scattering and Transmissometry (LISST) instruments, Sequoia Scientific, Inc.

### Install from Github:
```
# install.packages(devtools)
devtools::install_github(paste("r-quantities", c("units", "errors", "quantities"), sep="/"))
devtools::install_github("AlexCast/lisst")
```

### Provided functionality:
- Particle VSF retrieval from binary file;
- Particle number concentraton (1/L/µm) from ppm volume;
- PSD model fitting; 
- Diverse plot utilities; (in development...)
- Automatic data quality control; (in development...)
- Subseting by sample index, depth or time;
- Automatic units conversion and uncertainty track with the [quantities](https://github.com/r-quantities/quantities) package;

### Instrument models suported:
- LISST-100(X)
- LISST-200X (in development...)

### Note:
As the curent verion, is not yet possible to invert directly the raw binary data from LISST into particle size distribution (PSD) since the courtesy proprietary code provided by the manufacturer for the inversion is a source code in MATLAB P-code format, that can only be read by MATLAB (MathWorks, Inc).
