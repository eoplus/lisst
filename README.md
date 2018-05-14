# lisst

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip) [![Build Status](https://travis-ci.org/AlexCast/lisst.svg?branch=master)](https://travis-ci.org/AlexCast/lisst) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/AlexCast/lisst?branch=master&svg=true)](https://ci.appveyor.com/project/AlexCast/lisst) [![Coverage Status](https://img.shields.io/codecov/c/github/AlexCast/lisst/master.svg)](https://codecov.io/github/AlexCast/lisst?branch=master) [![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/lisst)](https://cran.r-project.org/package=lisst)

An R package to read, manipulate and visualize data from the Laser In-Situ Scattering and Transmissometry (LISST) instruments, Sequoia Scientific, Inc.

### Install from Github:
```
# install.packages(remotes)
remotes::install_github("AlexCast/lisst")
```

### Provided functionality:
- Subseting by sample, depth or time;
- Particle VSF retrieval from binary file;
- Particle number concentraton (1/L/µm) from ppm volume;
- PSD model fitting;
- Hovmöller diagram plotting per sample, depth or time;
- Diverse plot utilities; (in development...)
- Automatic data quality control; (in development...)
- Automatic units conversion with the [units](https://github.com/r-quantities/units) package;

### Instrument models suported:
- LISST-100(X);
- LISST-200X.

### Note:
As the curent verion, is not yet possible to invert directly the raw binary data from LISST into particle size distribution (PSD) since the courtesy proprietary code provided by the manufacturer for the inversion is a source code in MATLAB P-code format, that can only be read by MATLAB (MathWorks, Inc).
