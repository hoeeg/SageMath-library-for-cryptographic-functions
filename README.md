# SageMath Library for Cryptographic Functions

This project provides a SageMath library for working with cryptographic functions, with a particular focus on quadratic APN (Almost Perfect Nonlinear) functions.

The library aims to make research workflows involving APN families, representation conversions, and family membership testing easier to reproduce and extend.


## Features
- Conversion between principal representations of vectorial Boolean functions, univariate polynomial to / from
    - Truth table
    - Quadratic Matrix
    - Sequence
    - Bivariate
    - Trivariate
- Implementations of known infinite families of quadratic APN functions in univariate form
- Membership testing for identifying whether a function belongs to a known family


## Requirements
- SageMath 10.8 +
- Python 3.9 +


## Installation
### From PyPI (recommended)

In a Sage terminal write
```
pip install sagemath-cryptographic-functions-library
```

Or from a regular terminal write
```
sage --pip install sagemath-cryptographic-functions-library
```

### From GitHub
In a Sage terminal write
```
pip install git+https://github.com/hoeeg/SageMath_library_for_cryptographic_functions.git
```
Or in regular terminal write
```
sage --pip install git+https://github.com/hoeeg/SageMath_library_for_cryptographic_functions.git
```


## Don't have Sage?
You can try Sage online without installing anything locally from any of

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/sagemath/sage-binder-env/master
) &nbsp; [![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-Ready--to--Code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/sagemath/sage/tree/master
) &nbsp; [![Open in GitHub Codespaces](https://img.shields.io/badge/Open_in_GitHub_Codespaces-black?logo=github)](https://codespaces.new/sagemath/sage/tree/master)


For a local install, the [Sage Installation Guide](https://doc.sagemath.org/html/en/installation/index.html) walks through the available options.

More about Sage or how to download it, can be found in its [Documentation](https://doc.sagemath.org/html/en/index.html).


## Layout
```
├── sagemath_cryptographic_functions_library/
│   ├── __init__.py
│   ├── conversions.py    Representation conversions
│   ├── families.py       Known APN family constructions
│   ├── helpers.py        Finite field utility functions
│   └── membership.py     Family identification algorithms
├── README.md             This file
├── pyproject.toml        Build configuration
└── LICENCE               Copyright information
```
