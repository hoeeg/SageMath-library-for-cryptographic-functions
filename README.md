# SageMath Library for Cryptographic Functions

A SageMath library for cryptographic functions, focusing on quadratic APN (Almost Perfect Nonlinear) functions.

The library provides tools for:
- converting between different representations of a function
- generating known infinite families of quadratic APN polynomials in univariate form
- deciding whether a function belongs to these families

## Features
- Conversion routines between univariate polynomial representation and other common representations of vectorial Boolean functions:
    - truth tables
    - quadratic matrices
    - sequences
    - bivariate forms
    - trivariate forms
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

## Usage
The package can either be imported by module
```
from sagemath_cryptographic_functions_library import conversions, families, membership
families.family_4(3)
...
```
or by importing functions directly
```
from sagemath_cryptographic_functions_library import family_4
family_4(3)
...
```

## Using Sage Online
You can try Sage online without installing anything locally from any of

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/sagemath/sage-binder-env/master
) &nbsp; [![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-Ready--to--Code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/sagemath/sage/tree/master
) 

For a local install, the [Sage Installation Guide](https://doc.sagemath.org/html/en/installation/index.html) walks through the available options.

More about Sage or how to download it, can be found in its [Documentation](https://doc.sagemath.org/html/en/index.html).

## Layout
```
├── sagemath_cryptographic_functions_library/
│   ├── __init__.py
│   ├── conversions.py    Conversions between representations
│   ├── families.py       Known APN family constructions
│   ├── helpers.py        Finite field utility functions
│   └── membership.py     Membership testing algorithms
├── README.md             This file
├── pyproject.toml        Build configuration
└── LICENCE               Copyright information
```