# SageMath Library for Cryptographic Functions

This project provides a SageMath library for working with cryptographic functions, with a particular focus on quadratic APN (Almost Perfect Nonlinear) functions.

The library aims to make research workflows involving APN families, representation conversions, and family membership testing easier to reproduce and extend.


## Features
- Conversion between major representations of APN functions, univariate polynomial to / from  
    - Truth table
    - Quadratic Matrix
    - Sequence
    - Bivariate
    - Trivariate
- Implementations of known infinite families of quadratic APN functions in univariate form
- Membership testing for identifying whether a function belongs to a known family


## Requirements
- SageMath 10.8 +
- Python 3.x


## Installation
In sage terminal write
```
pip install git+https://github.com/hoeeg/SageMath-library-for-cryptographic-functions.git
```

Or in regular terminal write
```
sage --pip install git+https://github.com/hoeeg/SageMath-library-for-cryptographic-functions.git
```


## How to run


## Don't have Sage?
Those who are impatient may use prebuilt Sage available online from any of

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/sagemath/sage-binder-env/master
) &nbsp; [![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-Ready--to--Code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/sagemath/sage/tree/master
) &nbsp; [![Open in GitHub Codespaces](https://img.shields.io/badge/Open_in_GitHub_Codespaces-black?logo=github)](https://codespaces.new/sagemath/sage/tree/master)

without local installation. Otherwise read on.

The [Sage Installation Guide](https://doc.sagemath.org/html/en/installation/index.html)
provides a decision tree that guides you to the type of installation
that will work best for you. This includes building from source,
obtaining Sage from a package manager, using a container image, or using
Sage in the cloud.

More about Sage or how to download it, can be found in its [Documentation](https://doc.sagemath.org/html/en/index.html).


## Layout
```
├── library/
│   ├── convert.py        Representation conversions
│   ├── families.py       Known APN family constructions
│   ├── helpers.py        Finite field utility functions
│   └── membership.py     Family identification algorithms
├── README.md             This file
├── setup.py              Top-level configure script
└── LICENCE               Copyright information
```
