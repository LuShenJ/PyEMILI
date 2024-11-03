
## Introduction

PyEMILI is a python tool mainly used to assist in the line identifications of the spectrum, especially for the planetary nebulae (PNe). Its accuracy could achieve about 95% comparing with strict manual line identifications. Here are some features in the current version of PyEMILI:

* Identify the emission lines in many emission-line object, e.g., PNe, Wolf-Rayet star, and supernova remanent.
* Identify the absorption lines in the spectra of stars.
* Automaticaly read the 1D spectrum.
* Automaticaly fit the continuum.
* Automaticaly correct the radial velocity of the spectrum.
* Automaticaly find the spectral lines. (An interactive interface is available to manually revise or add spectral lines)

_**NOTE**_: The atomic transition database that PyEMILI uses is from [Atomic Line List v3.00b4](https://www.pa.uky.edu/~peter/newpage/index.html) developed by Peter A. M. Van Hoof. Wavelength coverage is **1000-20000** angstroms. **The molecular transitions is not included in current version.**

## Requirements

PyEMILI needs _Python 3_ compiler and `numpy`, `scipy`, `matplotlib`, `astropy`, `pandas`, and `numba` libraries.

## Installation

You could download the source code by:

```
git clone https://github.com/LuShenJ/PyEMILI.git
```

or just download the .ZIP file in the top of this repository.
Unzip the file and in the root directory, run this in the terminal:

```
python setup.py install
```

_**NOTE**_: The version of `numba` and `numpy` should be compatible with each other. Otherwise, you may get `ImportError: Numba needs NumPy 1.** or less`. Generally I recommend keeping the latest version.

## Manual

Usage of the code and examples are [here](./manual), including the description of the output file.

`/test/` directory includes the results of 3 samples (IC418, Hf 2-2, J0608) in paper.

## Troubleshooting

Any problems about the code and bugs or issues, please email <zjtu@bao.ac.cn> or <fangx@nao.cas.cn>.
