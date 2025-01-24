PyEMILI
=======

**PyEMILI** is a Python-based tool primarily designed to assist in the automated line identification of spectra, especially for emission-line objects such as planetary nebulae (PNe), H II regions, and Herbig-Haro (HH) objects, etc. Compared to manual line identifications in the literature, PyEMILI achieves approximately 90% agreement. The current version of PyEMILI offers the following features.

Key Features
------------

- **Automated Emission Line Identification**: Identify emission lines in various astronomical objects.
- **Automated Absorption Line Identification**: Identify absorption lines in stellar spectra.
- **Automatic 1D Spectrum Processing**:
  
  - **Continuum Fitting**: Automatically fits the continuum of the 1D spectrum.
  - **Line Detection**: Finds spectral lines in the spectrum. An interactive interface is available for manually reviewing, adding, or revising spectral lines.
  
.. note::
   PyEMILI uses the `Atomic Line List v3.00b4 <https://www.pa.uky.edu/~peter/newpage/index.html>`_ developed by Peter A. M. Van Hoof, covering wavelengths from **1000–20000 Å**. **Molecular transitions are not included in the current version**.

Requirements
------------

To run PyEMILI, you will need:

- **Python3** compiler
- Python libraries: `numpy`, `scipy`, `matplotlib`, `astropy`, `pandas`, and `numba`

Library Compatibility
~~~~~~~~~~~~~~~~~~~~~
Ensure that the versions of `numba` and `numpy` are compatible to avoid potential import errors, such as `ImportError: Numba needs NumPy 1.** or less`. It is recommended to use the latest versions of these libraries to avoid compatibility issues.

Installation
------------

You can install PyEMILI by downloading the source code:

Step 1: Download the code file from Releases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Step 2: Compile with pip
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In the root directory which contains the setup.py file, run:

.. code-block:: bash

   pip install .

Another Option
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Directly install using pip:

.. code-block:: bash

   pip install pyemili

Usage
-----

Detailed usage instructions and examples can be found in the `Manual <https://github.com/LuShenJ/PyEMILI/tree/main/manual>`_, including descriptions of output files and additional parameter configurations. Examples in `test <https://github.com/LuShenJ/PyEMILI/tree/main/test>`_ directory include results for three samples, IC 418, Hf 2-2, and J0608.

Basic Example
~~~~~~~~~~~~~
In general, you need a text file containing **at least the wavelengths and fluxes of the unidentified lines**. If you need to use the automatic line search function in PyEMILI, and generate a line list to be identified, refer to `Line_finding <https://github.com/LuShenJ/PyEMILI/blob/main/manual/Line_finding.md>`_.

.. code-block:: python

   from pyemili.Lines import Line_list
   import numpy as np

   # Load the line list
   hf22 = np.loadtxt('Hf2-2_linelist.txt',skiprows=1)

   hf22_out = Line_list(wavelength=hf22[:,0], # array of wavelengths
                        wavelength_error=10,  # set a uniform wavelength uncertainty of 10 km/s
                        flux=hf22[:,1])       # array of fluxes 

   # run the line identification process, set the output filename as 'Hf2-2', 
   # and use the preset nebular elemental abundance
   hf22_out.identify('Hf2-2',abun_type='nebula') 

Output Files
~~~~~~~~~~~~
After running `pyemili.Lines.Line_list.identify()`, two files ending with **'.dat'** and **'.out'** will be generated in the directory. The '.out' file contains complete candidate IDs of each input observed line, and '.dat' file contains primarily the A ranking candidate IDs for each line. More information can be found `here <https://github.com/LuShenJ/PyEMILI/blob/main/manual/Intro.md>`_.

Troubleshooting
---------------

If you encounter any issues with PyEMILI, e.g., installation problems, usage problems, or questions on line identification accuracy, please feel free to reach out to us.  

Contact:

- Email: `zjtu@bao.ac.cn <mailto:zjtu@bao.ac.cn>`_ or `fangx@nao.cas.cn <mailto:fangx@nao.cas.cn>`_

Alternatively, please open an issue on the GitHub repository, where we’ll be happy to assist.
