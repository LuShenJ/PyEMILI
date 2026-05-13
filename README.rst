PyEMILI
=======

**PyEMILI** (`Tu, Fang et al. 2025 <https://iopscience.iop.org/article/10.3847/1538-4365/adae00>`_) is a Python-based tool for automated spectral-line identification in astronomical spectra. It is designed primarily for emission-line objects, including planetary nebulae (PNe), H II regions, Herbig-Haro (HH) objects, and related sources, while also providing functionality for absorption-line identification and one-dimensional spectrum processing.

PyEMILI follows the general idea of EMILI, but implements the workflow in Python with updated atomic-data handling, automated ranking of candidate IDs, multiplet checks, and optional iterative refinement of ionization and velocity corrections. In comparison with manual line identifications in the literature, PyEMILI reaches approximately 90% agreement for tested emission-line spectra.

Key Features
------------

- **Automated emission-line identification**: identify observed emission lines and rank possible candidate IDs.
- **Automated absorption-line identification**: support line identification in stellar spectra.
- **Automatic 1D spectrum processing**:

  - **Continuum fitting**: fit the continuum of a one-dimensional spectrum.
  - **Line detection**: detect spectral lines and provide an interactive interface for manual inspection, addition, and revision.

- **Candidate ranking**: evaluate candidate IDs using wavelength agreement, predicted flux, ionization balance, and multiplet consistency.
- **Iterative identification**: optionally refine the ionization correction factors and velocity corrections using robust first-pass identifications.
- **Recombination-line fitting**: estimate electron temperature, electron density, and ionic abundance from recombination lines.

.. note::

   PyEMILI uses the `Atomic Line List v3.00b4 <https://www.pa.uky.edu/~peter/newpage/index.html>`_ developed by Peter A. M. van Hoof. The current atomic transition database covers wavelengths from **1000 to 20000 Angstrom**. Molecular transitions are not included in the current version.

Requirements
------------

PyEMILI requires Python 3 and the following Python packages:

- ``numpy``
- ``scipy``
- ``matplotlib``
- ``astropy``
- ``pandas``
- ``numba``
- ``tqdm``

Package Compatibility
~~~~~~~~~~~~~~~~~~~~~

Please make sure that the installed versions of ``numba`` and ``numpy`` are mutually compatible. If you encounter an error such as ``ImportError: Numba needs NumPy 1.** or less``, install a compatible pair of ``numpy`` and ``numba`` versions in a clean environment.

Installation
------------

PyEMILI is under active development. For the most up-to-date version, we recommend installing from the GitHub repository.

Install from GitHub
~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone https://github.com/LuShenJ/PyEMILI.git
   cd PyEMILI
   pip install .

Install from PyPI
~~~~~~~~~~~~~~~~~

A PyPI version is also available:

.. code-block:: bash

   pip install pyemili

If you plan to use the newest development version, please use the GitHub installation method above.

Usage
-----

Detailed usage instructions and examples are available in the `manual <https://github.com/LuShenJ/PyEMILI/tree/main/manual>`_ directory. Example data and results are provided in the `test <https://github.com/LuShenJ/PyEMILI/tree/main/test>`_ directory, including IC 418, Hf 2-2, and J0608.

Basic Example
~~~~~~~~~~~~~

A basic PyEMILI run requires a line list containing at least the observed wavelengths and fluxes of unidentified lines. If you need to generate such a line list automatically from a one-dimensional spectrum, see `Line_finding <https://github.com/LuShenJ/PyEMILI/blob/main/manual/Line_finding.md>`_.

.. code-block:: python

   import numpy as np
   from pyemili.Lines import Line_list

   # Load an observed line list.
   # The first column is wavelength in Angstrom, and the second column is flux.
   hf22 = np.loadtxt("Hf2-2_linelist.txt", skiprows=1)

   hf22_output = Line_list(
       wavelength=hf22[:, 0],
       wavelength_error=10,  # uniform 1-sigma wavelength uncertainty in km/s
       flux=hf22[:, 1],
   )

   # Run the line-identification process and write outputs with prefix "Hf2-2".
   # Here we use the built-in nebular abundance table.
   hf22_output.identify("Hf2-2", abun_type="nebula")

Input Format
~~~~~~~~~~~~

The main input parameters of ``pyemili.Lines.Line_list`` are:

- ``wavelength``: observed wavelengths in Angstrom.
- ``wavelength_error``: wavelength uncertainty. This can be provided as a single velocity uncertainty in km/s, a one-dimensional array of wavelength uncertainties in Angstrom, or a two-dimensional array giving asymmetric blue/red uncertainties for each line.
- ``flux``: observed line fluxes. Fluxes are commonly normalized to H beta = 1.
- ``ral_vel``: optional radial-velocity correction in km/s.
- ``flux_error``: optional flux uncertainties.
- ``snr``: optional signal-to-noise ratios, used in the output table.
- ``fwhm``: optional full width at half maximum values, used in the output table.

Important Identification Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The main parameters of ``Line_list.identify()`` include:

- ``filename``: prefix of the generated output files.
- ``icf``: ionization correction factors for the five energy bins. If not provided, PyEMILI uses default values or estimates them iteratively.
- ``v_cor``: velocity corrections for the five energy bins. If not provided, PyEMILI uses default values or estimates them iteratively.
- ``sigma``: search window in units of the wavelength uncertainty. The default value is ``5``.
- ``Ne``: electron density in cm^-3. The default value is ``10000``.
- ``Te``: electron temperature in K. The default value is ``10000``.
- ``I``: instrumental resolution in km/s used for multiplet checks. The default value is ``10``.
- ``abun_type``: abundance table. Built-in options are ``"solar"`` and ``"nebula"``. A custom abundance-file path can also be supplied.
- ``deplete``: optional depletion or enhancement of selected ions, for example ``"Fe 10 1 3"``.
- ``col_cor``: dilution factor for collisionally excited lines. The default value is ``0.1``.
- ``iteration``: if ``True``, PyEMILI performs a first-pass identification and then re-identifies lines with updated ``icf`` and ``v_cor``.
- ``erc_list``: if ``True``, PyEMILI writes an additional recombination-line fitting list.
- ``match_list``: optional matched-line table from a previous identification.

Output Files
~~~~~~~~~~~~

After running ``pyemili.Lines.Line_list.identify()``, PyEMILI generates two main output files:

- ``<filename>.out``: complete candidate IDs for each input observed line.
- ``<filename>.dat``: a short output table containing primarily the best-ranked candidate IDs.

If ``erc_list=True`` is used, PyEMILI also writes:

- ``<filename>_erc.dat``: a recombination-line fitting list.

More information about the output format can be found in the `manual introduction <https://github.com/LuShenJ/PyEMILI/blob/main/manual/Intro.md>`_.

Citation
--------

If you use PyEMILI in your research, please cite:

- `Tu, Fang et al. 2025 <https://iopscience.iop.org/article/10.3847/1538-4365/adae00>`_

Please also cite the original atomic data sources and relevant references when appropriate.

Troubleshooting and Contact
---------------------------

If you encounter installation problems, usage issues, or questions about line-identification accuracy, please open an issue on the GitHub repository or contact the developers.

Contact:

- Email: `zjtu@bao.ac.cn <mailto:zjtu@bao.ac.cn>`_
