"""
Recombination-Line Fitter
-----------------------------------------

Fit electron temperature (Te), electron density (Ne), and (optionally) ionic abundance
from recombination-line intensities using precomputed effective recombination coefficients.

Design (per user requirements):
- Class "Recom_Lines" loads the line list; all fitting options are passed to fit().
- Single public entry: .fit(...)
- Abundance is parameterized as A12 = 12 + log10(X/H^+); internally X/H^+ = 10**(A12-12).
- Normalization policy:
    * If fit_abundance is True  -> normalize to Hbeta ("Hbeta"). In this case, you NEED to input the flux also normalized to Hbeta.
    * If fit_abundance is False -> normalize within ion by its own max flux.


Inputs
------
A plain-text line list with columns:
    obs_wav  obs_flux  obs_fluxerr  ele  ion
e.g. 
    4085.11  6.70E-03   1.68E-03     O   II
whitespace-separated; lines starting with '#' are ignored.

You can simply obtain an input line list used for recombination line fitting by using pyemili.Lines.Line_list.identify(), 
and set erc_list=True. A file ends with 'erc.dat' will be generated, and can be strightly read by pyemili.recomb.Recom_Lines.

Dependencies: emcee, corner

"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import logging
import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
import emcee
import corner


# -----------------------
# Logging configuration
# -----------------------
logger = logging.getLogger("recom_lines")
if not logger.handlers:
    handler = logging.StreamHandler(sys.stdout)
    fmt = logging.Formatter("[%(levelname)s] %(message)s")
    handler.setFormatter(fmt)
    logger.addHandler(handler)
logger.setLevel(logging.INFO)


@dataclass
class FitResult:
    """Container for MCMC results.

    Attributes
    ----------
    flat_samples : np.ndarray
        Posterior samples after burn-in and thinning. Shape (Nsamples, Ndim).
    medians : Dict[str, float]
        Posterior medians for fitted parameters.
    intervals : Dict[str, Tuple[float, float]]
        16th and 84th percentiles for fitted parameters.
    labels : List[str]
        Parameter labels in the same order as flat_samples columns.
    fixed_params : Dict[str, float]
        Parameters held fixed (and their fixed values).
    figure_paths : Dict[str, str]
        File paths of saved figures (corner and optional data-model scatter).
    """

    flat_samples: np.ndarray
    medians: Dict[str, float]
    intervals: Dict[str, Tuple[float, float]]
    labels: List[str]
    fixed_params: Dict[str, float]
    figure_paths: Dict[str, str]


class Recom_Lines:
    """
    Recombination-line fitter; this class only loads the line list in __init__,
    and all fitting options are specified in .fit().

    Parameters
    ----------
    filename : str
        Path to the line list file. Must contain columns
        [obs_wav, obs_flux, obs_fluxerr, ele, ion].
    rootdir : Optional[str], default None
        Path to the directory containing pyemili/eff_reccoe. If None, assumed to be
        the parent directory of this file.

    Notes
    -----
    - Abundance parameterization is A12 = 12 + log10(X/H). Internally we convert to
      (X/H) via 10**(A12-12).
    - For H I, abundance is meaningless and `fit_abundance` is forced False.
    """

    EFF_IONS: Tuple[str, ...] = ("H I", "He I", "N II", "O II", "C II")

    def __init__(
        self,
        filename: str,
        rootdir: Optional[str] = None,
    ) -> None:
        # Store path for figure naming
        self.filename = filename

        # Paths
        if rootdir is None:
            rootdir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
        self.rootdir = rootdir


        with open(filename) as f:
            first_line = f.readline().strip()
        ncols = len(first_line.split())

        # set column names based on number of columns
        if ncols == 5:
            colnames = ["obs_wav", "obs_flux", "obs_fluxerr", "ele", "ion"]
        elif ncols == 6:
            colnames = ["obs_wav", "obs_flux", "obs_fluxerr", "ele", "ion", "lab_wav"]
        else:
            raise ValueError(f"Unexpected number of columns ({ncols}) in file: {filename}")

        # load file
        df = pd.read_table(
            filename,
            sep=r"\s+",
            names=colnames,
            comment="#",
        )
        # df = df[df.obs_flux != df.obs_fluxerr]  # drop degenerate rows
        df = df.sort_values("obs_flux").reset_index(drop=True)
        # normalize ion label like 'II]' -> 'II'
        df["ion"] = df["ion"].astype(str).str.replace(r"\]$", "", regex=True)
        self.df = df

        # Lazy-loaded coefficients
        self.coeTe: Dict[str, np.ndarray] = {}
        self.coeNe: Dict[str, np.ndarray] = {}
        self.Coes: Dict[str, np.ndarray] = {}
        self.Lines: Dict[str, np.ndarray] = {}

        # Runtime state set in fit()
        self.spec: Optional[str] = None
        self.normalize_to: Optional[str] = None
        self.vel: Optional[float] = None
        self.default_eleabuns: Optional[np.ndarray] = None
        self.data: Optional[pd.DataFrame] = None
        self.use_abundance_in_model: bool = False

    # -----------------------
    # Data & coefficient utils
    # -----------------------
    def _load_coefficients(self) -> None:
        """Load Te, Ne grids and emissivity tables from pyemili/eff_reccoe.

        Expects files:
            eff_Te.npy, eff_Ne.npy, and for each ion in EFF_IONS an array
            '{ION}_emiss_data.npy' and line list '{ION}_emiss_lines.dat'.
        """
        base = os.path.join(self.rootdir, "pyemili", "eff_reccoe")
        self.coeTe = np.load(os.path.join(base, "eff_Te.npy"), allow_pickle=True).item()
        self.coeNe = np.load(os.path.join(base, "eff_Ne.npy"), allow_pickle=True).item()

        self.Coes = {}
        self.Lines = {}
        for label in ("HI", "HeI", "NII", "OII", "CII"):
            self.Coes[label] = np.load(os.path.join(base, f"{label}_emiss_data.npy"), allow_pickle=True).astype(np.float64)
            self.Lines[label] = np.loadtxt(os.path.join(base, f"{label}_emiss_lines.dat"))

    def _build_interpolators(self) -> None:
        """Create RegularGridInterpolators for the selected ion and set Te/Ne bounds."""
        # maps like "O II" -> "OII"
        key = self.spec.replace(" ", "")
        Te_grid = np.log10(self.coeTe[key])
        Ne_grid = np.log10(self.coeNe[key])

        self.Te_grid = Te_grid
        self.Ne_grid = Ne_grid
        self.Te_min, self.Te_max = float(Te_grid.min()), float(Te_grid.max())
        self.Ne_min, self.Ne_max = float(Ne_grid.min()), float(Ne_grid.max())

        if self.Ne_log:
            if self.Ne_log < self.Ne_min or self.Ne_log > self.Ne_max:
                raise ValueError(f"Input log(Ne)={self.Ne_log:.2f} is outside the model grid range [{self.Ne_min:.2f}, {self.Ne_max:.2f}].")

        if self.Te_log:
            if self.Te_log < self.Te_min or self.Te_log > self.Te_max:
                raise ValueError(f"Input log(Te)={self.Te_log:.2f} is outside the model grid range [{self.Te_min:.2f}, {self.Te_max:.2f}].")
        # Full emissivity cube for this ion: [n_lines, n_Te, n_Ne]
        self.emiss_cube = self.Coes[key]

        # Interpolator over (line_index, logTe, logNe)
        idx = np.arange(self.emiss_cube.shape[0])
        self.emiss_interp = interpolate.RegularGridInterpolator(
            (idx, Te_grid, Ne_grid), self.emiss_cube, bounds_error=False, fill_value=None
        )

    def _build_reference_normalizer(self) -> None:
        """Prepare normalization to Hβ or to max-within-ion depending on setting."""
        # Hβ emissivity table
        H_Te = np.log10(self.coeTe["HI"])  # H I grids
        H_Ne = np.log10(self.coeNe["HI"])
        H_cube = self.Coes["HI"]
        Hbeta_idx = 21  # the H_beta index in database
        self.Hbeta_interp = interpolate.RegularGridInterpolator(
            (H_Te, H_Ne), H_cube[Hbeta_idx, :, :], bounds_error=False, fill_value=None
        )

    def _prepare_line_coeff_sums(self) -> None:
        """Sum emissivities for each observed feature over lab lines within a velocity window."""
        key = self.spec.replace(" ", "")
        lab_wav = self.Lines[key][:, 0]

        coeffs = []
        for _, row in self.data.iterrows():
            low = row.obs_wav * (1.0 - self.vel / self.c_kms)
            high = row.obs_wav * (1.0 + self.vel / self.c_kms)
            mask = (lab_wav >= low) & (lab_wav <= high)
            summed = self.Coes[key][mask].sum(axis=0)  # shape (n_Te, n_Ne)
            coeffs.append(summed)

        self.coe_sums = np.stack(coeffs, axis=0)  # shape (n_obs, n_Te, n_Ne)
        idx = np.arange(self.coe_sums.shape[0])
        self.coe_interp = interpolate.RegularGridInterpolator(
            (idx, self.Te_grid, self.Ne_grid), self.coe_sums, bounds_error=False, fill_value=None
        )

        # Observed vectors (with chosen normalization)
        obs = self.data["obs_flux"].to_numpy(dtype=float)
        err = self.data["obs_fluxerr"].to_numpy(dtype=float)

        if self.normalize_to == "within_ion_max":
            ref = obs.max()
            obs, err = obs / ref, err / ref
            self.obs_norm = obs
            self.err_norm = err
            self.use_abundance_in_model = False  # abundance cancels
        else:
            # Expect user to have already normalized to I(Hβ)=1. If not, the model
            # will still work but results interpret as relative to Hβ.
            self.obs_norm = obs
            self.err_norm = err
            # For metal ions, abundance is needed to link to Hβ emissivity
            self.use_abundance_in_model = (self.spec != "H I")


    # -----------------------
    # Probabilistic model
    # -----------------------
    def _model_vector(self, Te_log: float, Ne_log: float, A12: Optional[float]) -> np.ndarray:
        """Compute model vector for all observed lines at given parameters.

        Parameters
        ----------
        Te_log : float
            log10 electron temperature [K].
        Ne_log : float
            log10 electron density [cm^-3].
        A12 : Optional[float]
            Abundance 12+log10(X/H). If None or abundance is disabled, the model
            is not scaled by abundance.

        Returns
        -------
        np.ndarray
            Model intensities in the same normalization as self.obs_norm.
        """
        # Optional clipping for stability (for Hbeta-based)
        if self.normalize_to == "Hbeta" and self.spec != "H I":
            Te_log_Hbeta = float(np.clip(Te_log, np.log10(500.0), np.log10(30000.0)))
            Ne_log_Hbeta = max(Ne_log, 2.0)

        # Build indices for interpolating the summed emissivities of each observed line
        n = self.coe_sums.shape[0]
        idx = np.arange(n)
        pts = np.column_stack([idx, np.full(n, Te_log), np.full(n, Ne_log)])
        emiss = self.coe_interp(pts)  # shape (n,)

        # Normalization route
        if self.normalize_to == "within_ion_max":
            ref = emiss.max()
            model = emiss / ref
        else:
            Hbeta = float(self.Hbeta_interp([[Te_log_Hbeta, Ne_log_Hbeta]])[0])
            model = emiss / Hbeta


        # Apply abundance if needed
        if self.use_abundance_in_model and A12 is not None and self.spec != "H I":
            XH = 10.0 ** (A12 - 12.0)
            model = XH * model

        return model

    def _log_prior(self, Te_log: float, Ne_log: float, A12: Optional[float]) -> float:
        # Te/Ne inside grids
        if not (self.Te_min <= Te_log <= self.Te_max):
            return -np.inf
        if not (self.Ne_min <= Ne_log <= self.Ne_max):
            return -np.inf

        # Abundance prior: within ±2 dex around A12_init
        if self.fit_abundance and self.use_abundance_in_model and self.spec != "H I":
            assert self.A12_init is not None
            if not (self.A12_init - 2.0 <= A12 <= self.A12_init + 2.0):
                return -np.inf
        return 0.0

    def _log_prob(self, theta: np.ndarray) -> float:
        # Unpack according to which parameters are free
        i = 0
        Te_log = self.Te_log if not self.fit_Te else float(theta[i]); i += int(self.fit_Te)
        Ne_log = self.Ne_log if not self.fit_Ne else float(theta[i]); i += int(self.fit_Ne)
        A12 = None
        if self.fit_abundance and self.use_abundance_in_model and self.spec != "H I":
            A12 = float(theta[i]); i += 1

        lp = self._log_prior(Te_log, Ne_log, A12)
        if not np.isfinite(lp):
            return -np.inf

        model = self._model_vector(Te_log, Ne_log, A12)
        chi2 = np.sum((self.obs_norm - model) ** 2 / (self.err_norm ** 2))
        return lp - 0.5 * chi2

    # -----------------------
    # Public API
    # -----------------------
    def fit(
        self,
        ion: str,
        Te: Optional[float] = None,
        Ne: Optional[float] = None,
        fit_abundance: bool = True,
        A12: Optional[float] = None,
        eleabuns: Optional[np.ndarray | List[float]] = None,
        flux_threshold: float = 0.0,
        flux_normalize: float = 1.0,
        vel: float = 20.0,
        use_lab_wav: bool = False,
        nwalkers: int = 32,
        max_iter: int = 30000,
        progress: bool = True,
        seed: Optional[int] = None,
        save_prefix: Optional[str] = None,
        make_corner: bool = True,
        make_scatter: bool = False) -> FitResult:
        """Run MCMC and (optionally) make figures.

        Parameters
        ----------
        ion : str
            Ion to fit, e.g., "O II", "N II", "He I", "H I", "C II".
        Te, Ne : float or None
            Fix electron temperature (K) and/or density (cm^-3). If None, the parameter is sampled.
        fit_abundance : bool
            Whether to fit A12 = 12 + log10(X/H). If True, normalization is set to Hbeta; otherwise
            it is set to within-ion-max.
        A12 : float or None
            Initial abundance for the prior center when fitting abundance.
        eleabuns : array-like or None
            Ionic abundances (X/H) for [H I, He I, N II, O II, C II] to set the abundance prior center.
            Default values are [1.0, 8.5e-2, 6.8e-5, 4.9e-4, 2.7e-4].
        flux_threshold : float
            Drop observed features below this flux after any normalization.
        flux_normalize : float
            Divide observed flux and errors by this scalar before fitting.
        vel : float
            Velocity window in km/s to collect lab lines per observed feature.
        use_lab_wav : bool
            If True, use the last column (lab_wav) as the wavelength for matching; if False, use the first column (obs_wav). Default False.
        nwalkers : int
            Number of emcee walkers.
        max_iter : int
            Maximum sampler iterations; early-stops on convergence.
        progress : bool
            Show emcee progress bar.
        seed : int or None
            Random seed for initial positions.
        save_prefix : str or None
            Figure save prefix.
        make_corner : bool
            Whether to save the corner plot.
        make_scatter : bool
            Whether to save the data-model scatter plot at posterior medians.


        Returns
        -------
        FitResult
            Posterior summary and figure paths.
        """
        # ----- runtime configuration (moved from __init__) -----
        # 1) ion validation
        if not isinstance(ion, str):
            raise TypeError("Ion must be a string, e.g., 'O II'.")
        spec = " ".join(ion.split())
        if spec not in self.EFF_IONS:
            raise ValueError(f"Unknown ion '{ion}'. Valid: {self.EFF_IONS}")
        self.spec = spec

        # 2) basic run-time settings
        self.c_kms = 2.9979246e5
        self.vel = float(vel)

        # --- NEW: enforce normalization policy based on fit_abundance ---
        # For H I, abundance is meaningless -> force fit_abundance False
        fit_abundance = bool(fit_abundance) and (self.spec != "H I")
        self.fit_abundance = fit_abundance
        self.normalize_to = "Hbeta" if self.fit_abundance else "within_ion_max"

        # 3) default abundances for prior center
        if eleabuns is None:
            self.default_eleabuns = np.array([1.0, 8.5e-2, 6.8e-5, 4.9e-4, 2.7e-4], dtype=float)
        else:
            eleabuns = np.asarray(eleabuns, dtype=float)
            if eleabuns.shape != (5,):
                raise ValueError("eleabuns must be length-5 for [H I, He I, N II, O II, C II]")
            self.default_eleabuns = eleabuns

        # 4) lazy-load coefficients
        if not self.coeTe:
            self._load_coefficients()

        # 5) subset & (optionally) normalize the data for this spectrum
        df = self.df.copy()
        if flux_normalize != 1.0:
            df["obs_flux"] = df["obs_flux"] / flux_normalize
            df["obs_fluxerr"] = df["obs_fluxerr"] / flux_normalize
        ele, ion = self.spec.split()
        sub = df[(df["ele"] == ele) & (df["ion"] == ion)]
        sub = sub[sub["obs_flux"] >= flux_threshold].reset_index(drop=True)
        if len(sub) == 0:
            raise RuntimeError("No lines found for requested spectrum after thresholding.")
        self.data = sub

        # Optionally switch wavelength reference to laboratory wavelengths
        if use_lab_wav:
            if 'lab_wav' not in self.data.columns:
                raise KeyError(
                    f"The file does not contain a 'lab_wav' column, cannot switch to laboratory wavelength reference."
                )
            self.data = self.data.copy()
            self.data['obs_wav'] = self.data['lab_wav'].astype(float)
        # 6) flags and fixed values
        self.fit_Te = Te is None
        self.fit_Ne = Ne is None
        self.Te_log = None if self.fit_Te else float(np.log10(Te))
        self.Ne_log = None if self.fit_Ne else float(np.log10(Ne))

        if self.spec == "H I":
            self.A12_init = None
        else:
            ion_order = ["H I", "He I", "N II", "O II", "C II"]
            ix = ion_order.index(self.spec)
            XH_default = self.default_eleabuns[ix]
            A12_default = 12.0 + np.log10(XH_default)
            self.A12_init = A12_default if A12 is None else float(A12)

        # 7) per-spectrum prep
        self._build_interpolators()
        self._build_reference_normalizer()
        self._prepare_line_coeff_sums()

        # ----- sampling state -----
        rng = np.random.default_rng(seed)

        # Build initial position vector
        p0: List[float] = []
        labels: List[str] = []
        fixed: Dict[str, float] = {}

        if self.fit_Te:
            Te0 = 0.5 * (self.Te_min + self.Te_max)
            p0.append(Te0); labels.append(r"$\log T_e$")
        else:
            fixed[r"$\log T_e$"] = self.Te_log

        if self.fit_Ne:
            Ne0 = 0.5 * (self.Ne_min + self.Ne_max)
            p0.append(Ne0); labels.append(r"$\log N_e$")
        else:
            fixed[r"$\log N_e$"] = self.Ne_log

        if self.fit_abundance and self.use_abundance_in_model and self.spec != "H I":
            assert self.A12_init is not None
            p0.append(self.A12_init)
            ion_charge = {"I": "+", "II": "2+"}[self.spec.split()[1]]
            labels.append(fr"$12+\log(\mathrm{{{self.spec.split()[0]}}}^{{{ion_charge}}}/\mathrm{{H^+}})$")
        else:
            if self.spec != "H I":
                fixed[fr"$12+\log({self.spec.split()[0]}^{{+}}/H^+)$"] = (
                    self.A12_init if self.A12_init is not None else np.nan
                )

        p0 = np.asarray(p0, dtype=float)
        if p0.size == 0:
            raise RuntimeError("No free parameters left to fit. Enable at least one of Te/Ne/abundance.")

        pos = p0 + 1e-1 * rng.standard_normal(size=(nwalkers, p0.size))

        # Set up sampler
        sampler = emcee.EnsembleSampler(nwalkers, p0.size, self._log_prob)

        # Auto-stop after convergence
        max_n = int(max_iter)
        burnin_detected = 0
        old_tau = np.inf
        taus = []
        for _ in sampler.sample(pos, iterations=max_n, progress=progress):
            if burnin_detected == 0 and sampler.iteration % 100 == 0:
                tau = sampler.get_autocorr_time(tol=0)
                taus.append(np.mean(tau))
                converged = np.all(tau * 100 < sampler.iteration) and np.all(np.abs(old_tau - tau) / tau < 0.01)
                if converged:
                    burnin_detected = sampler.iteration
                    logger.info("Convergence detected at iteration %d", burnin_detected)
                old_tau = tau
            elif burnin_detected and sampler.iteration >= burnin_detected + 1000:
                break

        # Final thinning/burn-in
        tau = sampler.get_autocorr_time()
        burn = int(2 * np.max(tau))
        thin = max(1, int(0.5 * np.min(tau)))
        flat = sampler.get_chain(discard=burn, flat=True, thin=thin)

        # Summaries
        meds = np.percentile(flat, 50, axis=0)
        p16 = np.percentile(flat, 16, axis=0)
        p84 = np.percentile(flat, 84, axis=0)
        intervals = {labels[i]: (p16[i], p84[i]) for i in range(len(labels))}
        medians = {labels[i]: meds[i] for i in range(len(labels))}

        # Figures
        if save_prefix is None:
            base = os.path.splitext(os.path.basename(self.filename))[0]
            save_prefix = f"{base}_{self.spec.replace(' ', '')}"
        figs: Dict[str, str] = {}

        if make_corner:
            fig = corner.corner(
                flat, labels=labels, quantiles=[0.16, 0.5, 0.84], show_titles=True, title_kwargs={"fontsize": 11},
                truths=meds,
            )
            fig.tight_layout()
            path = f"{save_prefix}_corner.png"
            fig.savefig(path, dpi=200)
            fig.savefig(f"{save_prefix}_corner.pdf")
            plt.close(fig)
            figs["corner"] = path

        if make_scatter:
            # Build model at medians (mapping back to internal order)
            i = 0
            Te_log = self.Te_log if not self.fit_Te else float(meds[i]); i += int(self.fit_Te)
            Ne_log = self.Ne_log if not self.fit_Ne else float(meds[i]); i += int(self.fit_Ne)
            A12_out = None
            if self.fit_abundance and self.use_abundance_in_model and self.spec != "H I":
                A12_out = float(meds[i])
            model = self._model_vector(Te_log, Ne_log, A12_out)

            fig2 = plt.figure(figsize=(10, 5))
            x = np.arange(self.data.shape[0])
            plt.errorbar(x, self.obs_norm, yerr=self.err_norm, fmt="o", capsize=5, label="Obs")
            for j, wav in enumerate(self.data["lab_wav"].values):
                plt.annotate(f"{wav:.2f}", (x[j], self.obs_norm[j]), rotation=-45, rotation_mode="anchor")
            plt.plot(x, model, "*", ms=10, label="Model")
            plt.xlabel("Feature index")
            plt.ylabel("Relative intensity")
            plt.legend()
            path2 = f"{save_prefix}_scatter.png"
            fig2.tight_layout()
            fig2.savefig(path2, dpi=200)
            plt.close(fig2)
            figs["scatter"] = path2

        return FitResult(
            flat_samples=flat,
            medians=medians,
            intervals=intervals,
            labels=labels,
            fixed_params=fixed,
            figure_paths=figs,
        )


# Backward-compatible alias
RecomLines = Recom_Lines
__all__ = ["FitResult", "Recom_Lines", "RecomLines"]


if __name__ == '__main__':
    # Local test example (Windows path)
    ic = Recom_Lines(r'C:\Users\DELL\Desktop\PyEMILI_git\test\Abell46_erc.dat')


    # Case 1: fit abundance -> normalized to Hbeta automatically
    result = ic.fit(
        ion='O II',
        fit_abundance=True,              # normalized to Hbeta
        A12=None,
        vel=20.0,
        make_scatter=True,
        use_lab_wav=True
    )
    print('Test completed. Fit results:')
    for k, v in result.medians.items():
        try:
            print(f'{k}: {v:.3f}')
        except Exception:
            print(f'{k}: {v}')

    # Case 2 (commented): no abundance fit -> within-ion max normalization
    # result2 = ic.fit(
    #     spectrum='O II',                 $ normalized to max flux
    #     fit_abundance=False,           
    #     vel=20.0,
    #     make_scatter=True
    # )
