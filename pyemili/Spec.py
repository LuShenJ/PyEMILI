"""Several functions that may help generate the input line list.
"""

import numpy as np
import matplotlib
matplotlib.use("TkAgg")
matplotlib.rcParams.update({'font.size': 20})
from matplotlib import pyplot as plt
from scipy import signal
from datetime import datetime
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from inspect import signature
import sys
import os
from tqdm import tqdm

c = 2.9979246*10**5

def Spec_line_finding(filename, wavelength_unit='angstrom', ral_vel=0, length=100, \
    percentile=25, check_continuum=True, save_continuum=False,  vel_cor=False, snr_threshold=7, \
    prominence=6, check_lines=True, append=False, **kwargs):
    """ 
    An integrated function for spectral lines finding. The processes include: 
    * (1) Reading the spectrum from specified file.
    * (2) Fitting and then subtracting the continuum automatically.
    * (3) Correcting the radial velocity if needed.
    * (4) Finding the spectral lines.

    Only the prime parameters are listed here, more options can be specified by 
    the `**kwargs` parameter.

    Parameters
    ----------
    filename : str, file-like or `pathlib.Path`
        The spectral file to be read.
    wavelength_unit : str
        The unit of the wavelength. Two types are available: `nm` or `angstrom`.
        Default is `angstrom`.
    ral_vel : int, float, optional
        The radial velocity of the spectrum in units of km/s. Default is 0.
    length : float, optional
        The length of the moving window used to compute the continuum. A higher spectral 
        resolution generally needs a longer moving window. Default is 100 angstroms. This 
        Parameter is best set to 3-7 times the maximum full width of line.
    percentile : float, optional
        The percentile to compute the values of continuum in a moving window. Default is 25
        for pure emission line spectrum. 75 is suggested for absorption line spectrum. This 
        parameter must be between 0 and 100.
    check_continuum : bool, optional
        Whether to check the computed continuum in a plot. Default is True.
        If True, you will see a plot of the spectrum with the computed continuum and a
        command in the terminal. Follow what it says in the terminal, you can change the 
        testing parameters.
    save_continuum : bool, optional
        If True, save the plot of the spectrum with the computed continuum. Default is False.
    vel_cor : bool, optional
        If True, correct the radial velocity using the `correct_velocity` function.
        Default is False.
    snr_threshold : float, optional
        The minimum SNR value of the spectral line to be found. Default is 7.
    prominence : float, optional
        Required prominence of peaks. See details in `scipy.signal.peak_prominences`. The parameter input 
        here is the multiple of the continuum uncertainty. e.g., `prominence
        = 2` means 2 multiplied by the continuum uncertainties. Default is 6.
    check_lines : bool, optional
        If True, an interactive plot will be presented with the spectral lines automatically found. These 
        lines will be colored by blue or red in order to distinguish the boundaries of lines. Default is True.
    append : bool, optional
        If True, instead of overwriting the saved line list file, the line list will be added starting from 
        the last line of the line list file.
        NOTE: 
        * Place the cursor on the boundary of the line and press the keyboard 'X' to determine the boundaries 
          of the line you want to add. After pressing 'X' twice, the fluxes in covered wavelengths will 
          be fitted by Gaussian function. And the details of this line will be add in the output line 
          list if the fit is successful.
        * Place the cursor within the wavelength of line you want to delete and press the keyboard 'D' 
          to delete this line. Lines found automatically can also be deleted.

    Returns
    -------
    out : file, 2-D ndarray
        The output line list.
        * column 1 : center wavelength 
        * column 2 : wavelength error
        * column 3 : flux
        * column 4 : flux error
        * column 5 : FWHM (km/s)
        * column 6 : SNR
    """

    # Initialize the optional parameters in each function
    dic_readfile = {}
    dic_subtract_c = {}
    dic_c_vel = {}
    dic_find_line = {}

    # If optional parameters are input
    if len(kwargs) != 0:

        for key in kwargs:

            if key in list(signature(readfile).parameters):
                dic_readfile[key] = kwargs.get(key)
                break

            if key in list(signature(subtract_continuum).parameters):
                dic_subtract_c[key] = kwargs.get(key)
                break

            if key in list(signature(correct_velocity).parameters):
                dic_c_vel[key] = kwargs.get(key)
                break

            if key in list(signature(find_lines).parameters):
                dic_find_line[key] = kwargs.get(key)
                break

    # Split the path and the filename 
    basename = os.path.basename(filename)
    name = os.path.splitext(basename)[0]

    # Read the spectrum
    flux, waves, flux_err = readfile(filename, wavelength_unit=wavelength_unit, ral_vel=ral_vel, **dic_readfile)
    # flux = flux*1e14
    # Fit and subtract the continuum
    flux,con,con_std= subtract_continuum(flux, waves, check=check_continuum, save=save_continuum,\
                                         length=length, percentile=percentile, con_name=name, \
                                        **dic_subtract_c)

    # Correct the radial velocity if True
    if vel_cor:
        waves, v, vstd = correct_velocity(flux, waves, con_std, **dic_c_vel)
        print(f'Mean Radial Velocity: {v:.2f}')
        print(f'Standard deviation of Radial Velocity: {vstd:.2f}')

    # Find the spectral lines
    find_lines(flux, waves, con, con_std, flux_err, fl_snr_threshold=snr_threshold, prominence=prominence, \
                   show_graph=check_lines, linelist_name=name, append=append, **dic_find_line)




def readfile(file, wavelength_unit='angstrom', ral_vel=0):
    """
    Read the spectrum from a text/ASCII file.

    The input file must have at least three columns:
    * column 1 : wavelength (in Angstrom or nm, set by ``wavelength_unit``)
    * column 2 : flux
    * column 3 : flux uncertainty

    Parameters
    ----------
    file : file, str
        File or filename to read.
    wavelength_unit : str
        The unit of the wavelength. Two types are available: ``nm`` or ``angstrom``.
        Default is ``angstrom``.
    ral_vel : int, float, optional
        The radial velocity of the spectrum in units of km/s. This will shift all wavelength 
        points using the specified radial velocity. Default is 0. In units of km/s.

    Returns
    -------
    flux : ndarray
        The array of fluxes of the spectrum.
    wav : ndarray
        The array of wavelengths of the spectrum.
    flux_err : ndarray
        The array of flux uncertainties of the spectrum.
    """
    print('Reading spectrum from text/ASCII file.')

    try:
        spec = np.loadtxt(file)
    except Exception as e:
        print(e)
        print("Could not recognize this file. Please check that it is a plain text/ASCII file.")
        sys.exit()

    if spec.ndim != 2 or spec.shape[1] < 3:
        print("Input text file must have at least three columns: wavelength, flux, flux_err.")
        sys.exit()

    wav = spec[:, 0]
    flux = spec[:, 1]
    flux_err = spec[:, 2]

    if wavelength_unit == 'angstrom':
        wav = (c - ral_vel) * wav / c
    elif wavelength_unit == 'nm':
        wav = (c - ral_vel) * wav / c * 10

    print(f'Wavelength Range: {wav[0]:.3f} -- {wav[-1]:.3f}')

    return flux, wav, flux_err


def _subtract_continuum(flux, wavelength, length, percentile=25):
    """
    The helper function of `subtract_continuum`.
    """

    # Window length 
    length = length/(wavelength[1]-wavelength[0])
    half_len = int(length/2)
    continuum = np.zeros_like(flux)

    for i in range(len(flux)):

        if i < half_len:
            subspec = flux[0:i+half_len]

        elif i >= half_len and i <= len(continuum)-half_len:
            subspec = flux[i-half_len:i+1+half_len]

        else:
            subspec = flux[i-half_len:len(flux)-1]

        continuum[i] = np.percentile(subspec,percentile)

    return continuum


def subtract_continuum(flux, wavelength, percentile=25, multiple=3, length=100, check=True,\
    save=False, con_name=None):
    """
    Fit the continuum and do the subtraction. Also calculate the appropriate uncertainties 
    of the continuum in each wavelength point. 

    Parameters
    ----------
    flux : array_like
        The original fluxes of the spectrum.
    wavelength : array_like
        The wavelengths of the spectrum.
    percentile : float, optional
        The percentile to compute the values of continuum in a moving window. Default is 25
        for pure emission line spectrum. 75 is suggested for absorption line spectrum. This 
        parameter must be between 0 and 100.
    multiple : float, optional
        The multiple of the continuum values used to compute the uncertainties of 
        continuum. The larger this parameter is, the more flux points are used to compute.
        Default is 2. This parameter must be greater than 0.
    length : float, optional
        The length of the moving window used to compute the continuum. A higher spectral 
        resolution generally needs a longer moving window. Default is 100 angstroms. This 
        Parameter is best set to 3-7 times the maximum full width of line.
    check : bool, optional
        Whether to check the computed continuum in a plot. Default is True.
        If True, you will see a plot of the spectrum with the computed continuum and a
        command in the terminal. Follow what it says in the terminal, you can change the 
        testing parameters.
    save : bool, optional
        Whether save the plot of the spectrum with the computed continuum. Default is False.
    con_name : str, optional
        This is the title and filename of the plot if `check` and `save` are True respectively. 
        Otherwise, the current time will be used as the `con_name`.

    Returns
    -------
    subtracted_flux : ndarray
        The fluxes after subtracting the continuum.
    continuum : ndarray
        The array of continuum.
    continuum_unc : ndarray
        The array of continuum uncertainties.

    """

    # Compute the continuum of the spectrum
    continuum = _subtract_continuum(flux, wavelength, length, percentile)
    
    # If plot is available
    if check:

        # Set the `name` if not given
        if not con_name:
            con_name = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')
        yn = 'n'

        # Create interactive plot
        plt.ion()
        fig = plt.figure(figsize=(16,9))

        while yn == 'n':
            # plt.plot(wavelength,flux,'grey')
            plt.step(wavelength,flux, where='mid',color='black')
            plt.plot(wavelength,continuum,'--',color='r')
            # plt.title(con_name)
            plt.xlabel(r'Wavelength [$\rm{\AA}$]',fontsize=22)
            plt.ylabel('Relative Flux',fontsize=22)
            # plt.ylabel('Flux [$10^{-14}\\,\\rm{ergs\\,cm^{-2}\\,s^{-1}\\,\AA^{-1}}$]',fontsize=22)
            yn = ''
            yn = input("Press 'y' to finish or 'n' to reset the parameters: ")

            # When the input str is neither 'n' nor 'y'
            while yn != 'n' and yn != 'y':
                yn = ''
                yn = input("Press 'y' to finish or 'n' to reset the parameters: ")
            
            # When user thinks the current continuum isn't appropriate
            if yn == 'n':
                cpercentile = input("Enter 'percentile': ")

                if cpercentile:
                    percentile = float(cpercentile)

                clen = input("Enter 'length': ")

                if clen:
                    length = float(clen)

                # Re-compute the continuum
                continuum = _subtract_continuum(flux, wavelength, length, percentile)
            
            if yn == 'y' and save:
                fig.savefig(f'{con_name}_fitcon.png',dpi=120)
                plt.ioff()
            plt.clf()
        plt.ioff()
        plt.close()

    # The fluxes after subtracting the continuum
    subtracted_flux = flux - continuum

    continuum_unc = np.zeros_like(subtracted_flux)
    std_len = int((length/(wavelength[1]-wavelength[0])))

    for i in range(len(continuum_unc)):
        if i < std_len:
            subspec = subtracted_flux[0:i+std_len]

        elif i >= std_len and i <= len(subtracted_flux)-std_len:
            subspec = subtracted_flux[i-std_len:i+1+std_len]

        else:
            subspec = subtracted_flux[i-std_len:len(subtracted_flux)-1]

        median = np.percentile(abs(subspec),25)

        # Compute the continuum uncertainty of each wavelength point 
        continuum_unc[i] = np.std(subspec[(subspec<median*multiple)&\
                                          (subspec>-median*multiple)])

    return subtracted_flux, continuum, continuum_unc



def correct_velocity(flux, wavelength, continuum_unc, cv_wav_threshold=7, cv_snr_threshold=10):
    """
    Correct the radial velocity of the spectrum based on some specified lines from 'baselines.dat' 
    file. The specified lines should be the strongest lines around their wavelengths. Code will read
    the wavelengths from the file and find the peak around these wavelengths.

    Parameters
    ----------
    flux : array_like
        The fluxes of the spectrum.
    wavelength : array_like
        The wavelengths of the spectrum.
    continuum_unc : array_like
        The continuum uncertainties of the spectrum.
    cv_wav_threshold : float, optional
        The maximum wavelength difference to search the peak. Default is 7 angstroms.
    cv_snr_threshold : float, optional
        The minimum SNR value for a detected peak to be included in calculating radial velocity.
        Default is 10.
    
    Returns
    -------
    cor_wav : ndarray
        The array of wavelengths after correcting the radial velocity.
    mean_v : float
        The mean of the radial velocities.
    v_std : float
        The standard deviation of the radial velocities.
    """
    rootdir = os.path.abspath(os.path.join(os.path.dirname( \
                       os.path.abspath(__file__)), os.pardir))
    # Determine the length of Hann window
    num = int(9/(wavelength[1]-wavelength[0]))
    num = num+1 if num%2==0 else num
    win = signal.windows.hann(num)

    # Smooth the spectrum
    filter = signal.convolve(flux,win,mode='same')/sum(win)
        
    baselines = np.loadtxt(os.path.join(rootdir,'pyemili','Line_dataset','baselines.dat'))
    ob_wav = np.zeros_like(baselines)

    for i in range(len(baselines)):
        # Find if the base line is in the wavelength range
        if baselines[i] > wavelength[0] and baselines[i] < wavelength[-1]:
            condi = (wavelength >= baselines[i] - cv_wav_threshold) & \
                    (wavelength <= baselines[i] + cv_wav_threshold)

            subflux = filter[condi]
            subcon = continuum_unc[condi]
            snr = max(subflux)/subcon[np.argmax(subflux)]

            if snr >= cv_snr_threshold:
                # The observed wavelength is determined by the peak of fluxes
                ob_wav[i] = wavelength[condi][np.argmax(subflux)]

    # Compute the wavelength differences
    diff = c*(ob_wav[ob_wav!=0]-baselines[ob_wav!=0])/baselines[ob_wav!=0]

    mean_v = np.mean(diff)
    v_std = np.std(diff-mean_v)
    cor_wav = c*wavelength/(c+mean_v)


    return cor_wav , mean_v , v_std 



# The one/multi-Gaussian function
def multi_gauss(x, *p0):

    num = int(len(p0)/3)
    x = np.tile(x,(num,1))
    p0 = np.asarray(p0).reshape(-1,3)
    A = p0[:,2].reshape(-1,1)
    mu = p0[:,0].reshape(-1,1)
    sigma = p0[:,1].reshape(-1,1)
    y_fit = np.sum(A*np.exp(-(x-mu)**2/2/sigma**2),axis=0)

    return y_fit


# Fit the one/multi-Gaussian function with the data

def multi_gauss_fit(x, y, yerr, peak_index):
    """
    Fit one or multiple Gaussian profiles to the data using non-linear least squares,
    taking into account the flux uncertainties.

    Parameters
    ----------
    x : array_like
        Wavelength array.
    y : array_like
        Flux array (continuum-subtracted).
    yerr : array_like
        1-sigma uncertainties on the flux values at each wavelength point.
    peak_index : array_like
        Indices of the initial peak positions used to initialize the Gaussian centers.

    Returns
    -------
    popt : ndarray
        Optimized parameters of the Gaussian profiles, flattened as
        [mu1, sigma1, A1, mu2, sigma2, A2, ...].
    errs : ndarray
        1-sigma uncertainties on each parameter in ``popt``.
    """
    num = len(peak_index)
    p0 = []
    for i in range(num):
        mu = x[peak_index[i]]
        sigma = (x[-1] - x[0]) / 8 / num
        A = y[peak_index[i]] / 2
        p0.extend([mu, sigma, A])

    p0 = np.array(p0).reshape(-1, 3)

    yerr = np.asarray(yerr)
    if yerr.shape != y.shape:
        raise ValueError("yerr must have the same shape as y.")

    # Avoid zero or negative uncertainties
    if np.any(yerr <= 0):
        positive = yerr[yerr > 0]
        fallback = np.median(positive) if positive.size > 0 else 1.0
        yerr = np.where(yerr <= 0, fallback, yerr)

    popt, pcov = curve_fit(multi_gauss, x, y, p0=p0, sigma=yerr, absolute_sigma=True)
    errs = np.sqrt(np.diag(pcov))

    return popt, errs




def _get_fit_params(wavelength, subwav, subflux, continuum, subcon, suberr,
                    prominence, fwhm_threshold, announce=False):
    """
    A sub-function for ``find_lines``.
    Find the peaks in the given range of spectrum and fit with one/multi-Gaussian function,
    taking into account the flux uncertainties.

    Parameters
    ----------
    wavelength : array_like
        Full wavelength grid of the spectrum.
    subwav : array_like
        Wavelength grid in the local window around a candidate line.
    subflux : array_like
        Continuum-subtracted fluxes in the local window.
    continuum : array_like
        Continuum values on the full wavelength grid.
    subcon : array_like
        Continuum uncertainties in the local window.
    suberr : array_like
        Flux uncertainties in the local window (from the input spectrum).
    prominence : float
        Prominence threshold (in units of ``subcon``) for peak finding.
    fwhm_threshold : list of float
        Allowed FWHM range in km/s.
    announce : bool, optional
        If True, print diagnostic messages when fitting fails.

    Returns
    -------
    output : list
        Each element contains
        [mu, snr, fluxerr, margin_left, margin_right, y_fit, peak_flux, flux, fwhm, mu_err].
        For emission lines ``flux`` is the integrated line flux.
        For absorption lines ``flux`` is the equivalent width (EW).
    """

    output = []

    # Find peaks according to the prominence
    peaks = find_peaks(abs(subflux), prominence=subcon * prominence)[0]

    if len(peaks) == 0:
        if announce:
            print('Cannot find any peak. Try reducing the "prominence".')
        return 0

    # Combine continuum noise and input flux uncertainties in quadrature
    total_err = np.sqrt(suberr**2 + subcon**2)

    # Fit the profile using one/multi-Gaussian function
    try:
        popts, errs = multi_gauss_fit(subwav, subflux, total_err, peaks)
    except Exception as e:
        if announce:
            print(e)
            print('Fitting failed.')
        return 0

    # Number of Gaussian profiles
    num = int(len(popts) / 3)

    for i in range(num):
        # Fitting parameters for each Gaussian profile
        popt = popts[3 * i:3 * (i + 1)]
        perr = errs[3 * i:3 * (i + 1)]

        mu, sigma, A = popt
        mu_err, sigma_err, A_err = perr

        snr = float(subflux[peaks[i]] / subcon[peaks[i]])

        # Boundaries of profile (±5σ around the center)
        margin_left = np.argmin(np.abs(wavelength - (mu - 5 * sigma)))
        margin_right = np.argmin(np.abs(wavelength - (mu + 5 * sigma)))
        peak_flux = subflux[peaks[i]]

        fwhm = abs(2.355 * sigma) / mu * c

        if fwhm > fwhm_threshold[1] or fwhm < fwhm_threshold[0]:
            if announce:
                print('FWHM out of range.')
            else:
                continue

        y_fit = multi_gauss(wavelength[margin_left:margin_right + 1], *popt)

        if len(y_fit) == 0:
            continue

        # If it's absorption profile -> compute equivalent width (EW)
        if A < 0:
            sub_c = continuum[margin_left:margin_right + 1]
            delta_x = np.diff(subwav).mean()
            flux = np.sum((y_fit / sub_c) * delta_x)  # EW
        else:
            # Emission line: integrated line flux
            flux = A * np.sqrt(2 * np.pi) * sigma

        # Propagate uncertainties from A and sigma; they already include the flux errors
        if A == 0 or sigma == 0:
            fluxerr = np.nan
        else:
            fluxerr = abs(flux) * np.sqrt((A_err / A) ** 2 + (sigma_err / sigma) ** 2)

        output.append([mu, snr, fluxerr, margin_left,
                       margin_right, y_fit, peak_flux, flux, fwhm, mu_err])

    return output




def find_lines(flux, wavelength, continuum, continuum_unc, flux_err,
    linelist_name=None, fl_snr_threshold=7, prominence=6, show_graph=True,
    fwhm_threshold=[5, 400], append=False):
    """
    Find spectral lines based on ``scipy.signal.find_peaks``.

    Parameters
    ----------
    flux : array_like
        The continuum-subtracted fluxes of the spectrum.
    wavelength : array_like
        The wavelengths of the spectrum.
    continuum : array_like
        The continuum of the spectrum.
    continuum_unc : array_like
        The continuum uncertainties of the spectrum.
    flux_err : array_like
        The 1-sigma uncertainties of the original flux at each wavelength point.
    linelist_name : str, optional
        Name of the file that saves the line list. If None, the current time will be used as the ``linelist_name``.
    fl_snr_threshold : float, optional
        The minimum SNR value of the spectral line to be found. Default is 7.
    prominence : float, optional
        Required prominence of peaks. See details in ``scipy.signal.peak_prominences``. The parameter input 
        here is the multiple of the continuum uncertainty. e.g., ``prominence = 2`` means 2 multiplied by
        the continuum uncertainties. Default is 6.
    show_graph : bool, optional
        If True, an interactive plot will be presented with the spectral lines automatically found. The 
        lines found automatically will be colored by blue or red in order to distinguish the boundaries of 
        lines. Default is True.
        NOTE: 
        * Place the cursor on the boundary of the line and press the keyboard 'X' to determine the boundaries 
          of the line you want to add. After pressing 'X' twice, the fluxes in covered wavelengths will 
          be fitted by Gaussian function. And the details of this line will be add in the output line 
          list if the fit is successful.
        * Place the cursor within the wavelength of line you want to delete and press the keyboard 'D' 
          to delete this line. Lines found automatically can also be deleted.
    fwhm_threshold : list of float, optional
        The FWHM range of the fitted spectral lines; out-of-range lines will be excluded. The list length 
        of ``fwhm_threshold`` must be 2, the first being the lower limit and the second the upper limit.
        Default is [8, 200] km/s.

    Returns
    -------
    out : 2-D ndarray
        The output line list.
        * column 1 : center wavelength 
        * column 2 : wavelength error
        * column 3 : flux (emission) or equivalent width (absorption)
        * column 4 : flux/EW error
        * column 5 : FWHM (km/s)
        * column 6 : SNR
    """

    # Extra wavelength points on each boundary of the line
    move_step = int(0.5 / (wavelength[1] - wavelength[0]) + 1) + 1

    # Threshold condition of the SNR
    condi = (flux / continuum_unc > fl_snr_threshold) | (flux / continuum_unc < -fl_snr_threshold)

    # Extract all boundaries (index of True)
    margin_bool = np.diff(condi) 

    # Check if the first wavelength point or the last wavelength point is a line boundary
    if condi[0]:
        margin_bool = np.insert(margin_bool, 0, 1)
    else:
        margin_bool = np.insert(margin_bool, 0, 0)

    if condi[-1]:
        margin_bool = np.insert(margin_bool, -1, 1)
    else:
        margin_bool = np.insert(margin_bool, -1, 0)

    # Reshape the indexes of boundaries
    margin_index = np.argwhere(margin_bool).reshape(-1, 2)

    # Each boundary moves a ``move_step``
    margin_index[:, 0] = margin_index[:, 0] - move_step 
    margin_index[:, 1] = margin_index[:, 1] + move_step 

    # Ensure the indexes do not exceed the wavelength range of spectrum
    r_beyond = np.sum((margin_index >= len(wavelength)) == True)
    l_beyond = np.sum((margin_index < 0) == True)
    margin_index[margin_index < 0] = np.arange(l_beyond)
    margin_index[margin_index >= len(wavelength)] =         np.linspace(len(wavelength) - r_beyond, len(wavelength) - 1, num=r_beyond)

    # Combine the lines if indexes overlap
    for i in range(len(margin_index) - 1):
        if margin_index[i, 1] >= margin_index[i + 1, 0]:
            margin_index[i, 1] = margin_index[i + 1, 1]
            margin_index[i + 1, 0] = margin_index[i, 0]

    for i in range(len(margin_index)):
        l_condi = margin_index[:, 0] == margin_index[i, 0]
        margin_index[l_condi, 1] = max(margin_index[l_condi, 1])

    margin_index = np.unique(margin_index).reshape(-1, 2)
    margin = wavelength[margin_index]

    if show_graph:
        fig = plt.figure(figsize=(16, 9))
        plt.step(wavelength, flux, where='mid', color='black')
        plt.title(linelist_name)
        plt.suptitle("Press 'x' to fit the line, 'd' to remove line, 'a' to integrate between two boundaries",fontsize=14)
        plt.xlabel(r'Wavelength [$\rm{\AA}$]', fontsize=22)
        plt.ylabel('Relative Flux', fontsize=22)

    output = []
    lines = []
    colors = ['r', 'b']
    num = 0

    for i in range(len(margin_index)):
        # For each group of line boundaries
        condi = (wavelength >= margin[i, 0]) & (wavelength <= margin[i, 1])
        subwav = wavelength[condi]
        subflux = flux[condi]
        subcon = continuum_unc[condi]
        suberr = flux_err[condi]

        # Get the fitting parameters
        out = _get_fit_params(wavelength, subwav, subflux, continuum, subcon, suberr,
                              prominence, fwhm_threshold)

        if out != 0:
            for j in out:
                output.append(j)

    # Plot the found lines
    if show_graph:
        for i in range(len(output)):
            num += 1
            auto_line = plt.plot(wavelength[output[i][3]:output[i][4] + 1], output[i][5],
                                 '--', color=colors[num % 2])
            if output[i][1] > 0:
                auto_line1 = plt.text(output[i][0], max(output[i][5]) * 1.03,
                                      f'{output[i][0]:.3f}', rotation=90,
                                      horizontalalignment='center', verticalalignment='bottom')
            else:
                auto_line1 = plt.text(output[i][0], min(output[i][5]) * 1.03,
                                      f'{output[i][0]:.3f}', rotation=90,
                                      horizontalalignment='center', verticalalignment='top')
            lines.append([wavelength[output[i][3]], wavelength[output[i][4]],
                          auto_line, auto_line1, output[i][0]])

    edge = []
    edge_int = []

    # Functions of interactive interface
    def key_event(event):

        # Add extra lines
        if event.key == 'x':

            if len(edge) != 2:
                print(f'x: {event.xdata:.3f}')
                edge.append(event.xdata)

            if len(edge) == 2:
                if edge[0] == edge[1]:
                    print('Wavelength range is too small.')
                    edge.clear()
                else:
                    edge.sort()
                    ldx = np.argmin(abs(wavelength - edge[0]))
                    rdx = np.argmin(abs(wavelength - edge[1]))
                    x = wavelength[ldx:rdx + 1]
                    y = flux[ldx:rdx + 1]
                    print(f'Total observed flux (original): {sum(y) * (np.diff(x).mean()):.3e}')
                    subcon = continuum_unc[ldx:rdx + 1]
                    suberr = flux_err[ldx:rdx + 1]

                    out = _get_fit_params(wavelength, x, y, continuum, subcon, suberr,
                                          prominence, fwhm_threshold, announce=True)

                    if out != 0:
                        print('Successful line fitting.')
                        for j in out:
                            line = plt.plot(wavelength[j[3]:j[4] + 1], j[5], '--',
                                            color='green', linewidth=2)
                            if j[1] > 0:
                                line1 = plt.text(j[0], max(j[5]) * 1.03, f'{j[0]:.3f}',
                                                 rotation=90, horizontalalignment='center',
                                                 verticalalignment='bottom')
                            else:
                                line1 = plt.text(j[0], min(j[5]) * 1.03, f'{j[0]:.3f}',
                                                 rotation=90, horizontalalignment='center',
                                                 verticalalignment='top')
                            output.append(j)
                            lines.append([wavelength[j[3]], wavelength[j[4]],
                                          line, line1, j[0]])

                            if j[7] > 0:
                                print(f'Line center: {j[0]:.3f}. Line fitted flux : {j[7]:.3e}.')
                            else:
                                print(f'Line center: {j[0]:.3f}. Equivalent Width : {j[7]:.3e}.')

                    edge.clear()
                    plt.show()


        # Integrate flux directly between two boundaries
        if event.key in ['a', 'A']:

            if event.xdata is None:
                return

            if len(edge_int) != 2:
                print(f'a: {event.xdata:.3f}')
                edge_int.append(event.xdata)

            if len(edge_int) == 2:
                if edge_int[0] == edge_int[1]:
                    print('Wavelength range is too small.')
                    edge_int.clear()
                else:
                    edge_int.sort()
                    ldx = np.argmin(abs(wavelength - edge_int[0]))
                    rdx = np.argmin(abs(wavelength - edge_int[1]))

                    x = wavelength[ldx:rdx + 1]
                    y = flux[ldx:rdx + 1]
                    subcon = continuum_unc[ldx:rdx + 1]
                    suberr = flux_err[ldx:rdx + 1]
                    total_err = np.sqrt(suberr**2 + subcon**2)

                    if len(x) < 2:
                        print('Too few points in the selected region.')
                        edge_int.clear()
                        return

                    delta_x = np.diff(x).mean()

                    # Directly integrated line flux (continuum-subtracted)
                    int_flux = np.sum(y * delta_x)
                    int_err = np.sqrt(np.sum((total_err * delta_x)**2))

                    # Determine center, width and SNR
                    idx_peak = np.argmax(np.abs(y))
                    mu = x[idx_peak]
                    mu_err = delta_x / 2.0
                    fwhm = abs(edge_int[1] - edge_int[0]) / mu * c if mu != 0 else 0.0

                    # For emission use integrated flux; for absorption compute EW
                    if int_flux >= 0:
                        flux_val = int_flux
                        flux_err_val = int_err
                    else:
                        sub_c = continuum[ldx:rdx + 1]
                        ew = np.sum(y / sub_c * delta_x)
                        ew_err = np.sqrt(np.sum(((total_err / sub_c) * delta_x)**2))
                        flux_val = ew
                        flux_err_val = ew_err

                    snr = flux_val / flux_err_val if flux_err_val > 0 else 0.0
                    peak_flux = np.max(y) if flux_val >= 0 else np.min(y)

                    # Build a "fit" profile as the actual data segment
                    y_fit = y.copy()

                    new_entry = [mu, snr, flux_err_val, ldx, rdx,
                                y_fit, peak_flux, flux_val, fwhm, mu_err]
                    output.append(new_entry)

                    # Plot the integrated region and label
                    line = plt.plot(x, y_fit, '--', color='magenta', linewidth=2)
                    if flux_val > 0:
                        text_y = max(y_fit) * 1.03
                        va = 'bottom'
                    else:
                        text_y = min(y_fit) * 1.03
                        va = 'top'

                    line1 = plt.text(mu, text_y, f'{mu:.3f}', rotation=90,
                                    horizontalalignment='center', verticalalignment=va)

                    lines.append([wavelength[ldx], wavelength[rdx], line, line1, mu])

                    if flux_val > 0:
                        print(f'Line center: {mu:.3f}. Integrated line flux : {flux_val:.3e} ± {flux_err_val:.3e}.')
                    else:
                        print(f'Line center: {mu:.3f}. Equivalent Width : {flux_val:.3e} ± {flux_err_val:.3e}.')

                    edge_int.clear()
                    plt.show()
        # Delete lines
        if event.key == 'd':
            left = np.array([i[0] for i in lines])
            right = np.array([i[1] for i in lines])
            lcenter = (left + right) / 2
            del_ix = np.argmin(abs(event.xdata - lcenter))

            if event.xdata < left[del_ix] or event.xdata > right[del_ix]:
                print("Can't find a line coverd by cursor")
            else:
                print(f"Line {lines[int(del_ix)][4]:.3f} has been deleted.")
                lines[int(del_ix)][2][0].remove()
                lines[int(del_ix)][3].remove()
                del lines[int(del_ix)]
                del output[int(del_ix)]
                plt.show()

    if show_graph:
        fig.canvas.mpl_connect('key_press_event', key_event)
        plt.show()

    line_c = [i[0] for i in output]
    line_s = [i[1] for i in output]
    line_r = [i[2] for i in output]
    line_f = [i[7] for i in output]
    line_w = [i[8] for i in output]
    line_m = [i[9] for i in output]

    print(f'Total of {len(line_c)} lines found.')

    out = np.array(np.stack((line_c, line_m, line_f, line_r, line_w, line_s), axis=-1),
                   dtype=np.float64)
    out = out[out[:, 0].argsort()]

    if not linelist_name:
        linelist_name = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')

    if append:
        with open(f'{linelist_name}_linelist.txt', 'a') as f:
            np.savetxt(f, out,
                       fmt='%-12.3f %-12.3f %-12.2e %-12.2e %-12.2f %-12.1f')  
    else:
        np.savetxt(f'{linelist_name}_linelist.txt', out,
                   fmt='%-12.3f %-12.3f %-12.2e %-12.2e %-12.2f %-12.1f',
                   header='wavelength \t wav_err \t flux \t flux_err \t FWHM \t snr')




if __name__ == "__main__":

    Spec_line_finding(r"C:\Users\DELL\Desktop\PyEMILI1\V1405Cas_2021May3.txt",prominence=6)
