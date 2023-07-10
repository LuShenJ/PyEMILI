"""Several functions that may help generate the input line list.
"""

import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
from astropy.io import fits
from scipy import signal
from datetime import datetime
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from inspect import signature
import sys
import os

c = 2.9979246*10**5

def Spec_line_finding(filename, wavelength_unit='angstrom', ral_vel=0, length=100, \
    percentile=25, check_continuum=True, save_continuum=False,  vel_cor=False, snr_threshold=7, \
    prominence=4, check_lines=True, **kwargs):
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
        The radial velocity of the spectrum in the unit of km/s. Default is 0.
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
        = 2` means 2 multiplied by the continuum uncertainties. Default is 4.
    check_lines : bool, optional
        If True, an interactive plot will be presented with the spectral lines automatically found. These 
        lines will be colored by blue or red in order to distinguish the boundaries of lines. Default is True.
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
        * column 2 : flux 
        * column 3 : FWHM 
        * column 4 : SNR
        * column 5 : chi-square value
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
    flux,waves = readfile(filename, wavelength_unit=wavelength_unit, ral_vel=ral_vel, **dic_readfile)

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
    find_lines(flux, waves, con, con_std, fl_snr_threshold=snr_threshold, prominence=prominence, \
                   show_graph=check_lines, linelist_name=name, **dic_find_line)



def readfile(file, wavelength_unit='angstrom', \
             checklist=['CRVAL1','CRPIX1','CTYPE1','CDELT1','CD1_1'],\
             ral_vel = 0):
    """
    Read the spectrum from 'FITS' file or text file. If input is text file, the column 1 should
    be wavelengths, column 2 should be fluxes.

    Parameters
    ----------
    file : file, str
        File or filename to read.
    wavelength_unit : str
        The unit of the wavelength. Two types are available: `nm` or `angstrom`.
        Default is `angstrom`.
    checklist: list, optional
        The names of key arguments that used to generate the spectrum in a 'FITS' file.
        * 'CRVAL1': The wavelength value at reference pixel.
        * 'CRPIX1': The reference pixel. Default is 1.
        * 'CTYPE1': The wavelength type. Default is linear.
        * 'CDELT1': The prime argument of wavelength dispersion.
        * 'CD1_1' : The second argument of wavelength dispersion.

        NOTE: Only modify this parameter if you cannot read the spectrum.
    ral_vel : int, float, optional
        The radial velocity of the spectrum in the unit of km/s. This will shift all wavelength 
        points using the specified radial velocity. Default is 0.
    
    Returns
    -------
    flux : ndarray
        The array of fluxes of the spectrum.
    wav : ndarray
        The array of wavelengths of the spectrum.
    """

    if file.endswith('.fits') or file.endswith('.fit') or file.endswith('.fits.gz'):

        print('Try Finding Spectrum From FITS File.')
        spec = fits.open(file)
        flux = spec[0].data
        length = len(flux)
        checklist2 = [checklist in spec[0].header for checklist in checklist]


        if not checklist2[0]:
            print("Couldn't Find Wavelength Values At Reference Pixel."+ \
                 f"\nNo Keyword  {checklist[0]}.")
            sys.exit()

        else:
            wavs = spec[0].header[checklist[0]]


        if not checklist2[3]:
            if checklist[4]:
                disper = spec[0].header[checklist[4]]

            else:
                print("Couldn't Find Wavelength Dispersion"+ \
                     f"No Keyword {checklist[3]} & {checklist[4]}.")
                sys.exit()

        else:
            disper = spec[0].header[checklist[3]]


        if not checklist2[1]:
            print("Warning: Couldn't Find Reference Pixel"+ \
                 f"No Keyword {checklist[1]} Setting To 1.0.")
            repi = 1.0

        else:
            repi = spec[0].header[checklist[1]]


        if not checklist2[2]:
            print("Warning: Couldn't Find Wavelength Type" + \
                 f"No Keyword {checklist[2]} Setting To Linear.")
            wav = np.linspace(wavs-(repi-1)*disper,wavs+(length-repi)*disper,num=length)

        elif 'LIN' in spec[0].header[checklist[2]] or 'Lin' in spec[0].header[checklist[2]]:
            wav = np.linspace(wavs-(repi-1)*disper,wavs+(length-repi)*disper,num=length)

        else:
            # wav = np.array([wavs*np.exp((i-repi)*disper/wavs) for i in range(1,length+1)])
            wav = np.linspace(wavs-(repi-1)*disper,wavs+(length-repi)*disper,num=length)

        
        if wavelength_unit=='angstrom':
            wav = (c-ral_vel)*wav/c

        elif wavelength_unit=='nm':
            wav = (c-ral_vel)*wav/c*10


        


    else:
        print('Try Finding Spectrum From Text File.')

        try:
            spec = np.loadtxt(file)
            flux = spec[:,1]
            wav = spec[:,0]

            if wavelength_unit=='angstrom':
                wav = (c-ral_vel)*wav/c

            elif wavelength_unit=='nm':
                wav = (c-ral_vel)*wav/c*10
        except:
            print("Could Not Recognize This File.\n"+ \
            "Please Check The Format of File or Whether The FITS-format File Ends With '.fits'. ")
            sys.exit()

    print(f'Wavelength Range: {wav[0]} -- {wav[-1]}')

    return flux , wav



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


def subtract_continuum(flux, wavelength, percentile=25, multiple=2, length=100, check=True,\
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
            plt.plot(wavelength,flux)
            plt.plot(wavelength,continuum,'--',color='r')
            plt.title(con_name)
            plt.xlabel('Wavelength [$\\rm{\AA}$]')
            plt.ylabel('Relative Flux')
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


def Chi2(y, yfit, n, A):
    """
    The chi-squared distribution with 1 degree of freedom. This is used to check the degree of 
    fit of the spectral line and the Gaussian function. The lower this value is , the better the
    fit is.
    """

    return sum(((y-yfit)/A)**2)/n



# The Gaussian function
def gauss(x, mu, sigma, A):

    return A*np.exp(-(x-mu)**2/2/sigma**2)


# Fit the gaussian function with the data
def gauss_fit(x, y, peak_index):

    mu = x[peak_index]
    sigma = (x[-1] - x[0])/8
    A = y[peak_index]/2

    popt = curve_fit(gauss,x,y,p0=[mu,sigma,A])[0]
    y_fit = gauss(x,*popt)

    return popt, y_fit, Chi2(y,y_fit,len(y)-2,A)


def find_lines(flux, wavelength, continuum, continuum_unc, linelist_name=None, fl_snr_threshold=7, \
    prominence=4, show_graph=True, Chi2_threshold=0.2):
    """
    Find spectral lines based on `scipy.find_peaks`.

    Parameters
    ----------
    flux : array_like
        The fluxes of the spectrum.
    wavelength : array_like
        The wavelengths of the spectrum.
    continuum : array_like
        The continuum of the spectrum.
    continuum_unc : array_like
        The continuum uncertainties of the spectrum.
    linelist_name : str, optional
        Name of the file that saves the line list. If None, the current time will be used as the `linelist_name`.
    fl_snr_threshold : float, optional
        The minimum SNR value of the spectral line to be found. Default is 7.
    prominence : float, optional
        Required prominence of peaks. See details in `scipy.signal.peak_prominences`. The parameter input 
        here is the multiple of the continuum uncertainty. e.g., `prominence
        = 2` means 2 multiplied by the continuum uncertainties. Default is 4.
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
    Chi2_threshold : float, optional
        The maximum chi-squared value for a fitted line to be included in the output line list. Default
        is 0.2.

    Returns
    -------
    out : file, 2-D ndarray
        The output line list.
        * column 1 : center wavelength 
        * column 2 : flux 
        * column 3 : FWHM 
        * column 4 : SNR
        * column 5 : chi-square value
    """

    # Extra wavelength points on each boundary of the line
    move_step = int(0.5/(wavelength[1]-wavelength[0]) + 1) + 1

    # Threshold condition of the SNR
    condi = (flux/continuum_unc >fl_snr_threshold)|(flux/continuum_unc <-fl_snr_threshold)

    # Extract all boundaries (index of True)
    margin_bool = np.diff(condi) 

    # Check if the first wavelength point or the last wavelength point is a line boundary
    if condi[0] :
        margin_bool = np.insert(margin_bool,0,1)

    else:
        margin_bool = np.insert(margin_bool,0,0)

    if condi[-1] :
        margin_bool = np.insert(margin_bool,-1,1)

    else:
        margin_bool = np.insert(margin_bool,-1,0)

    # Reshape the indexes of boundaries
    margin_index = np.argwhere(margin_bool).reshape(-1,2)

    # Each boundary moves a `move_step`
    margin_index[:,0] = margin_index[:,0] - move_step 
    margin_index[:,1] = margin_index[:,1] + move_step 

    # plt.plot(wavelength,flux,'grey')
    # plt.vlines(wavelength[margin_index],ymin=-1,ymax=1)
    # plt.show()

    # Ensure the indexes do not exceed the wavelength range of spectrum
    r_beyond = np.sum((margin_index>=len(wavelength)) == True)
    l_beyond = np.sum((margin_index<0) == True)
    margin_index[margin_index<0] = np.arange(l_beyond)
    margin_index[margin_index>=len(wavelength)] = \
    np.linspace(len(wavelength)-r_beyond,len(wavelength)-1,num=r_beyond)

    # Combine the lines if indexes overlap
    for i in range(len(margin_index)-1):
        if margin_index[i,1] >= margin_index[i+1,0]:
            margin_index[i,1] = margin_index[i+1,1]
            margin_index[i+1,0] = margin_index[i,0]

    for i in range(len(margin_index)):
        l_condi = margin_index[:,0]==margin_index[i,0]
        margin_index[l_condi,1] = max(margin_index[l_condi,1])

    margin_index = np.unique(margin_index).reshape(-1,2)
    margin = wavelength[margin_index]

    if show_graph:
        fig = plt.figure(figsize=(16,9))
        main = plt.plot(wavelength,flux,'grey')

    output = []
    wavindice = []
    for i in range(len(margin_index)):

        # For each group of line boundaries
        condi = (wavelength>=margin[i,0])&(wavelength<=margin[i,1])
        subwav = wavelength[condi]
        subflux = flux[condi]
        subcon = continuum_unc[condi]
        
        # Find each peak according to the prominence
        peaks ,properties= find_peaks(abs(subflux),prominence=subcon*prominence)

        # If just one peak is found
        if peaks.size == 1:

            try:
                # Fit the profile by Gaussian function
                popt, y_fit, res = gauss_fit(subwav,subflux,int(peaks))

                # If the chi-square is too large
                if res > Chi2_threshold:
                    continue

                line_center = popt[0]
                snr = float(subflux[peaks]/subcon[peaks])
                margin_left = margin_index[i,0]
                margin_right = margin_index[i,1]
                peak_flux = subflux[peaks]

                if sum(subflux) < 0:
                    sub_c = continuum[condi]
                    sumFlux = sum(subflux/sub_c)*(subwav[1]-subwav[0])

                else:
                    sumFlux = sum(subflux)

                fwhm = abs(2.355*popt[1])
                output.append([line_center,snr,res,margin_left,\
                    margin_right,y_fit,peak_flux,sumFlux,fwhm])

            except:#RuntimeError:
                continue

        # If more than one peak is found        
        elif peaks.size > 1:
            
            # For each peak, determine its boundaries.
            bases = np.unique(np.hstack([properties['left_bases'],properties['right_bases']]))

            for k,j in zip(peaks,range(len(bases)-1)):

                try:
                    left = bases[j]
                    right = bases[j+1]
                    popt, y_fit, res = gauss_fit(subwav[left:right+1],subflux[left:right+1],k-left)

                    if res > Chi2_threshold:
                        continue

                    line_center = popt[0]
                    snr = float(subflux[k]/subcon[k])
                    margin_left = margin_index[i,0]+left
                    margin_right = margin_index[i,0]+right
                    peak_flux = subflux[k]

                    if sum(subflux[left:right+1]) < 0:
                        sub_c = continuum[left:right+1]
                        sumFlux = sum(subflux[left:right+1]/sub_c)*(subwav[1]-subwav[0])

                    else:
                        sumFlux = sum(subflux[left:right+1])

                    fwhm = abs(2.355*popt[1])
                    output.append([line_center,snr,res,margin_left,\
                                   margin_right,y_fit,peak_flux,sumFlux,fwhm])


                except:#RuntimeError:
                    continue
        
        else:
            continue



        # Save the flux points if fit is successful
        wavindice += list(range(margin_index[i,0],margin_index[i,1]+1))

    # Delete the flux and wavelength points of these lines, the rest points are treated as continuum
    newwav = np.delete(wavelength,wavindice)
    newflux = np.delete(flux,wavindice)

    if len(newwav) == 0:
        print('The input SNR is too small!')
        sys.exit()


    delnum = []
    c = ['r','b']
    num = 0
    lines = []
    # Re-determine the SNR of each line found
    for i in range(len(output)):
        length = int(5/(wavelength[1]-wavelength[0]))
        left_b = max(0,np.argmin(abs(wavelength[output[i][3]]-newwav))-length)
        right_b = min(len(newwav)-1,np.argmin(abs(wavelength[output[i][4]]-newwav))+length)
        con_std = np.std(newflux[left_b:right_b])
        newsnr = float(output[i][6]/con_std)
        output[i][1] = newsnr

        if newsnr <= fl_snr_threshold-1 and newsnr >= -fl_snr_threshold+1:
            delnum.append(i)

        else:
            if show_graph:
                num += 1
                auto_line = plt.plot(wavelength[output[i][3]:output[i][4]+1],output[i][5],\
                                    '--',color=c[num%2])
                auto_line1 = plt.text(output[i][0],output[i][6]*1.1,\
                                    f'{output[i][0]:.3f}',rotation=90)
                lines.append([wavelength[output[i][3]],wavelength[output[i][4]],\
                              auto_line,auto_line1,output[i][0]])


    output = np.array(output,dtype=object)
    output = np.delete(output,delnum,0).tolist()


    #output [line_center, snr, res, left_indice, right_indice, y_fit, peakflux, totalflux, fwhm]
    edge = []

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
                    ldx = np.argmin(abs(wavelength-edge[0]))
                    rdx = np.argmin(abs(wavelength-edge[1]))
                    x = wavelength[ldx:rdx+1]
                    y = flux[ldx:rdx+1]
                    peakdx = int((rdx+ldx)/2)
                    lineflux = sum(y)

                    if lineflux <= 0:
                        peaky = min(y)
                        lineflux = sum(y/continuum[ldx:rdx+1])*(x[1]-x[0])

                    else:
                        peaky = max(y)

                    snr = float(peaky/continuum_unc[ldx:rdx+1][np.where(y==peaky)])

                    try:
                        popt, y_fit, res = gauss_fit(x,y,peakdx-ldx)
                        line = plt.plot(x,y_fit,'--',color='orange')
                        line1 = plt.text(popt[0],peaky*1.1,f'{popt[0]:.3f}',rotation=90)
                        lcenter = popt[0]
                        output.append([lcenter,snr,res,ldx,rdx,y_fit,max(y), \
                                       lineflux,abs(2.355*popt[1])])
                        if lineflux > 0:
                            print(f'Successful line fitting. Line center: {lcenter:.3f}.'+ \
                                  f'Line flux : {lineflux:.3e}.')
                        else:
                            print(f'Successful line fitting. Line center: {lcenter:.3f}.'+ \
                                  f'Equivalent Width : {lineflux:.3e}.')

                    except:
                        lcenter = x[peakdx-ldx]
                        print(f'Failed to fit from range {wavelength[ldx]}--{wavelength[rdx]}.'+ \
                              f'Taking the middle {lcenter:.3f} as the line center.'+ \
                              f'Line flux : {lineflux:.3e}.')
                        output.append([lcenter,snr,0,ldx,rdx,None,max(y),lineflux,0])
                        line = plt.plot(x,y,'--',color='orange')
                        line1 = plt.text(lcenter,peaky*1.1,f'{lcenter:.3f}',rotation=90)
 

                    lines.append([wavelength[ldx],wavelength[rdx],line,line1,lcenter])
                    edge.clear()
                    plt.show()
        
        # Delete lines
        if event.key == 'd':
            left = np.array([i[0] for i in lines])
            right = np.array([i[1] for i in lines])
            lcenter = (left+right)/2
            del_ix = np.argmin(abs(event.xdata-lcenter))

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

    print(f'Total of {len(line_c)} lines found.')

    out = np.array(np.stack((line_c,line_f,line_w,line_s,line_r),axis=-1),dtype=np.float64)
    out = out[out[:,0].argsort()]
    
    if not linelist_name:
        linelist_name = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')

    np.savetxt(f'{linelist_name}.txt',out,fmt='%-12.3f %-12.2e %-12.2f %-12.1f %-12.4f')
    


if __name__ == "__main__":

    Spec_line_finding(r"C:\Users\DELL\Desktop\PyEMILI1\V1405Cas_2021May3.fits",prominence=6)
