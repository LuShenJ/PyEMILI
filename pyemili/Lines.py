"""
Spectral line identifier. The most important function in PyEMILI. This line identification 
algorithm follows the general idea of EMILI (Sharpee et al. 2003), where many details have 
been corrected and updated.
"""

import numpy as np
import pandas as pd
import sys
from tqdm import tqdm
from pyemili.numba_func import *
import os 
import numba as nb

class Line_list(object):

    def __init__(self, wavelength, wavelength_error, flux, ral_vel=None, snr=None, fwhm=None):
        """
        Create an input line list that needs to be identified.

        Parameters
        ----------
        wavelength : array_like
            Input array of wavelengths of unidentified lines in the unit of Angstroms.
        wavelength_error : int or float or 1-D or 2-D array_like
            Input array of 1 sigma wavelength uncertainties corresponds to `wavelength`.
            This parameter can be input in 3 different ways: 

            * If input is a `int` or `float` type, this value should be in the unit of km/s,
              and will be used to determine the wavelength uncertainties on both side of observed
              lines. e.g., `10` means 10 km/s wavelength uncertainty for each side of input `wavelength`. 
            * If input is a `1-D array-like` type, the values should be in the unit of km/s,
              the first value will be used to determine the wavelength uncertainty on blue (left) side, 
              another values will be used to determine the red (right) side. e.g., `[-10,30]` means 
              10 km/s wavelength uncertainty for the blue (left) side, and 30 km/s wavelength 
              uncertainty for the red (right) side.
            * If input is a `2-D array-like` type with the shape of (n,2), these values should be in the 
              unit of Angstroms. It means you specify the wavelength uncertainty on each side of each line.

        flux : array_like
            Input array of fluxes corresponds to `wavelength`.
        ral_vel : int, float, optional
            The systematic radial velocity (typically calculated using H beta) of the input line list in 
            the unit of km/s. Default is 0.
        snr : array_like, optional
            Input array of signal-to-noise ratios corresponds to `wavelength`. This parameter
            will only be used to show along with the identification results.
        fwhm : array_like, optional
            Input array of FWHMs (full width at half maximum) corresponds to `wavelength`. This parameter
            will only be used to show along with the identification results.
        """

        # Convert Input data to ndarray and sort
        self.wav = np.array(wavelength)

        self.obs_flux = np.array(flux)

        # Symbols for elements
        self.elesym = ['H','He','Li','Be','B','C','N','O','F','Ne',\
                    'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca',\
                    'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',\
                    'Ga','Ge','As','Se','Br','Kr']

        # Roman numerals used for marking ionization state
        self.ionstat = ['I','II','III','IV','V','VI','VII','VIII','IX','X',\
            'XI','XII','XIII','XIV','XV','XVI','XVII','XVIII','XIX','XX',\
            'XXI','XXII','XXIII','XXIV','XXV','XXVI','XXVII','XXVIII','XXIX','XXX',\
            'XXXI','XXXII','XXXIII','XXXIV','XXXV','XXXVI']

        self.rootdir = os.path.abspath(os.path.join(os.path.dirname( \
                       os.path.abspath(__file__)), os.pardir))
        # Array used for output
        self.elesymarr = np.array(self.elesym,dtype='<U2')

        self.ionstatarr = np.array(self.ionstat,dtype='<U7')
        
        # Threshold value of each energy bin, e.g., [0,13.6,24.7,54.5,100] means 5 energy bins are 
        # 0-13.6, 13.6-24.7, 24.7-54.5, 55-100, >100. The unit is eV.
        self.energybin = [0,13.6,24.7,55,100]

        # Speed of light in the unit of km/s.
        self.c = 2.9979246*1.0E+5

        # Boltzmman constant in the unit of eV/K.
        self.k = 8.617333262e-5

        if ral_vel:
            if isinstance(ral_vel,(float,int)):
                self.wav = (self.c-ral_vel)*self.wav/self.c
            else:
                print("Invalid type of 'ral_vel', set to 0.")


        # Initialize the arrays of wavelength uncertainties, snr, and fwhm
        self.waverr = self.__init_para(wavelength_error,'Wavelength_Error')

        # Parameter for checking the proper wavelength uncertainties
        self.waverr_pro = 0

        self.waverr_init = wavelength_error

        self.snr = self.__init_para(snr,'snr')

        self.fwhm = self.__init_para(fwhm,'fwhm')



    def identify(self, filename, icf=None, v_cor=None, sigma=5, Ne=10000, Te=10000, \
        I=10, deplete=None, abun_type='solar', iteration=True, match_list=None, eff_coe=True):
        """
        The prime function of this class. Identify each line in the input line list. One file with 
        complete candidates of each line and another file with the best candidates of each line will 
        be generated.

        Parameters
        ----------
        filename : str, file-like or `pathlib.Path`
            Name of generated files.
        icf : array_like, optional
            The percentage of each energy bin. The array length of this parameter must be 5 (the 
            same as energy bin). These values are percentages of particular ions of that elements 
            take up of the entire abundance for that element. Default values are 
            `[0.01,0.5,0.4,0.1,0.00001]`. NOTE: If `iteration` is `True`, the `icf` will be 
            calculated automatically based on the first iteration. If both `iteration` is `True` 
            and `icf` is specified, code will prefer the specified `icf`.
        v_cor : array_like, optional
            The velocity correction of each energy bin. The array length of this parameter must be 
            5 (the same as energy bin). Default values are `[0,0,0,0,0]`. NOTE: If `iteration` is 
            `True`, `v_cor` will be calculated automatically. If both `iteration` is `True` and 
            `v_cor` is specified, code will prefer the specified `v_cor`.
        sigma : int, optional
            The maximum times of the wavelength uncertainty used for line identifications. e.g., 
            `sigma` = 5 means candidates within 5 times wavelength uncertainty of observed line 
            will be included.
        Ne : int, optional
            The electron density in the unit of cm^-3. Default value is 10000.
        Te : int, optional
            The electron temperature in the unit of K. Default value is 10000.
        I : int, optional
            The instrumental resolution of the spectrum in the unit of km/s. Multiplet lines whose 
            wavelengths are within `I` will be treated as one multiplet line in the "Multiplet Check" part.
        deplete : str, optional
            To reduce or enhance the abundance of certain ions after the ionic abundances have 
            been calculated using the `icf` values. This parameter has up to 4 arguments specified 
            in the following order:

              1. The element to deplete or enhance.
              2. The amount to deplete or enhance (negative value enhance).
              3. The lower end of the range of ionization states to deplete or enhance.
              4. The upper end of the range of ionization states to deplete or enhance.

            Different arguments should be separated with space. Different ions should be separated
            with comma.
            
            Examples:
              * `Fe 10 1 3` This input type will reduce the abundances of Fe I, Fe II, and Fe III 
                by a factor of 10.
              * `Fe -10 1` This input type will enhance the abundance of Fe I by a factor of 10.
              * `Fe 20, Ti 10` This input type will reduce the abundances of all ion of Fe and Ti 
                by a factor of 20 and 10 respectively.

        abun_type : str, file-like or `pathlib.Path`, optional
            The input elemental abundance table. Here are 2 options: 'solar' (M. Asplund et al. 2009),
            'nebula' (Osterbrock & Ferland 2006). Default is 'solar'. Besides, you can use the filename 
            of other abundance file as input.
        iteration : bool
            If True, code will extract a sub line list with relatively robust identifications made 
            in the first iteration to calculate the new `icf` and `v_cor` models and re-identify each 
            line. Default is True.
        match_list : str, optional
            The filename of the input match line table. This is a obsolete argument. It is in conflict
            with `iteration` parameter and may be deleted in future version. See (Sharpee et al. 2003) 
            for details.
        """

        
        self.name = filename

        self.sigma = sigma

        self.Ne = Ne

        self.I = I

        if iteration:
            self.loop = 2
        else:
            self.loop = 1
        
        # The dilution factor for collisional excitation term
        self.col_cor = 0.1

        self.match_list = match_list

        self.icf = None

        self.icfuc = False

        self.v_cor = None

        self.vcoruc = False

        # Initialize the ranking A candidates list (used in last iteration)
        self.Aterm = []
        self.Awave = []
        self.Anum = []
        self.extract_done = False
        # Initialize iteration parameter
        self.i = 0

        print('Initializing database')

        # The ionization energy of each ion
        self.ionene = np.load(os.path.join(self.rootdir,'pyemili','Line_dataset','ioni_energy.npy')) 

        # The index of the energy bin that each ion locates in
        self.ele_binindex = np.load(os.path.join(self.rootdir,'pyemili','Line_dataset','elebin_index.npy')) 

        # The complete atomic transition database
        self.linedb = np.load(os.path.join(self.rootdir,'pyemili','Line_dataset','Linedb.npy'))

        # Compress the atomic transition database
        uplmt = self.wav[-1]+self.waverr[-1,1]*self.wav[-1]*3*sigma/self.c
        lowlmt = self.wav[0]+self.waverr[0,0]*self.wav[0]*3*sigma/self.c      
        self.linedb = self.linedb[(self.linedb[:,0]>=lowlmt)&(self.linedb[:,0]<=uplmt)]
        self.linedb = np.asfortranarray(self.linedb)
        # The radiative recombination coefficients table
        self.RR = np.loadtxt(os.path.join(self.rootdir,'pyemili','recom','RR_rate.dat'),skiprows=1)

        # The dielectronic recombination coefficients table
        self.DR = np.loadtxt(os.path.join(self.rootdir,'pyemili','recom','DR_rate.dat'),skiprows=1) 

        # The elemental abundances table
        if abun_type == 'solar':
            self.abunfile = pd.read_table(os.path.join(self.rootdir,'pyemili','abundance','abun_solar.dat'),\
                                delim_whitespace=True,names=['id','abun'],comment='#')
            print('Use Solar Abundance')

        elif abun_type == 'nebula':
            self.abunfile = pd.read_table(os.path.join(self.rootdir,'pyemili','abundance','abun_nebula.dat'),\
                                delim_whitespace=True,names=['id','abun'],comment='#')
            print('Use Nebula Abundance')

        elif type(abun_type) == str:
             self.abunfile = pd.read_table(abun_type,\
                                delim_whitespace=True,names=['id','abun'],comment='#') 
             print(f"Use Specified Abundance From '{abun_type}'")    

        else:
            print("Invalid input 'abun_type'.")
            sys.exit()

        self.Te = Te

        self.eff_coe = eff_coe

        if self.eff_coe:

            effcoe_ion = ['HI','HeI','HeII','CII','NII','OII','NeII']
            self.effcoe_num = np.array([[1,1],[2,1],[2,2],[6,2],[7,2],[8,2],[10,2]]) - 1
            # Initialize the ranges of Te and Ne of each specie
            self.HI_Ne = 10**np.linspace(2,14,num=13)
            self.HI_Te = 1e2*np.array([5,10,30,50,75,100,125,150,200,300])

            self.HeI_Ne = 10**np.linspace(1,14,num=14)
            self.HeI_Te = 1e3*np.arange(5,26)

            self.HeII_Ne = 10**np.linspace(2,14,num=13)
            self.HeII_Te = 1e2*np.array([5,10,30,50,75,100,125,150,200,300,500,1000])

            self.CII_Ne = np.array([1e4])
            # self.CII_Te = 10**np.linspace(2,4.6,num=27)
            self.CII_Te = np.array([500,750,1000,1250,1500,2000,2500,3500,5000,7500,\
                           10000,12500,15000,20000])

            self.OII_Ne = 10**np.linspace(2,5,num=16)
            self.OII_Te = 10**np.linspace(2,4.4,num=25)

            self.NII_Te = np.array([126,158,200,251,316,398,501,631,794,1000,1260,1580,\
                            2000,2510,3160,3980,5010,6310,7940,10000,12600,15800,20000])
            self.NII_Ne = 10**np.linspace(2,6,num=41)

            self.NeII_Ne = 10**np.array([4,6])
            self.NeII_Te = 1e2*np.array([10,20,30,50,75,100,125,150,200])
            
            # Load the effective coefficients
            effcoes = nb.typed.List()
            effcoe_ix = []

            for i in effcoe_ion:

                exec(f"self.{i}coe = np.load(os.path.join(self.rootdir,'pyemili','eff_reccoe','{i}_emiss_data.npy'), \
                     allow_pickle=True).astype(np.float64)")

                exec(f"self.{i}_ix = np.array([np.argmin(abs(self.{i}_Te-self.Te)),np.argmin(abs(self.{i}_Ne-self.Ne))])")

                exec(f"effcoes.append(self.{i}coe)")

                exec(f"effcoe_ix.append(self.{i}_ix)")

            self.effcoes = effcoes
            self.effcoe_ix = np.array(effcoe_ix)

        # The table of details of each energy level of each ion
        self.ion_tab = pd.read_csv(os.path.join(self.rootdir,'pyemili','Line_dataset','ionlevel.dat'), \
                        usecols=[0,1,4,6,7],names=['ele','stat','config','term','termnum'])

        # Initialize output string of each candidate line
        self.tab_termnum =  self.ion_tab.termnum.values

        self.tab_ele = self.ion_tab.ele.values

        self.tab_stat = self.ion_tab.stat.values

        self.tab_term = self.ion_tab.term.values.astype(str)

        self.tab_conf = self.ion_tab.config.values.astype(str)

        # Initialize the icf and v_cor models
        self.__init_icfvcor(match_list,icf,v_cor)

        # Calculate the ionic abundances of each element
        self.abun = self.__calcu_abun()

        # Deplete or enhance the abundances of certain ions
        self.deplete = deplete

        self.__deplete(deplete)

        # Calculate a absolute predicted flux of H beta
        self._H_beta_flux(1)

        # If there is absorption line, also calculate the predicted EW of H beta
        if min(self.obs_flux) < 0:

            self._H_beta_flux(0)

        print('Initializing Done')

        open(f'{self.name}.out','w')
        open(f'{self.name}.dat','w')
        f1 = open(f'{self.name}.out','a')
        f2 = open(f'{self.name}.dat','a')

        for self.i in range(self.loop):

            print(f'\nNumber of Iterations: {self.i+1}')

            if self.i == self.loop-1 :

                    self._extract_Alist()
                    print(self._icf_v_out(self.icf,self.v_cor))
                    self.extract_done = True
                    out_sum = 'PyEMILI Output File\n'+\
                            '------------------------------------\n'+\
                            f'Input Matched List: {self.match_list}\n'+\
                            f'Results List: {self.name}\n'+\
                            f'Short Results List: {self.name}.dat\n'+\
                            f'Electron Temp.: {self.Te}\n'+\
                            f'Electron Density: {self.Ne}\n'+\
                            f'Inst. Resolution: {self.I}\n'
                    if self.waverr_type != 2:
                        out_sum += f'Input Wavelength uncertainty (1 sigma): {self.waverr_init} km/s\n'
                    
                    if self.waverr_pro != 0:
                        out_sum += f'Proper Wavelength uncertainty (1 sigma): {self.waverr_pro:.3f} km/s\n'

                    f1.write(out_sum)
                    f1.write('\n')
                    f1.write(self._icf_v_out(self.icf,self.v_cor))
                    f1.write('\n')
                    f1.write('\n\n')
                    self.write(f1,f2)
                    break

            else:
                self._firstrun()
            

                


    def _H_beta_flux(self,typeh):
        """
        Initialize the predicted flux or equivalent width of H beta.
        """
        
        if typeh == 0:
            self.H_abs = self.abs_flux_formu(self.abun[0,0],8.420e+06,82259.11,8,4861.325)

        if typeh == 1:
            # voidin = np.zeros((1,9))
            # RR_H = np.array([[8.318e-11 ,0.7472 ,2.965 ,700100.0 ,0.0 ,0.0]])
            # self.H_beta = self.em_flux_formu(self.abun[0,0],self.abun[0,1],1,\
            #     102823.9,4,8,8.420e+06,RR_H,voidin,voidin,4861.325,np.array([1]),0.0141341,1)[0]
            self.H_beta_eff = self.abun[0,1]*self.HIcoe[21,self.HI_ix[0],self.HI_ix[1]]/6.626e-27/2.99e18*1e14



    def __init_para(self,para,paraname):
        """
        Initialize `snr`, `fwhm`, and `wavelength_error` parameters.
        """

        if paraname == 'snr' or paraname == 'fwhm':    

            if isinstance(para,(np.ndarray,list,tuple,pd.Series)):
                return np.array(para)

            # Default values for `snr` and `fwhm` are zero
            elif para is None:
                return np.zeros_like(self.wav)

            else:
                print(f"Invalid type of '{paraname}'.")
                sys.exit()
        
        if paraname == 'Wavelength_Error':

            waverr = np.array(para)

            # If input type is `int` or `float`
            if waverr.ndim == 0:
                self.waverr_type = 0
                return np.tile([-para,para],len(self.wav)).reshape(-1,2)

            # If input type is 1-D array-like
            if waverr.ndim == 1 and len(para) == 2:
                self.waverr_type = 1
                return np.tile(para,len(self.wav)).reshape(-1,2)

            # If input type is 2-D array-like
            elif waverr.ndim == 2 and len(para) == len(self.wav):
                self.waverr_type = 2
                return waverr/np.tile(self.wav,(2,1)).T*self.c

            else:
                print(f"Invalid type of '{paraname}'.")
                sys.exit()
    


    def __init_icfvcor(self,match_list,icf,v_cor):
        """
        Initialize the `icf` and `v_cor` parameters.
        """

        if icf is not None and (not isinstance(icf,(np.ndarray,list,tuple)) or len(icf) != 5):
            print(f'Invalid Manually Specified ICF value: {icf}')
            sys.exit()

        if v_cor is not None and (not isinstance(v_cor,(np.ndarray,list,tuple)) or len(v_cor) != 5):
            print(f'Invalid Manually Specified Velocity Matrix: {v_cor}')
            sys.exit()

        if match_list is not None and (icf is not None or v_cor is not None):
            print('WARNING: Both Match File And Manually Specified ICF Or/And Velocity Matrix Detected.'+\
            'The Manually Specified ICF Or/And Velocity Matrix Will Be Used.')

            if icf is not None:
                self.icf = np.array(icf)
                self.icfuc = True

            if v_cor is not None:
                self.v_cor = np.array(v_cor)
                self.vcoruc = True

        if match_list is None:

            if icf is None:
                print('Using Default ICF')
                self.icf = np.array([0.01,0.5,0.4,0.1,0.0001])
                self.icf = self.icf/sum(self.icf)
                
            else:
                self.icf = np.array(icf)/sum(icf)
                self.icfuc = True
            
            if v_cor is None:
                print('Using Default Velocity Matrix')
                self.v_cor = np.array([0,0,0,0,0],dtype=np.float64)

            else:
                self.v_cor = np.array(v_cor,dtype=np.float64)
                self.vcoruc = True

        else:
            self.match_lines(match_list)



    def _deplete(self,deplete):
        """
        Deplete or enhance the abundances of certain ions. See details in `Line_list.identify`.
        """

        deplete = deplete.split(' ')
        ele = self.elesym.index(deplete[0])
        factor = int(deplete[1])

        if factor == 0:
            print('Factor cannot be zero as a dividend')
            print("Please check your 'deplete' parameter")
            sys.exit()

        # Negative value enhances
        if factor < 0:
            factor = 1/abs(factor)

        # If the 3rd and 4th arguments are not given
        if len(deplete) == 2:
            self.abun[ele,:] = self.abun[ele,:]/factor

        # If the 4th argument is not given
        if len(deplete) == 3:
            ion = int(deplete[2])
            self.abun[ele,ion-1] = self.abun[ele,ion-1]/factor

        # If complete arguments are given
        if len(deplete) == 4:
            lowion = int(deplete[2])-1
            upion = int(deplete[3])
            self.abun[ele,lowion:upion] = self.abun[ele,lowion:upion]/factor    



    def __deplete(self,deplete):
        """
        The main process of `deplete`.
        """

        if deplete is None:
            pass

        elif isinstance(deplete,str):

            # Separate the command by comma
            if ',' in deplete:
                deplete = deplete.split(',')
                for i in deplete:
                    self._deplete(i)
            
            else:
                self._deplete(deplete)

        else:
            print("Invalid type of 'deplete'.")
            sys.exit()



    def __calcu_abun(self):
        """
        Calculate the abundances of each ion of each element based on the determined 
        `icf` and the input table of elemental abundance.
        """

        abun = np.zeros_like(self.ionene)
        fabun = self.abunfile.abun
        fid = self.abunfile.id.values
        
        # 'i' is the sequence number of element symbol 
        for i in range(len(self.elesym)):

            # 'j' is the sequence number of ionization state
            for j in range(len(self.elesym)+1):
                # if i == 7 and j == 3:
                #     breakpoint()
                if j - i >= 2:
                    break

                # Find the bin index that the ion locates
                binix = self.ele_binindex[i,j]

                if binix != 4:
                    
                    # Find the bin index of next higher ionized ion
                    binix_u = self.ele_binindex[i,j+1]
                
                else:
                    binix_u = 4
                
                # If the indexes are the same, just multiply the `icf` of that bin
                if binix == binix_u:
                    
                    abun[i,j] = fabun.iloc[fid == self.elesym[i]].item()* \
                                self.icf[binix]
                    
                # If bin index of next higher ionized ion is larger, calculate the
                # mean value of icfs of the two bins.
                elif binix_u > binix:
                    
                    abun[i,j] = fabun.iloc[fid == self.elesym[i]].item()* \
                                (0.5*self.icf[binix]+ \
                                sum(self.icf[binix+1:binix_u])+ \
                                0.5*self.icf[binix_u])


        # Specify values for special ions
        abun[0,1] = self.icf[1]*fabun.iloc[fid =='H'].item()
        abun[1,1] = self.icf[2]*fabun.iloc[fid =='He'].item()
        abun[1,2] = self.icf[3]*fabun.iloc[fid =='He'].item()

        return abun



    def match_lines(self,match_list=None,match_table=None):
        """
        Match the identified lines to the atomic transition database and calculate more accurate 
        models of `icf` and `v_cor` based on these identified lines. This function is used when 
        the input `match_list` exists or when `iteration` is set to True.
        """

        # File contains the lines that needs to be matched
        # match = np.loadtxt(os.path.join(self.rootdir,'pyemili','Line_dataset','match_lines.dat')) 

        # If `match_list` exists
        if match_table is None:
            lines = pd.read_table(match_list,delim_whitespace=True,\
                names=['obs_wav','lab_wav','ele','ion','flux'])

            
        # If `iteration` is set to True
        else:
            lines = pd.DataFrame(match_table,columns=['obs_wav','lab_wav','ele',\
                'ion','flux','pre_flux','effcoe'])
            
        # the sequence numbers of each line's element symbol and ionization state \
        # and transition type.
        ele = lines.ele.map( \
            lambda x: self.elesym.index(x[1:]) if x[0] == '[' \
                 else self.elesym.index(x)).rename('ele')
        ion = lines.ion.map( \
            lambda x: self.ionstat.index(x[:-1]) if x[-1]==']' \
                 else self.ionstat.index(x)).rename('ion')
        tt = lines.apply( \
            lambda x: 1 if x.ion[-1]==']' and x.ele[0]!='[' \
                 else 2 if x.ion[-1]==']' and x.ele[0]=='[' \
                 else 0,axis=1).rename('tt')

        # ion_prog = pd.concat([ion,tt,lines.flux],axis=1).apply(lambda x:\
        #     int(x.ion) if x.tt>0 or x.flux<0 else int(x.ion+1),axis=1)


        # Match the H beta if it exists
        H_beta = lines[(ele==0)&(ion==0)&(lines.effcoe==22)] #
        # Determine the H beta flux and normalized
        if len(H_beta.flux) != 0:
            
            if len(H_beta.flux) != 1:
                H_beta = H_beta[H_beta.flux==H_beta.flux.max()]
                
            if H_beta.flux.item() != 1:
                print('H beta flux is not normalized. About to normalize the flux.')
                lines.flux = lines.flux/H_beta.flux.item()
                if H_beta.flux.item() > 0:
                    self.obs_flux = self.obs_flux/H_beta.flux.item()
        # If do not match H beta successfully, consider the nearest line as H beta
        else:
            H_betaflux = self.obs_flux[np.argmin(abs(self.wav-4861.325))]
            H_beta = self.wav[np.argmin(abs(self.wav-4861.325))]
            print('H beta is not included in matched list.')
            print(f'Consider H_beta as {H_beta}, observed flux:{H_betaflux}')
            if H_betaflux != 1:
                lines.flux = lines.flux/H_betaflux
                if H_beta.flux.item() > 0:
                    self.obs_flux = self.obs_flux/H_betaflux


        tbin = pd.concat([ele,ion,tt],axis=1).apply( \
            lambda x: self.ele_binindex[x.ele,x.ion+1] if x.tt==0 \
                 else self.ele_binindex[x.ele,x.ion],axis=1)

        
        # Modify the dillution factor of collisional excitation term
        if len(lines[tt==2]) >= 5:
            self.col_cor = self.col_cor/np.median(lines[tt==2].pre_flux/lines[tt==2].flux)
                    
        # If user did not specify the `icf` parameter
        if not self.icfuc:

            # Set the minimum and maximum values of each bin
            min_val = [1.0E-3]*4 +[1.0E-4]
            max_val = [0.5]*4 + [0.3]

            # Check if there is any line in bin 1 and set icf
            if sum(tbin==0) != 0:
                flx1 = lines.flux[tbin==0].max()
                ix1 = 1.0E-1 if flx1 > 1.0E-2 \
                 else 1.0E-3 if flx1 < 1.0E-4 \
                 else 1.0E-2

            else:
                ix1 = min_val[0]

            # Check if there is any line in bin 5 and set icf
            if sum(tbin==4) != 0:
                flx5 = lines.flux[tbin==4].max()
                if flx5 > 0:
                    ix5 = flx5/10 + 1.0E-4
                else:
                    ix5 = min_val[4]

                if ix5 > max_val[4]:
                    ix5 = max_val[4]

            else:
                ix5 = min_val[4]

            # Match all He I lines with effective recombination coefficients
            HeI_obsflux = abs(lines[(ele==1)&(ion==0)&(lines.effcoe!=0)].flux)
            HeI_preflux = lines[(ele==1)&(ion==0)&(lines.effcoe!=0)].pre_flux

            # If no He I lines found, set the ratio as 0
            if len(HeI_obsflux) != 0:
                ixr1 = np.median(HeI_obsflux/HeI_preflux)
            else:
                ixr1 = 0

            # Match all He II lines with effective recombination coefficients
            HeII_obsflux = abs(lines[(ele==1)&(ion==1)&(lines.effcoe!=0)].flux)
            HeII_preflux = lines[(ele==1)&(ion==1)&(lines.effcoe!=0)].pre_flux

            # If no He I lines found, set the ratio as 0
            if len(HeII_obsflux) != 0:
                ixr2 = np.median(HeII_obsflux/HeII_preflux/ixr1)
            else:
                ixr2 = 0


            remain=1-ix1-ix5

            if not ixr1 and not ixr2:
            # Either everything is fairly low ionization
                ix4 = min_val[3]
                ix3 = min_val[2]
                ix2 = remain-ix3-ix4
            # Or everthing is rather high ionization
                if ix5 == 0.3: 
                    ix2= min_val[1]

            if ixr1 and not ixr2:
            # Probably everything is lower ionization bin 1-3
                ix4 = min_val[3]
                ix2 = (remain-ix4)/(1+ixr1)

                if ix2 < min_val[1]:
                     ix2 = min_val[1]

                ix3 = ixr1*ix2

                if ix3 < min_val[2]:
                    ix3=min_val[2]

            if ixr1 and ixr2:
            # Split everything a bit more evenly
                ix2 = remain/(1+ixr1+ixr1*ixr2)

                if ix2 < min_val[1]:
                    ix2 = min_val[1]

                ix3 = ix2*ixr1

                if ix3 < min_val[2]:
                    ix3 = min_val[2]

                ix4 = ix3*ixr2

                if ix4 < min_val[3]:
                    ix4 = min_val[3]

            if not ixr1 and ixr2:
            # Bit higher ionization
                ix2 = min_val[1]
                ix4 = (remain-ix2)/((1/ixr2)+1)

                if ix4 < min_val[3]:
                    ix4 = min_val[3]

                ix3 = ix4/ixr2

                if ix3 < min_val[2]:
                    ix3 = min_val[2]

            self.icf = np.array([ix1,ix2,ix3,ix4,ix5])

            
        rvcor = (lines.obs_wav - lines.lab_wav)*self.c/lines.lab_wav

        # # Check the wavelength uncertainty
        # sigma1 = np.std(rvcor[rvcor>=0])
        # sigma2 = np.std(rvcor[rvcor<0])
        # # self.waverr = self.__init_para([-sigma2,sigma1],'Wavelength_Error')
        # self.waverr_pro = [-sigma2,sigma1]
        self.waverr_pro = np.std(rvcor)

        # If user did not specify the `v_cor` parameter
        if not self.vcoruc:

            rv_out = [99,99,99,99,99]
            bin_index = np.unique(tbin.values)

            for i in bin_index:
                rv_out[i] = np.mean(rvcor[tbin==i])

            for i in range(len(rv_out)):

                if rv_out[i] == 99 and i == 0:
                    for j in range(1,len(rv_out)):
                        if rv_out[i+j] != 99:
                            rv_out[i] = rv_out[i+j]
                            break

                elif rv_out[i] == 99:
                    for j in range(1,i+1):
                        if rv_out[i-j] != 99:
                            rv_out[i] = rv_out[i-j]
                            break

            self.v_cor = np.array(rv_out)


        self.abun = self.__calcu_abun()
        self.__deplete(self.deplete)

        
        ## A dangerous function to modify the elemental abundances.
        # for i,x in lines.groupby([ele,ion]):
        #     print((x.flux/x.pre_flux)[x.effcoe!=0].mean())
        #     print((x.flux/x.pre_flux)[x.effcoe!=0].std())
        #     print((x.flux/x.pre_flux)[x.effcoe!=0].median())

        #     if sum(x.effcoe!=0) == 0:
        #         continue

        #     elif sum(x.effcoe!=0) == 1:
        #         multiple = (x.flux/x.pre_flux)[x.effcoe!=0].item()/2

        #     else:  
        #         multiple = (x.flux/x.pre_flux)[x.effcoe!=0].median()
            
        #     self.abun[i[0]] = multiple*self.abun[i[0]]

 

        self._H_beta_flux(1)
        if min(self.obs_flux) < 0:
            self._H_beta_flux(0)
            

    def linesubframe(self,wavl,wavelerr,sigma):
        """
        Extract all candidate lines from the atomic transition database based on input
        wavelength, wavelength uncertainty and sigma values.
        """
        
        # All lines' wavelengths in the atomic transition database
        wav = self.linedb[:,0]

        # The basic type of sigma that can be input
        if isinstance(sigma,(int,float)):
            upperl = wavl+wavelerr[1]*wavl*sigma/self.c
            lowerl = wavl+wavelerr[0]*wavl*sigma/self.c
            linesubframe = self.linedb[(wav>=lowerl)&(wav<=upperl)]

        # An internal function used to find the possible strongest lines between the input sigma \
        # and double sigma. e.g., input sigma is 5, this will search between 5-10.
        elif isinstance(sigma,(list,np.ndarray,tuple)) and len(sigma) == 2:
            lowerl1 = wavl+wavelerr[0]*wavl*sigma[1]/self.c
            lowerl2 = wavl+wavelerr[0]*wavl*sigma[0]/self.c
            upperl1 = wavl+wavelerr[1]*wavl*sigma[0]/self.c
            upperl2 = wavl+wavelerr[1]*wavl*sigma[1]/self.c
            linesubframe = self.linedb[((wav>=lowerl1)&(wav<lowerl2))| \
                                       ((wav>upperl1)&(wav<=upperl2))]

        else:
            print("\nInvalid type of 'sigma'.")    
            sys.exit()
        
        if len(linesubframe) == 0:
            print('\nWavelength error is too low!')
            sys.exit()

        return linesubframe



    def calcu_wavlscore(self,linesubframe,wavl,wavlerr):
        """
        Calculate the score of "Wavelength Agreement" part.
        """

        self.wavscore = np.zeros(linesubframe.shape[0])
        ele_binindex = np.zeros(linesubframe.shape[0]).astype(int)
        permitted = linesubframe[:,3]==0
        forbidden = linesubframe[:,3]!=0
        ele = linesubframe[:,1].astype(int)
        ionnum = linesubframe[:,2].astype(int)

        # Locate which energy bin each candidate line is in
        ele_binindex[permitted] = self.ele_binindex[ele[permitted]-1,ionnum[permitted]]
        ele_binindex[forbidden] = self.ele_binindex[ele[forbidden]-1,ionnum[forbidden]-1]

        # For different candidate line in different energy bin, correct with different velocities.
        self.wav_cor = wavl/(1+(self.v_cor[ele_binindex]/self.c))

        # Calculate the differences between each candidate line and observed line
        self.wavdiff = (self.wav_cor-linesubframe[:,0])*self.c/linesubframe[:,0]

        # Score each candidate line based on the wavelength differences
        self.wavscore[self.wavdiff< 0] = (-self.wavdiff[self.wavdiff< 0]/(wavlerr[1])).astype(int)
        self.wavscore[self.wavdiff>=0] = (-self.wavdiff[self.wavdiff>=0]/(wavlerr[0])).astype(int)
        # wavscore[wavscore>=4] = 4

        

        return self.wavscore



    def calcu_flux(self,linesubframe,RR,DR,obs_flux):
        """
        Calculate the flux or equivalent width of each candidate line.
        """

        # If the observed line is an emission line
        if obs_flux >= 0:
            Aki = linesubframe[:,10].copy()
            cc = np.ones(linesubframe.shape[0])
            spec = linesubframe[:,3]
            dd = np.ones(linesubframe.shape[0])
            rr = np.ones(linesubframe.shape[0])
            upene = linesubframe[:,9]
            lowene = linesubframe[:,8]
            jk = linesubframe[:,5]
            ji = linesubframe[:,4]
            ele = linesubframe[:,1].astype(int) - 1
            ionnum = linesubframe[:,2].astype(int) - 1
            br = linesubframe[:,7]
            eff_ix = linesubframe[:,15].astype(int)

            # set the transition probability for those lines without references
            Aki[(linesubframe[:,3]==0)&(linesubframe[:,10]==0)] = 1e4
            Aki[(linesubframe[:,3]==1)&(linesubframe[:,10]==0)] = 10
            Aki[(linesubframe[:,3]==2)&(linesubframe[:,10]==0)] = 1e-5

            # Correction factor for the intercombination lines
            cc[(linesubframe[:,3]==1)&(lowene>3000)] = 3e-4
            cc[(linesubframe[:,3]==1)&(lowene<=3000)] = 1


            # If the ion is neutral de-value the collisional term, as neutral ions have coulomb \
            # repulsion working against excitation
            dd[linesubframe[:,2]!=1] = (linesubframe[:,2][linesubframe[:,2]!=1]-1)
            dd[(linesubframe[:,3]==0)] = 0
            # dd[(linesubframe[:,3]==0)&(linesubframe[:,10]==0)] = 0

            # For forbidden line, set the recombination term to 0
            rr[(linesubframe[:,3]==2)] = 0

            br[(linesubframe[:,1]==1)] = br[(linesubframe[:,1]==1)]*100


            RRs = np.tile([0,0,1,1,0,0],(linesubframe.shape[0],1)).astype(np.float64)
            pc = np.zeros((linesubframe.shape[0],9)).astype(np.float64)
            pE = np.zeros((linesubframe.shape[0],9)).astype(np.float64)

            Rcoe(RRs,pc,pE,linesubframe,RR,DR,rr)

            flux = self.em_flux_formu(self.abun[ele,ionnum],self.abun[ele,ionnum+1],cc,\
                upene,ji,jk,Aki,RRs,pc,pE,linesubframe[:,0],linesubframe[:,1],br,dd,rr)/self.H_beta_eff


        
            if sum(eff_ix) != 0:
                effcoe_cal(flux, eff_ix, ele, ionnum, self.abun, self.effcoes, self.effcoe_num, \
                        self.effcoe_ix, self.H_beta_eff)

        # If the observed line is an absorption line
        else:
            Aki = linesubframe[:,10].copy()
            Aki[(linesubframe[:,3]==0)&(linesubframe[:,10]==0)] = 1e4
            spec = linesubframe[:,3]
            flux = np.zeros_like(spec)

            # Set the flux of intercombination lines and collisionally excited lines to zero.
            ele = linesubframe[:,1][spec==0].astype(int) - 1
            ionnum = linesubframe[:,2][spec==0].astype(int) - 1
            Aki = Aki[spec==0]
            jk = linesubframe[:,5][spec==0]
            lowene = linesubframe[:,8][spec==0]
            wavl = linesubframe[:,0][spec==0]
            flux[spec==0] = self.abs_flux_formu(self.abun[ele,ionnum],Aki,lowene,jk,wavl)/self.H_abs
        
        return flux


    def em_flux_formu(self,abun1,abun2,cc,upene,ji,jk,Aki,RRs,pc,pE,wavl,Z,br,dd,rr):
        """
        Calculate the predicted flux of emission lines.

        Notes
        -----
        abun1 : the ion density of collisional excitation term.
        abun2 : the ion density of recombination term.
        cc : correction factor for intercombination lines in collisional excitation term.
        upene : the upper level energy.
        ji : statistic weight of the lower level.
        jk : statistic weight of the upper level.
        Aki : transition probability.
        RRs : array of parameters of radiative recombination coefficients.
        pc & pE : array of parameters of dielectronic recombination coefficients.
        wavl : wavelength in Angstrom.
        Z : ionization state from which line arises.
        br : factor that convert the total recombination coefficient to effective recombination 
             coefficient.
        dd : correction factor for neutral ions.
        """

        # The collisional excitation term
        q21 = 8.63e-6/((self.Te**0.5)*jk)
        q12 = 8.63e-6*np.exp(-upene*1240/(self.k*self.Te*1e7))/(((self.Te)**0.5)*ji)
        new_col = self.col_cor*cc*dd*abun1*q12/(1+self.Ne*q21/Aki)

        # The radiative recombination coefficient
        alpharr = RRs[:,0]*(np.sqrt(self.Te/RRs[:,2])*((1+np.sqrt(self.Te/RRs[:,2]))**(1-RRs[:,1]-\
                  RRs[:,4]*np.exp(-RRs[:,5]/self.Te)))*((1+np.sqrt(self.Te/RRs[:,3]))**(1+RRs[:,1]+\
                  RRs[:,4]*np.exp(-RRs[:,5]/self.Te))))**-1

        # The older general formula for the calculation of radiative recombination coefficient. \
        # This formula is used when the above one isn't available.
        alpharr[(RRs[:,0]==0)&(rr!=0)] = 1e-13*(Z[(RRs[:,0]==0)&(rr!=0)]+1)*\
                (self.Te/10000/(Z[(RRs[:,0]==0)&(rr!=0)]+1))**(-0.7)

        # The dielectronic recombination coefficient
        alphadr = np.sum(pc*np.exp(-pE/self.Te),axis=1)/(self.Te**(3/2))

        # The total recombination term
        new_rec = abun2*(alphadr+alpharr)*br

        return (new_col + new_rec)*1e14/wavl



    def abs_flux_formu(self,abun,Aki,lowene,jk,wavl):
        """
        Calculate the equivalent width of absorption lines. See `Line_list.em_flux_formu` for 
        details of each parameter.
        """
        Ni = abun*np.exp(-lowene*1240*1e-7/(self.k*self.Te))
        Wik = Ni*jk*(wavl)**2*Aki

        return Wik



    def calcu_fluxscore(self,linesubframe,obs_flux):
        """
        Calculate the score of "Predicted Template Flux" part.
        """
        # Calculate fluxes of candidate lines
        self.flux = self.calcu_flux(linesubframe,self.RR,self.DR,obs_flux)

        fluxscore = np.zeros(linesubframe.shape[0])

        # Find the maximum flux within the candidate lines
        maxf = self.flux.max()

        # For those lines with extremely low predicted fluxes, set to a certain value
        self.flux[self.flux<=1e-20] = 1e-24

        # Score each candidate line
        fluxscore = np.int64((np.log(maxf/self.flux)/np.log(10)))

        # If the highest flux has a reference of effective recombination coefficient,
        # increase its weight by subtracting 1
        # if linesubframe[np.argmax(self.flux),15] != 0:
        #     fluxscore[np.argmax(self.flux)] = fluxscore[np.argmax(self.flux)] - 1
        # For those lines with fluxes lower than 1e-5*max, all are set to 5 and ready \
        # to remove.
        fluxscore[fluxscore>=5] = 5

        return fluxscore



    def calcu_mulscore(self,wavl,pobs_flux,num):
        """
        Calculate the score of "Multiplet Check" part.
        """

        length = self.lineframe.shape[0]

        # Initialize the number of possibly observable multiplet members.
        possi_num = np.zeros(length).astype(int)

        # Initialize the number of detected multiplet members.
        detect_num = np.zeros(length).astype(int)

        # Initialize the score of "Multiplet Check".
        mulscore = np.zeros(length)

        # the notation to indicate whether all multiplet lines are within the parameter `I`. 
        self.wav_nota2 = np.zeros(length,dtype='<U4')
        self.wav_nota2[:] = '|'

        self.score_nota = np.zeros_like(mulscore,dtype='<U1')
        
        # save the wavelengths and wavelength differences of the multiplet lines
        self.candidate = np.zeros(length,dtype='<U50')

        # Boundaries on both sides of the wavelength of all input lines
        lowerls = self.wav*(self.c+self.waverr[:,0]*self.sigma)/self.c
        upperls = self.wav*(self.c+self.waverr[:,1]*self.sigma)/self.c

        if self.extract_done:
            Aterm = np.concatenate(\
                    (self.Aterm[:sum(self.Anum[:num])],self.Aterm[sum(self.Anum[:num+1]):]))
            Awave = np.concatenate(\
                    (self.Awave[:sum(self.Anum[:num])],self.Awave[sum(self.Anum[:num+1]):]))
        else:
            Aterm = np.zeros(1)
            Awave = np.zeros(1)
        # The main function used to search for possible detected multiplet lines from \
        # the input line list
        mul_check(self.wav_nota2,self.linedb,self.c,self.I,self.k,self.Te,self.Ne,\
            self.lineframe,self.wavdiff,self.wav_cor,self.wav,wavl,self.ele_binindex,\
            self.v_cor,self.waverr,self.obs_flux,lowerls,upperls,pobs_flux,\
            self.candidate,detect_num,possi_num,self.extract_done,Aterm,\
            Awave,self.score_nota,self.effcoes,self.effcoe_ix,self.effcoe_num)

        self.possi_num = possi_num
        self.detect_num = detect_num
        
        # Criteria of "Multiplet Check" part
        mulscore[((detect_num==2)&(possi_num==2))| \
                 (detect_num>2)] = 0

        mulscore[((detect_num==1)&(possi_num==1))| \
                 ((possi_num>2)&(detect_num==2))] = 1

        mulscore[((detect_num==0)&(possi_num==0))| \
                 ((possi_num>1)&(detect_num==1))] = 2

        mulscore[(detect_num==0)&(possi_num<=3)&(possi_num>0)] = 3

        mulscore[(detect_num==0)&(possi_num>3)] = 4

        # If it's 0/0 with the highest predicted flux, increse weight
        hf00 = (detect_num==0)&(possi_num==0)&(self.flux==max(self.flux))
        mulscore[hf00] = mulscore[hf00] - 1 

        # If the flux of observed line is lower than 1e-4 H_beta, do not score the \
        # rest conditions 
        if self.i != 0 and pobs_flux < 1e-4 and pobs_flux > 0:

            mulscore[(detect_num==0)&(possi_num>2)] = 3

        self.score_nota[self.score_nota!='^'] = ''
        # If it is the last iteration
        # if self.extract_done:
        #     self.score_nota = np.zeros_like(mulscore,dtype='<U1')
        #     # Remove the multiplet line itself
        #     Aterm = np.concatenate(\
        #             (self.Aterm[:sum(self.Anum[:num])],self.Aterm[sum(self.Anum[:num+1]):]))
        #     Awave = np.concatenate(\
        #             (self.Awave[:sum(self.Anum[:num])],self.Awave[sum(self.Anum[:num+1]):]))
            
        #     # Chech whether other multiplet lines are in ranking A list
        #     cond = np.isin(self.lineframe[:,12],Aterm)

        #     self.score_nota[cond] = '*'   
        #     self.score_nota[~cond] = ' ' 

        return mulscore



    def calcu_tatolscore(self,wavl,wavlerr,obs_flux,num,complete=True):
        """
        Calculate the total score of each candidate line.
        """

        # Generate the line subframe
        self.lineframe = self.linesubframe(wavl,wavlerr,self.sigma)

        # If it needs complete identification
        if complete:

            # Generate the peripheral line subframe
            exlineframe = self.linesubframe(wavl,wavlerr,[self.sigma,2*self.sigma])
            exflux = self.calcu_flux(exlineframe,self.RR,self.DR,obs_flux)
            exlineframe = exlineframe[np.argmax(exflux)]
            exflux = exflux[np.argmax(exflux)]

            fluxscore = self.calcu_fluxscore(self.lineframe,obs_flux)

            # Remove the candidate lines with score 5 of flux 
            self.lineframe = np.vstack((self.lineframe[fluxscore != 5],exlineframe))
            self.flux = np.append(self.flux[fluxscore != 5],exflux)
            fluxscore = np.append(fluxscore[fluxscore != 5],0)
            # Sort the candidate lines by template flux. Here is used for keep the largest 
            # template flux of duplicate lines which will be merged.
            self.lineframe[:-1] = self.lineframe[:-1][self.flux[:-1].argsort()[::-1]]
            fluxscore[:-1] = fluxscore[:-1][self.flux[:-1].argsort()[::-1]]
            self.flux[:-1] = np.sort(self.flux[:-1])[::-1] 

        
        else:
            fluxscore = self.calcu_fluxscore(self.lineframe,obs_flux)
            self.lineframe = self.lineframe[fluxscore != 5]
            self.flux = self.flux[fluxscore != 5]
            fluxscore = fluxscore[fluxscore != 5]


        # Initialize self.wavdiff
        _ = self.calcu_wavlscore(self.lineframe,wavl,wavlerr)

        mulscore = self.calcu_mulscore(wavl,obs_flux,num)
        wavscore = self.calcu_wavlscore(self.lineframe,wavl,wavlerr)

        # Calculate the total score
        totalscore = mulscore + fluxscore + wavscore

        if self.i == self.loop - 1:
            totalscore[-1] = 99

        totalscore[totalscore<0] = 0

        return totalscore       



    def _ion_notation(self,ele,ionnum,trans):
        """
        Generate the string of ion.
        """

        notation = self.elesym[ele-1] + ' ' + self.ionstat[ionnum-1]

        # `trans` 1 means intercombination line
        if trans == 1:
            notation = notation + ']'

        # `trans` 2 means collisionally excited line
        elif trans == 2:
            notation = '[' + notation + ']'

        return notation


    def ion_notation(self,lineframe):
        """
        Main process of Generating string of ion.
        """

        if lineframe.ndim == 2:
            ele = lineframe[:,1].astype(int)
            ionnum = lineframe[:,2].astype(int)
            trans = lineframe[:,3].astype(int)
            output = []

            for i,j,k in zip(ele,ionnum,trans):
                output.append(self._ion_notation(i,j,k))

            return np.array(output,dtype=str)

        else:
            return self._ion_notation(int(lineframe[1]),int(lineframe[2]),int(lineframe[3]))


    def _icf_v_out(self,icf,v_cor):
        """
        Print the `icf` and `v_cor` parameters in output file.
        """

        b = ' '
        return '\nICF Values: Bin/%'+f'{b:10s}'+'Velocity Structure: Bin/Vel(km/s)\n'\
                f'ix 1:  {str(icf[0]*100):6.6s}'+f'{b:14s}'+f'irvcor 1:  {v_cor[0]:.3f}\n'+\
                f'ix 2:  {str(icf[1]*100):6.6s}'+f'{b:14s}'+f'irvcor 2:  {v_cor[1]:.3f}\n'+\
                f'ix 3:  {str(icf[2]*100):6.6s}'+f'{b:14s}'+f'irvcor 3:  {v_cor[2]:.3f}\n'+\
                f'ix 4:  {str(icf[3]*100):6.6s}'+f'{b:14s}'+f'irvcor 4:  {v_cor[3]:.3f}\n'+\
                f'ix 5:  {str(icf[4]*100):6.6s}'+f'{b:14s}'+f'irvcor 5:  {v_cor[4]:.3f}\n'



    def _firstrun(self):
        """
        Extract a sub line list with relatively robust line identifications and re-run the
        identification process.
        """

        print(self._icf_v_out(self.icf,self.v_cor))
        output = []
        
        with tqdm(total=len(self.wav)) as pbar:
            pbar.set_description('SubProcessing:')

            for line,lineerr,obs_flux,num in zip(self.wav,self.waverr,self.obs_flux,range(len(self.wav))):
                # if line == 4471.657:
                #     breakpoint()
                
                # Calculate scores of all candidate lines
                score = self.calcu_tatolscore(line,lineerr,obs_flux,num,complete=False)
                mmin = np.min(score)
                # Find the candidate lines with minimum score
                minframe = self.lineframe[score==mmin]
                ion_nota = self.ion_notation(minframe)
                ele,ion = ion_nota[0].split()
                pre_flux = self.flux[score==mmin]
                ionna = np.unique(np.stack((minframe[:,1],minframe[:,2]),axis=-1),axis=0)
                # Find the second minimum score
                if len(np.unique(score)) > 1:
                    submin = np.unique(score)[1]

                else:
                    submin = 9

                # If number of candidate lines with minimum score is 1, and \
                # the second minimum score minus the minimum score is greater than 1
                if sum(score==mmin) == 1 and submin - mmin >= 2:
                    output.append([line,minframe[:,0][0],ele,ion,obs_flux,pre_flux[0],minframe[:,15][0]])

                # If number of candidate lines with minimum score is not 1, but \
                # they are all from the same ion
                elif sum(score==mmin) != 1 and len(ionna) == 1:
                    # Remove the intercombination lines
                    if len(np.unique(ion_nota)) != 1:
                        spec = np.bincount(minframe[:,3].astype(int)).argmax()
                        pre_flux = pre_flux[minframe[:,3]==spec]
                        minframe = minframe[minframe[:,3]==spec]
                    # elecheck = self.lineframe[:,1]!=self.lineframe[np.argmin(score),1]
                    # ioncheck = self.lineframe[:,2]!=self.lineframe[np.argmin(score),2]
                    # The scores of other ions

                    lab_wav = sum(minframe[:,5]*minframe[:,0])/sum(minframe[:,5])
                    pre_flux = sum(minframe[:,5]*pre_flux)/sum(minframe[:,5])     
                    output.append([line,lab_wav,ele,ion,obs_flux,pre_flux,sum(minframe[:,15] != 0)])

                
                # If number of candidate lines with minimum score is not 1, but \
                # they are all in the same energy bin
                elif sum(score==mmin) != 1 and submin - mmin >= 2:
                    enrbin = np.array([self.energybin]*sum(score==mmin))
                    eleenr = self.ionene[np.int8(minframe[:,1])-1,np.int8(minframe[:,2])-1]
                    enediff =  enrbin - np.tile(eleenr,(5,1)).T

                    if len(np.unique(np.sum(enediff<0,axis=1))) == 1:
                        output.append([line,minframe[:,0][0],ele,ion,obs_flux,pre_flux[0],0])

                    else:
                        pbar.update(1)
                        continue
                        
                else:
                    pbar.update(1)
                    continue

                pbar.update(1)

        # Re-calculate the `icf` and `v_cor` parameters
        self.match_lines(match_table=output)

        

    def _extract_Alist(self):
        """
        Extract a list with all ranking A candidates to further check whether other multiplet lines 
        of the tested transition are ranked as A. 
        """
        with tqdm(total=len(self.wav)) as pbar:

            pbar.set_description('Extracting A ranking list')

            for wavl,wavlerr,obs_flux,num in zip(self.wav,self.waverr,self.obs_flux,range(len(self.wav))):

                score = self.calcu_tatolscore(wavl,wavlerr,obs_flux,num,complete=False)
                mmin = np.min(score)

                # Find the candidate lines with minimum score
                minframe = self.lineframe[(score==mmin)]

                # Save the identification code of the multiplet and the number of ranking A candidates
                self.Aterm.extend(minframe[:,12].tolist())
                self.Awave.extend(minframe[:,0].tolist())
                self.Anum.append(len(minframe))

                pbar.update(1)

            
            self.Aterm = np.array(self.Aterm).reshape(-1)
            self.Awave = np.array(self.Awave).reshape(-1)
            self.Anum = np.array(self.Anum)
  

            

    def write(self,f1,f2):
        """
        Read the results of line identifications and write to the output files.
        """

        with tqdm(total=len(self.wav)) as pbar:
            pbar.set_description('MainProcessing:')
            for line,lineerr,num,obs_flux in zip(self.wav,self.waverr,range(len(self.wav)),self.obs_flux):
                if line == 3218.2:
                    breakpoint()
                score = self.calcu_tatolscore(line,lineerr,obs_flux,num).astype(int)
                # self.lineframe[:-1] = self.lineframe[:-1][self.flux[:-1].argsort()[::-1]]
                # self.flux[:-1] = np.sort(self.flux[:-1])[::-1] 

                # Sort each line and give the index values, remove duplicate lines
                index = np.unique(np.stack((self.lineframe[:-1,0], \
                                            self.lineframe[:-1,12], \
                                                ),axis=-1),axis=0,return_index=True)[1]

                # Append the index of peripheral line 
                index = np.append(index,len(self.lineframe)-1)

                lineinfo =  f'Number {num+1} \tObserved line:\t{line:.5f} \t {obs_flux:.2E}\t'+ \
                            f'SNR:{np.around(self.snr[num],2)}\tFWHM:{np.around(self.fwhm[num],2)}\n'

                # Sort each column by the index and generate the output arrays
                score = score[index]
                self.score_nota = self.score_nota[index]
                length = len(score)
                self.possi_num = self.possi_num[index]
                self.detect_num = self.detect_num[index]
                self.lineframe = self.lineframe[index]
                self.wavscore = self.wavscore[index]
                ranklist = np.zeros(length,dtype='<U1')
                termlist = np.zeros(length,dtype='<U23')
                configlist = np.zeros(length,dtype='<U71')
                multinota = np.zeros(length,dtype='<U6')
                idi = np.zeros(length,dtype='<U10')
                wav_nota1 = np.zeros(length,dtype='<U1')
                flux_nota = np.zeros(length,dtype='<U1')
                scores = score.copy()
                # scores[scores<=2] = 0
                uniquescore = np.unique(scores)

                # Generate the strings
                term(self.lineframe,self.tab_termnum,self.tab_ele,self.tab_stat,self.tab_term,\
                     self.tab_conf,uniquescore,scores,self.possi_num,self.detect_num,self.elesymarr,\
                     self.ionstatarr,self.wavscore,termlist,configlist,ranklist,multinota,idi,\
                     wav_nota1,flux_nota,obs_flux)

                self.flux = self.flux[index]
                flux = ["{:.2E}".format(i) for i in self.flux]
                self.wav_cor = self.wav_cor[index]
                self.wav_nota2 =  np.array(self.wav_nota2)[index]
                self.wavdiff = self.wavdiff[index]
                self.candidate = self.candidate[index]

                # Combine the output arrays
                out = np.stack((wav_nota1, \
                                self.wav_cor, \
                                self.wav_nota2, \
                                self.lineframe[:,0], \
                                idi, \
                                termlist,\
                                flux_nota,\
                                flux, \
                                multinota, \
                                score, 
                                self.score_nota, \
                                ranklist, \
                                self.wavdiff, \
                                self.candidate,\
                                configlist,\
                                self.lineframe[:,4],\
                                self.lineframe[:,5],\
                                ),axis=-1)

                out = out[out[:,9].astype(int).argsort()]
                Aout = out[(out[:,11]=='A')|(out[:,10]=='^')]
                Astr = ''.join(f'{i[4]}  {i[3]:.9s}, ' for i in Aout)


                # Write the complete output file
                np.savetxt(f1,out,fmt='%-2s%-9.8s%-3s%-13.9s%-11s%-23s%-1s%-11.8s%-6s%-6s%-1s%-5s%-10.8s' \
                                     +'%-54s%-71s%-5s%-5s',header=lineinfo,footer='\n\n',comments='')

                # Write the brief output file 
                f2.write(f'{line:8.2f}  {obs_flux:9.2E}  |  {Astr} \n')
                pbar.update(1)




if __name__ == "__main__":
    hf22 = np.loadtxt('../test/Hf2-2_linelist.txt',skiprows=1)
    hf22_out = Line_list(wavelength=hf22[:,0],wavelength_error=10,flux=hf22[:,1],snr=hf22[:,2],fwhm=hf22[:,3])    
    hf22_out.identify('Hf2-2',abun_type='nebula')
    # J0608 = pd.read_table('../test/J0608_linelist.txt',delim_whitespace=True)
    # J0608_out = Line_list(wavelength=J0608.wave_cor.values,wavelength_error=30,flux=J0608.F.values,snr=J0608.snr.values,fwhm=J0608.fwhm.values)
    # J0608_out.identify('J0608_2',Te=30000,abun_type='../test/abun_WC.dat')
    # ic418 = np.loadtxt('../test/ic418_linelist.txt',skiprows=1)
    # ic418_out = Line_list(wavelength=ic418[:,0],wavelength_error=10,flux=ic418[:,3],snr=ic418[:,5],fwhm=ic418[:,4])
    # ic418_out.identify('ic418_2',abun_type='nebula')