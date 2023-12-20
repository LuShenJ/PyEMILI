from pyemili.Lines import Line_list
from pyemili.Spec import Spec_line_finding
import numpy as np
import pandas as pd



# # run IC418 line identifications
# ic418 = np.loadtxt('ic418_linelist.txt',skiprows=1)
# ic418_out = Line_list(wavelength=ic418[:,0],wavelength_error=10,flux=ic418[:,3],snr=ic418[:,5],fwhm=ic418[:,4])
# ic418_out.identify('ic418',abun_type='nebula')


# # run Hf 2-2 line identifications
# hf22 = np.loadtxt('Hf2-2_linelist.txt',skiprows=1)
# hf22_out = Line_list(wavelength=hf22[:,0],wavelength_error=10,flux=hf22[:,1],snr=hf22[:,2],fwhm=hf22[:,3])
# hf22_out.identify('Hf2-2',abun_type='nebula')



# # run J0608 line identifications
# J0608 = pd.read_table('J0608_linelist.txt',delim_whitespace=True)
# J0608_out = Line_list(wavelength=J0608.wave_cor.values,wavelength_error=30,flux=J0608.F.values,snr=J0608.snr.values,fwhm=J0608.fwhm.values)
# J0608_out.identify('J0608',icf=[0.01,0.5,0.4,0.1,0.00001],Te=30000,abun_type='abun_WC.dat',deplete='O 50 3')

# linefind = Spec_line_finding('J0608_03may2019MagE.fits',length=50,percentile=30,snr_threshold=7,prominence=6)
# linefind = Spec_line_finding('Hf2_2_346_30m2as_0027_3614_3671_m_ss.txt',length=4,percentile=25,snr_threshold=7,prominence=7)
# linefind = Spec_line_finding('Hf2_2_346_30m2as_0030_3701_3759_m_ss.txt',length=4,percentile=25,snr_threshold=7,prominence=7)