import numpy as np
import pandas as pd
import os
from scipy import interpolate
import matplotlib.pyplot as plt
import emcee
import corner
import sys

# from scipy.optimize import least_squares

class Recom_Lines(object):

    def __init__(self, filename, eleabuns=None, flux_normalize=1):
        """
            Create a Recom_Lines object that contains the recombination lines that need 
            to be fitted with effective recombination coefficients.
        """
        self.filename = filename
        
        self.rootdir = os.path.abspath(os.path.join(os.path.dirname( \
                       os.path.abspath(__file__)), os.pardir))

        self.effcoe_ion = ['HI','HeI','NII','OII','CII']

        if eleabuns == None:
            self.eleabuns = np.array([1,0.85E-01,0.68E-04,0.49E-03,0.27E-03])
        elif isinstance(eleabuns,(list,np.ndarray)) and len(eleabuns) == 5:
            self.eleabuns = eleabuns
        else:
            raise AttributeError("Wrong type of `eleabuns`.")
        
        self.coeTe = np.load(os.path.join(self.rootdir,'pyemili','eff_reccoe',\
                                'eff_Te.npy'),allow_pickle=True).item()
        self.coeNe = np.load(os.path.join(self.rootdir,'pyemili','eff_reccoe',\
                                'eff_Ne.npy'),allow_pickle=True).item()

        self.Coes = {}
        self.Lines = {}
        for i in self.effcoe_ion:

            self.Coes[f'{i}'] = np.load(os.path.join(self.rootdir,'pyemili','eff_reccoe',\
                                    f'{i}_emiss_data.npy'),allow_pickle=True).astype(np.float64)

            self.Lines[f'{i}'] = np.loadtxt(os.path.join(self.rootdir,'pyemili','eff_reccoe',\
                                    f'{i}_emiss_lines.dat'))
            
        self.c = 2.9979246*1.0E+5

        self.data = pd.read_table(filename,delim_whitespace=True, \
                    names=['obs_wav','obs_flux','obs_fluxerr','ele','ion','lab_wav'],comment='#')

        self.data = self.data[self.data.obs_flux!=self.data.obs_fluxerr]
        self.data = self.data.sort_values('obs_flux')
        self.ion = self.data.ion.apply(lambda x: x if x[-1] != ']' else x[:-1])

        if isinstance(flux_normalize,(float,int)) and flux_normalize != 1:
            self.data.obs_flux = self.data.obs_flux/flux_normalize
            self.data.obs_fluxerr = self.data.obs_fluxerr/flux_normalize



    def inter(self,length,Te,Ne,coe):
        """Create a interpolator for all the lines with specific spectrum in the line table.
        The interpolation is linear.
        """
        if length == 0:
            return interpolate.RegularGridInterpolator((Te,Ne),coe)
        
        return interpolate.RegularGridInterpolator((np.arange(length),Te,Ne),coe)


    #prior distribution
    def log_prior(self,Te,Ne,abun):

        if self.paras_valid[2]:
            if np.log10(self.eleabun/100) <= abun <= np.log10(self.eleabun*100) and \
                self.minTe <= Te <= self.maxTe and \
                self.minNe <= Ne <= self.maxNe:
                    return 0.0
        
        else:
            if  self.minTe <= Te <= self.maxTe and \
                self.minNe <= Ne <= self.maxNe:
                    return 0.0
            # if 3000 <= Te <= 4000 and \
            #    2000 <= Ne <=4000:
        
        return -np.inf

    def return_theta(self,theta):
        
        defval = self.default*-(self.paras_valid-1)

        defval[self.paras_valid==1] += theta

        return defval

    #likelihood function
    def log_probability(self,theta):

        Te,Ne,abun = self.return_theta(theta)

        lp = self.log_prior(Te,Ne,abun)
        if not np.isfinite(lp):
            return -np.inf
        

        # indexes of selected Te and Ne
        paras = np.stack((np.arange(len(self.coe)).reshape((-1,1)),\
                        np.tile(Te,(len(self.coe),1)),\
                        np.tile(Ne,(len(self.coe),1))),\
                        axis=-1).reshape((-1,3))
        
        if self.maxinter == self.Hbetacoe and not self.isHI:
            if np.log10(500) > Te:
                Te = np.log10(500)
            elif Te > np.log10(30000):
                Te = np.log10(30000)
            if Ne < 2:
                Ne = 2
            # Te = 8900
            # Ne = 10000

        model = 10**(abun)*self.intermec([paras])[0]/self.maxinter([Te,Ne])[0]

        obs = self.subdata.obs_flux.values/self.obsmax
        obserr2 = (self.subdata.obs_fluxerr.values/self.obsmax)**2 # + (self.flux_threshold)**2

        return  -0.5 * np.sum(((obs - model) ** 2 / obserr2) ) #+ np.log(2*np.pi*obserr2)

    
    def _get_modelvalue(self, Te,Ne,abun):
        Te = np.log10(Te)
        Ne = np.log10(Ne)
        
        paras = np.stack((np.arange(len(self.coe)).reshape((-1,1)),\
                          np.tile(Te,(len(self.coe),1)),\
                          np.tile(Ne,(len(self.coe),1))),axis=-1).reshape((-1,3))
        
        if self.maxinter == self.Hbetacoe and not self.isHI:
            if np.log10(500) > Te:
                Te = np.log10(500)
            elif Te > np.log10(30000):
                Te = np.log10(30000)
            if Ne < 2:
                Ne = 2
            # Te = 8900
            # Ne = 10000

        model = abun*self.intermec([paras])[0]/self.maxinter([Te,Ne])[0]

        return model
    

    def check_TeNe_range(self,Te,Ne,spectrum):

        if Te:
            if self.maxTe < self.Te or self.minTe > self.Te:
                print("Electron temperature out of limit range")
                print(f'{spectrum[0]} {spectrum[1]} Te range: {10**self.minTe:.1f}'+\
                        f'-{10**self.maxTe:.1f}')
                sys.exit()

        if Ne:
            if self.maxNe < self.Ne or self.minNe > self.Ne:
                print("Electron density out of limit range")
                print(f'{spectrum[0]} {spectrum[1]} Ne range: {10**self.minNe:.1f}'+\
                        f'-{10**self.maxNe:.1f}')
                sys.exit()


    def mcmc_run(self,Spectrum,Te=None,Ne=None,abun=True,flux_threshold=0,\
                 suptitle=True,plotscatter=False,set_limits=False,vel=20):#set_limits=False,rpd=False,save=False,

        self.paras_valid = np.zeros(3)
        self.flux_threshold = flux_threshold
        self.suptitle = suptitle

        if Te:
            self.Te = np.log10(Te)
        else:
            self.Te = 0
            self.paras_valid[0] = 1

        if Ne:
            self.Ne = np.log10(Ne)
        else:
            self.Ne = 0
            self.paras_valid[1] = 1

        
        if abun or (self.Te and self.Ne):
            self.paras_valid[2] = 1

        if sum(self.paras_valid) == 0:
            print("None of the parameter will be fitted.")
            sys.exit()

        
        self.isHI = False


        if not isinstance(Spectrum,list):
            exit("Spectrum must be a list")
        
        if not np.all([i.replace(' ','') in self.effcoe_ion for i in Spectrum]):
            exit("Unknown spectrum exists")
        
        
        self.Hbetacoe = self.inter(0,np.log10(self.coeTe['HI']),
                                   np.log10(self.coeNe['HI']),
                                    self.Coes['HI'][21,:,:])


        for i in range(len(Spectrum)):

            self.spectrum = Spectrum[i].split()
            self.spec = ''.join(self.spectrum)

            if self.spec == 'HI':
                self.isHI = True
                self.paras_valid[2] = 0
                
            self.minTe = min(np.log10(self.coeTe[self.spec]))
            self.maxTe = max(np.log10(self.coeTe[self.spec]))
            self.minNe = min(np.log10(self.coeNe[self.spec]))
            self.maxNe = min(max(np.log10(self.coeNe[self.spec])),5)
            
            self.check_TeNe_range(self.Te,self.Ne,self.spectrum)

            ix = self.effcoe_ion.index(self.spec)
            self.eleabun = self.eleabuns[ix]
            
            self.subdata = self.data[(self.data.ele==self.spectrum[0])&(self.ion==self.spectrum[1])]
            self.subdata = self.subdata[self.subdata.obs_flux>=flux_threshold]
            self.subdata = self.subdata.reset_index(drop=True)

            if len(self.subdata) == 0:
                print('Cannot find any lines for specific spectrum. Please check your parameters setting.')
                sys.exit()

            if not self.paras_valid[2] and not self.isHI:
                relele = ['OII','NII','HeI']
                rellab_wav = [4649.14,5679.56,4471.47]
                ercs = [6438,2992,373]
                lab_wav = rellab_wav[relele.index(self.spec)]
                erc = ercs[relele.index(self.spec)]
                imax = np.argmin(abs(self.subdata.lab_wav-lab_wav))
                self.obsmax = self.subdata.iloc[imax].obs_flux
                self.maxinter = self.inter(0,np.log10(self.coeTe[self.spec]),
                                           np.log10(self.coeNe[self.spec]),
                                           self.Coes[self.spec][erc,:,:])
            else:
                self.maxinter = self.Hbetacoe
                self.obsmax = 1
            # self.subdata = self.subdata[self.subdata.obs_fluxerr/self.subdata.obs_flux<=0.2]
            # self.subdata = self.subdata[self.subdata.obs_flux>=1e-3]



            flat_samples = self._mcmc_para(vel=vel)


            flat_samples = 10**(flat_samples)


            

            # if rpd:
            #     y =  self.subdata.obs_flux/self.obsmax
            #     yerr = self.subdata.obs_fluxerr/self.obsmax
            #     y_model = self._get_modelvalue(flat_samples)
            #     pre_len = len(self.subdata)
            #     self.subdata = self.subdata[abs(y_model-y)<=abs(y-yerr)]
            #     if len(self.subdata) != pre_len:
            #         self.subdata = self.subdata.reset_index(drop=True)
            #         flat_samples = self._mcmc_para()

            if set_limits:
                lower = np.percentile(flat_samples,16,axis=0)
                mid = np.percentile(flat_samples,50,axis=0)
                upper = np.percentile(flat_samples,84,axis=0)
                lower_l = mid - 3*(mid-lower)
                upper_l = mid + 3*(upper-mid)
                internal = (flat_samples>=lower_l)&(flat_samples<=upper_l)
                intersec = np.prod(internal,axis=1)
                flat_samples = flat_samples[intersec!=0]

            self.plot_corner(flat_samples,flux_threshold)
            
            if plotscatter:
                if self.paras_valid[2]:
                    flat_samples[:,-1] = flat_samples[:,-1]/10**-self.abunorder
                self.plot_fig(flat_samples)



    def _mcmc_para(self,vel):

        lab_wav = self.Lines[self.spec][:,0]
        self.coe = np.zeros((len(self.subdata),\
                            self.Coes[self.spec].shape[1],\
                            self.Coes[self.spec].shape[2]))
        
        for i,j in self.subdata.iterrows():
            low_lim = j.lab_wav*(1-vel/self.c)           
            up_lim = j.lab_wav*(1+vel/self.c)     
            self.coe[i] =  np.sum(self.Coes[self.spec][(lab_wav>=low_lim)&(lab_wav<=up_lim)],axis=0)



        # effcoe_ix = self.subdata['erc'].values - 1
        # self.coe = self.Coes[self.spec][effcoe_ix]
        self.intermec = self.inter(len(self.coe),
                                   np.log10(self.coeTe[self.spec]),
                                   np.log10(self.coeNe[self.spec]),self.coe)
        

        

        theta0 = np.array([3,3,np.log10(self.eleabun)])*self.paras_valid
        theta0 = theta0[theta0!=0]

        self.default = np.array([self.Te,self.Ne,-(self.paras_valid[2]-1)]).astype(np.float64)


        # res_lsq = least_squares(self.log_probability, [3,3,np.log10(self.eleabun)],bounds=([self.minTe,self.minNe,np.log10(self.eleabun/100)],[self.maxTe,self.maxNe,np.log10(self.eleabun*100)]),max_nfev=100000)
        # sigma0 = [300,300,0.49E-04]
        pos = theta0 + 0.1*np.random.randn(16, len(theta0))
        # pos = np.random.normal(loc=theta0,scale=sigma0,size=(20,3))
        paras_name = np.array(['Te','Ne','abun'])
        self.fitparalist = [i for i in paras_name[self.paras_valid==1]]
        print(f"Fitting parameters:{self.fitparalist}")
        deparaname = paras_name[self.paras_valid==0]
        self.deparaname = [f'{deparaname[i]}={self.default[self.paras_valid==0][i]}'for i in range(len(deparaname))]
        print(f"Specific parameters:{self.deparaname}")
        # Set up the backend
        # filename = "IC4776OII.h5"
        # backend = emcee.backends.HDFBackend(filename)
        # Don't forget to clear it in case the file already exists
        nwalkers, ndim = pos.shape       
        # backend.reset(nwalkers, ndim)
        sampler = emcee.EnsembleSampler(nwalkers, ndim, self.log_probability)#,backend=backend)
        max_n = 30000
        # We'll track how the average autocorrelation time estimate changes
        index = 0
        autocorr = np.empty(max_n)
        burnin = 0
        # This will be useful to testing convergence
        old_tau = np.inf
        # Now we'll sample for up to max_n steps
        for sample in sampler.sample(pos, iterations=max_n, progress=True):
            # Only check convergence every 100 steps
            if burnin == 0:
                if sampler.iteration % 500:
                    continue

                # Compute the autocorrelation time so far
                # Using tol=0 means that we'll always get an estimate even
                # if it isn't trustworthy
                tau = sampler.get_autocorr_time(tol=0)
                autocorr[index] = np.mean(tau)
                index += 1

            # Check convergence
                converged = np.all(tau * 100 < sampler.iteration)
                converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
                if converged:
                    # state = sample
                    burnin = sampler.iteration
                    # break
                    print('\nConvergence')
                old_tau = tau
            
            else:
                if sampler.iteration == burnin+1000:
                    break

        tau = sampler.get_autocorr_time()
        burnin = int(2 * np.max(tau))
        thin = int(0.5 * np.min(tau))
        flat_samples = sampler.get_chain(discard=burnin, flat=True, thin=thin)

        return flat_samples
    
    def plot_corner(self,flat_samples,flux_threshold):

        labels = []
        if self.paras_valid[0]:
            labels.append(r"$Te$")
        if self.paras_valid[1]:
            labels.append(r"$Ne$")
        if self.paras_valid[2]:
            self.abunorder = int(np.log10(np.percentile(flat_samples[:,-1],50))) - 1
            flat_samples[:,-1] = flat_samples[:,-1]*10**-self.abunorder
            labels.append(rf"$Abun (\times 10^{{{self.abunorder}}})$")

        figure = corner.corner(
        flat_samples,
        labels=labels,
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True,
        title_kwargs={"fontsize": 12},
        range=[1,(3000,20000),1]
        )
        if self.suptitle:
            self.suptitle = ','.join([i for i in self.deparaname if i.startswith('Te') or i.startswith('Ne')])
            self.suptitle += f' flux_threshold={flux_threshold:.1e}'
            figure.suptitle(self.suptitle,fontsize=8)
            figure.tight_layout()
        # plt.figure()
        figure.savefig(f"{self.filename[:-4]}_{self.spec}_{''.join(self.fitparalist)}.png")
        # if save:
        #     np.savetxt(f"{self.filename[:-4]}_{self.spec}_{''.join(self.fitparalist)}_samples.txt",flat_samples)


    def plot_fig(self,flat_samples):

        plt.figure(figsize=(16,9))
        x = np.arange(len(self.coe))
        y =  self.subdata.obs_flux/self.obsmax
        plt.errorbar(x, y, yerr = self.subdata.obs_fluxerr/self.obsmax,fmt='o',capsize=7,markersize=5,label='Obs')
        # plt.scatter(x, y, 8, label='obs')
        for i,wav in enumerate(self.subdata.lab_wav.values):
            plt.annotate(f'{wav:.2f}', (x[i], y[i]),rotation=-45,rotation_mode='anchor')
        Te,Ne,abun = self.return_theta(np.percentile(flat_samples,50,axis=0))
        plt.plot(np.arange(len(self.coe)),self._get_modelvalue(Te,Ne,abun),'*',markersize=12,label='Model')
        
        plt.legend()
        plt.xlabel('Number (n)')
        plt.ylabel('Relative intensity')
        if self.suptitle:
            plt.suptitle(self.suptitle,fontsize=12)

        plt.savefig(f"{self.filename[:-4]}_{self.spec}_{''.join(self.fitparalist)}_scatter.png")





if __name__ == '__main__':

    ic = Recom_Lines('IC4776_v2_erc.dat')
    ic.mcmc_run(['O II'],Te=3560,Ne=2850,abun=True)   
