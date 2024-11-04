""" All functions here use `numba` to accelerate computations,
    and are used by `Line_list` class. 
"""
import numpy as np 
import numba as nb



# Generate each recombination coefficient of candidate line from the 'RR_rate2.dat' \
# and 'DR_rate2.dat' files
@nb.njit(parallel=True,cache=True)
def Rcoe(RRs,pc,pE,linesubframe,RR,DR,rr):

    for i in nb.prange(len(linesubframe)):

        # Check if it's forbidden line
        if rr[i] == 0:
            continue
        # Match the certain ion
        cond1 = (RR[:,0]==linesubframe[i,1])&(RR[:,1]==linesubframe[i,2])

        # Check if the parameters of this certain ion exists
        if sum(cond1) == 0:
            pass

        # Save the parameters of radiative recombination coefficients
        else:
            RRs[i] = RR[cond1][0][2:8]

        # Do the same thing to dielectronic recombination coefficients
        cond2 = (DR[:,0]==linesubframe[i,1])&(DR[:,1]==linesubframe[i,2])

        if sum(cond2) == 0:
            pass

        else:
            pc[i] = DR[cond2][0][2:11]
            pE[i] = DR[cond2][0][11:20]

            
@nb.njit(parallel=True,cache=True)
def effcoe_cal(flux, eff_ix, ele, ionnum, abun,  effcoes, effcoe_num, \
              effcoe_ix, H_beta_eff):
    
    cond = (eff_ix != 0)
    for i in nb.prange(sum(cond)):

        for j in range(len(effcoe_num)):

            if ele[cond][i] == effcoe_num[j,0] and ionnum[cond][i] == effcoe_num[j,1]:

                eff = effcoes[j]
                ix = effcoe_ix[j]
                break
            else:
                continue

        flux[np.argwhere(cond)[i]] = abun[ele[cond][i],ionnum[cond][i]+1]*eff[eff_ix[cond][i]-1,ix[0],ix[1]]/ \
                6.626e-27/2.99e18*1e14/H_beta_eff

        

# The main process foe checking the multiplet lines.
@nb.njit(parallel=True,cache=True)
def mul_check(wav_nota2, linedb, c, I, k, Te, Ne, lineframe, \
    wavdiff, wav_cor, wav, wavl, ele_binindex, v_cor, \
    waverr, obs_flux, lowerls, upperls, pobs_flux, candidate, \
    detect_num, possi_num, extract_done, Aele, Aion, Alowterm, \
    Aupterm, Awave, Awavdiff, score_nota, effcoes, effcoe_ix, effcoe_num):

    for i in nb.prange(len(lineframe)):

        # Match all multiplet lines from the atomic transition database
        multiplet = linedb[(linedb[:,1]==lineframe[i,1])&(linedb[:,2]==lineframe[i,2])&\
                           (linedb[:,12]==lineframe[i,12])&(linedb[:,13]==lineframe[i,13])]
        # multiplet = multiplet[(multiplet[:,12]==lineframe[i,12])&(multiplet[:,13]==lineframe[i,13])]
        # if lineframe[i,0] == 5346.02:
        #     breakpoint()

        # If just only one multiplet line exist
        if len(multiplet) == 1:
            multi = multiplet[0]
            cond = 0
            # If it's H I, check if Alist has H I lines from the same lower level but higher upper level
            if multi[1] == 1 and multi[2] == 1:
                cond = sum((multi[1]==Aele)&(multi[2]==Aion)&(multi[12]==Alowterm)&(multi[13]<Aupterm)) > 0

            # If it's He II, check if Alist has He II lines from the same lower level but higher upper level
            if multi[1] == 2 and multi[2] == 2:
                cond = sum((multi[1]==Aele)&(multi[2]==Aion)&(multi[12]==Alowterm)&(multi[13]<Aupterm)) > 0
            
                    
            if cond:
                # Average wave_diff of the same ions with the same lower level
                Adiff = np.mean(Awavdiff[(multi[1]==Aele)&(multi[2]==Aion)&(multi[12]==Alowterm)])
                #v_cor for recombination line
                vcor_index = ele_binindex[int(multi[1])-1,int(multi[2])]
                lamcor = wavl/(1+(v_cor[vcor_index])/c)
                wvdiv = c*(lamcor-multi[0])/multi[0]
                if Adiff >= 0:
                    usevel = waverr[np.argmin(np.abs(wav-multi[0])),0]
                else:
                    usevel = waverr[np.argmin(np.abs(wav-multi[0])),1]

                if np.abs(wvdiv-Adiff) <= np.abs(usevel):
                    score_nota[i] = '^'

            continue

        # If the wavelengths of all multiplet lines are within the input parameter `I`
        elif (multiplet[-1,0]-multiplet[0,0])*c/multiplet[0,0] <= I:
            # Calculate weighted wavelengths
            lineframe[i,0] = sum(multiplet[:,0]*multiplet[:,5])/sum(multiplet[:,5])
            # Calculate the new wavelength difference
            wavdiff[i] = (wav_cor[i]-lineframe[i,0])*c/lineframe[i,0]
            # Mark an `*` in the output file
            wav_nota2[i] = '| *'
        
        else:
            mul_wav = multiplet[:,0]
            # leftbd = wavl*(c+wavlerr[0]*sigma)/c
            # rightbd = wavl*(c+wavlerr[1]*sigma)/c
            # internel = (mul_wav>=leftbd)&(mul_wav<=rightbd)

            # Determine the wavelength bounds
            wav_range = ((mul_wav>=np.min(wav))&(mul_wav<=np.max(wav)))
            # Match the multiplet lines in appropriate wavelength range and transition types
            multiplet = multiplet[wav_range& \
                                # ~internel
                                #  (multiplet[:,4]>=lineframe[i,4]-2)& \
                                #  (multiplet[:,5]>=lineframe[i,5]-2)& \
                                 (multiplet[:,11]<=lineframe[i,11])]

            # internel = np.argwhere((lineframe[i,0]>lowerls)&(lineframe[i,0]<upperls)==True)
            # if len(internel) !=0:

            # Determine the index of the observed line
            internel = np.argwhere(wav==wavl)[0][0]
            # else:
            #     internel = -1

            # Remove the candidate line itself in multiplet lines
            multiplet = multiplet[multiplet[:,0]!=lineframe[i,0]]
            candimul = ''

            if len(multiplet) != 0:

                mul_rangen = np.zeros(len(multiplet))

                for z in nb.prange(len(multiplet)):
                    # Match each multiplet line, check whether it is in one wavelength range of all input \
                    # observed lines (except the line that is being tested)
                    condi = (multiplet[z,0]>lowerls)&(multiplet[z,0]<upperls)

                    # If one wavelength range is matched
                    if sum(condi) == 1:

                        # Check wether it is in the range of the line being tested
                        if np.argwhere(condi==True)[0][0] == internel:
                            mul_rangen[z] = -2
                        else:
                            mul_rangen[z] = np.argwhere(condi==True)[0][0]
                    
                    # If more than one wavelength range is matched
                    elif sum(condi) > 1:

                        if np.argmin(np.abs(wav-multiplet[z,0])) == internel:
                            mul_rangen[z] = -2
                        
                        # Match it to the closest one
                        else:
                            mul_rangen[z] = np.argmin(np.abs(wav-multiplet[z,0]))
                    else:
                        mul_rangen[z] = -1

                # Organize the array of index
                bins = np.unique(mul_rangen[(mul_rangen!=-1)&(mul_rangen!=-2)])
                possi_num[i] = len(bins)+sum(mul_rangen==-1)
                # possinum = 0
                detenum = 0

                for j in nb.prange(len(bins)):
                    condi = mul_rangen==bins[j]
                    multis = multiplet[condi]

                    for multi in multis:
                        wwav = wav[int(bins[j])]
                        wobs_flux = obs_flux[int(bins[j])]

                        if (pobs_flux > 0 and wobs_flux > 0) or (pobs_flux < 0 and wobs_flux < 0):
                            # possinum += 1
                            fluxratio = pobs_flux/wobs_flux
                            
                            # Check whether effective recombination coefficient exists
                            if multi[14] != 0 and lineframe[i,14] != 0:
                                cond = (effcoe_num[:,0]==multi[1]-1)&(effcoe_num[:,1]==multi[2]-1)
                                # match the ion
                                effix = effcoe_ix[cond][0]
                                # match the index
                                ix = np.argwhere(cond)[0][0]
                                # Calculate the ratio of the effective recombination coefficient
                                ajkrat = effcoes[ix][int(lineframe[i,14])-1][effix[0]][effix[1]]/ \
                                         effcoes[ix][int(multi[14])-1][effix[0]][effix[1]]

                                # Check the real flux ratio and the theoretical ratio
                                if fluxratio > 10*ajkrat or fluxratio < (0.1)*ajkrat:
                                    continue

                            # Check whether the reference of transition probability exists
                            elif    lineframe[i,10] != 0 \
                                    and multi[10] != 0 \
                                    and lineframe[i,5] <= 99 \
                                    and lineframe[i,11] == multi[11]:

                                # Calculate a theoretical ratio
                                if lineframe[i,3] == 2:
                                    q12_l = 8.63e-6*np.exp(-lineframe[i,9]*1240/(k*Te*1e7))/(((Te)**0.5)*lineframe[i,4])
                                    q12_m = 8.63e-6*np.exp(-multi[9]*1240/(k*Te*1e7))/(((Te)**0.5)*multi[4])
                                    q21_l = 8.63e-6/((Te**0.5)*lineframe[i,5])
                                    q21_m = 8.63e-6/((Te**0.5)*multi[5])
                                    ajkrat = q12_l/(1+Ne*q21_l/lineframe[i,10])/(q12_m/(1+Ne*q21_m/multi[10]))
                                else:
                                    ajkrat = lineframe[i,10]*lineframe[i,5]/(multi[10]*multi[5])

                                # Check the real flux ratio and the theoretical ratio
                                if fluxratio > 10*ajkrat or fluxratio < (0.1)*ajkrat:
                                    continue
                            
                            # If these two multiplet lines are from different transition types, \
                            # relax the restriction
                            elif lineframe[i,11]>multi[11]:

                                if fluxratio > 10**4:
                                    continue
                                
                            else: 
                                if fluxratio < 0.1 or fluxratio > 10:
                                    continue
                        # If these two lines are not the same type (em. or abs.)
                        elif wobs_flux >= 0 and pobs_flux < 0:
                            continue
                        
                        elif wobs_flux < 0 and pobs_flux >= 0:
                            continue

                        else:
                            # possinum += 1
                            pass

                        # Check the agreement of wavelength differences
                        if multi[3] == 0:
                            vcor_index = ele_binindex[int(multi[1])-1,int(multi[2])]

                        else:
                            vcor_index = ele_binindex[int(multi[1])-1,int(multi[2])-1]

                        lamcor = wwav/(1+(v_cor[vcor_index])/c)
                        wvdiv = c*(lamcor-multi[0])/multi[0]

                        if wavdiff[i]>=0:
                            usevel = waverr[int(bins[j]),0]

                        else:
                            usevel = waverr[int(bins[j]),1]

                        if np.abs(wvdiv-wavdiff[i]) > np.abs(usevel):
                            continue

                        # Save the detected multiplet lines
                        if detenum <= 2 and extract_done:
                            # check for Alist if there is any multiplet of this candidate
                            if sum((multi[0]==Awave)&(multi[1]==Aele)&(multi[2]==Aion)&\
                                   (multi[12]==Alowterm)&(multi[13]==Aupterm)) == 1:
                                ast = '^'
                                score_nota[i] = '^'

                            else:
                                ast = ''

                            wi = int(multi[0])
                            ws = int((multi[0]-wi)*1000)
                            di = int(wvdiv)
                            ds = int(abs(wvdiv-di)*10)
                            if ws < 100:
                                candimul += f'{ast}{wi}.0{ws}  {di}.{ds} |'
                            else:
                                candimul += f'{ast}{wi}.{ws}  {di}.{ds} |'
                        detenum += 1

                        break

                if detenum != 0:
                    candidate[i] = candimul
                detect_num[i] = detenum
                # possi_num[i] = possinum



# The process for generate the strings in output files.
@nb.njit(parallel=True,cache=True)
def term(lineframe,ion_tnum,ion_ele,ion_stat,ion_term,ion_config,\
    score,scores,possi_num,detect_num,elesym,ionstat,wavscore,\
    termlist,configlist,ranklist,multinota,idi,wav_nota1,flux_nota,\
    obs_flux):

    for i in nb.prange(len(lineframe)):
        cond1 = (lineframe[i,1] == ion_ele) & (lineframe[i,2] == ion_stat)

        ion_term1 = ion_term[cond1]
        ion_config1 = ion_config[cond1]
        ion_tnum1 = ion_tnum[cond1]

        condlowt = lineframe[i,12] == ion_tnum1
        condupt = lineframe[i,13] == ion_tnum1

        # Transition terms
        termlist[i] = ion_term1[condlowt][0]+'-'+ion_term1[condupt][0]

        # Electron configurations
        configlist[i] = ion_config1[condlowt][0]+'-'+ion_config1[condupt][0]

        # Rankings
        ranklist[i] = 'A' if scores[i]==score[0] \
                 else 'B' if scores[i]==score[1] \
                 else 'C' if scores[i]==score[2] \
                 else 'D' if scores[i]==score[3] \
                 else ''      

        # Notations for "Multiplet Check"
        multinota[i] = f'{possi_num[i]}/{detect_num[i]}'

        # Notations for ions
        notation = elesym[int(lineframe[i,1])-1] + ' ' + ionstat[int(lineframe[i,2])-1]

        if lineframe[i,3] == 1:
            idi[i] = notation + ']'
        elif lineframe[i,3] == 2:
            idi[i] = '[' + notation + ']'
        else:
            idi[i] = notation

        # Notations for wavelength differences
        wav_nota1[i] = '+' if wavscore[i] < 1 else ''

        # Notations for whether reference of transition probability exists
        flux_nota[i] = '*' if lineframe[i,10] == 0 else '~' if lineframe[i,14] != 0 and obs_flux> 0 else ''


