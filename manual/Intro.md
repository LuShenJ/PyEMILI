
## Introduction

PyEMILI is a Python-based automatic line identification tool designed for analyzing weak lines in high-resolution spectra. It builds upon the foundation of EMILI, the first-generation line identification program developed by [B. Sharpee et al. (2003)](https://ui.adsabs.harvard.edu/abs/2003ApJS..149..157S/abstract). PyEMILI incorporates many updates and corrections to improve upon the original EMILI, addressing various details to enhance its accuracy and performance for modern spectral analysis.

PyEMILI is currently equipped to handle spectral data within the wavelength range of 1000–20000 Å, utilizing the comprehensive [Atomic Line List v3.00b4](https://www.pa.uky.edu/~peter/newpage/index.html), which contains approximately 900,000 atomic transitions in this wavelength range. Additionally, PyEMILI utilizes transition probabilities sourced from the [Kurucz Line List](http://kurucz.harvard.edu/linelists.html), ensuring accurate line identifications and reliable results for the analysis of emission and absorption features.

A key enhancement in the current version of PyEMILI is the introduction of the `Spec` module. This module significantly enhances the program’s ability to automatically or manually search for a large number of spectral lines in deep spectra. By using the `Spec` module, users can efficiently analyze complex spectral data and identify emission features that might otherwise be missed in manual analysis. More information can be found in [Automatic line Searching](./Line_finding.md).

### Workflow

![workflow](./pic/PyEMILI_workflow.png)

### Line Identification Method

PyEMILI is designed not to model spectra, but to assist in the identification of weak spectral lines by proposing multiple potential identifications (IDs) and ranking them based on a scoring system. This scoring evaluates how well each candidate line satisfies certain predefined criteria. Each spectral line identified by PyEMILI is assessed and scored based on three distinct components: **wavelength agreement (W)**, **predicted flux (F)**, and **multiplet check (M)**. The overall identification score for each putative ID is defined by the sum of these three components.

#### **1. Wavelength Agreement (W)**
The **wavelength agreement** score is based on how closely the observed line’s wavelength matches the wavelength of a potential identification after correcting for systemic velocity. The residual difference between the observed line's wavelength and the wavelength of the putative ID is expressed as a multiple of the wavelength uncertainty ($\sigma$, needed in the *Input Line List*) in the observed line's wavelength. This residual difference, denoted as $\Delta \lambda$ (in km/s), determines the **W** score according to the following conditions:

- **W = 0** if $\Delta\lambda \leq 1\sigma$
- **W = 1** if $1\sigma < \Delta\lambda \leq 2\sigma$
- **W = 2** if $2\sigma < \Delta\lambda \leq 3\sigma$
- **W = 3** if $3\sigma < \Delta\lambda \leq 4\sigma$
- **W = 4** if $4\sigma < \Delta\lambda \leq 5\sigma$

Lines with a closer match (lower $\Delta\lambda$) will receive a lower **W** score, indicating a better match.

#### **2. Predicted Template Flux (F)**
The **predicted template flux** in PyEMILI is scored based on the expected flux contribution from **radiative recombination**, **dielectronic recombination**, and **collisional excitation** processes. For radiative and dielectronic recombination, the flux is represented by:


$$I_\mathrm{R} \propto L_{\rm R} = N_{\rm e}N_{i+1}C(\alpha_{\rm RR} + \alpha_{\rm DR})h\nu,$$


where $N_{\rm e}$ is the electron density, $N_{i+1}$ is the ion density, $\alpha_{\rm RR}$ and $\alpha_{\rm DR}$ are the radiative and dielectronic recombination coefficients, respectively, and $h\nu$ is the photon energy.

For collisional excitation, the flux is represented as:

$$I_\mathrm{C} \propto L_{\rm C} = DN_{\rm e}N_{1}q_{12}h\nu\frac{1}{1 + N_{\rm e}q_{21}/A_{21}},$$

where $N_{1}$ is the number density populated in the lower energy state, $q_{12}$ and $q_{21}$ are the collisional excitation and de-excitation rates, and $A_{21}$ is the spontaneous transition probability.

The **predicted template flux** for each line, whether **$I_{\rm R}$** or **$I_{\rm C}$**, is calculated based on the transition type:
- **Permitted transitions** are assumed to form primarily through recombination.
- **Forbidden transitions** are assumed to arise from collisional excitation.
- For **intercombination transitions** (or semi-forbidden transitions), the predicted template flux is a sum of recombination and collisional excitation contributions, with the collisional excitation component diluted by a factor of 100. This factor is empirically derived from the comparisons of numerous PNe and H II regions samples.

The **predicted flux** score is based on the expected intensity of the spectral line. For each putative ID within a $5\sigma$ wavelength difference (surviving the wavelength agreement assessment), the predicted flux ($I$) is calculated as mentioned above. The highest predicted flux of an unidentified line is denoted as $I_{\text{max}}$. Each ID is scored based on how its predicted flux compares to $I_{\text{max}}$, following these criteria:

- **F = 0** if $I \geq 0.1I_{\text{max}}$
- **F = 1** if $I \geq 0.01I_{\text{max}}$
- **F = 2** if $I \geq 0.001I_{\text{max}}$
- **F = 3** if $I \geq 0.0001I_{\text{max}}$

Putative IDs with predicted fluxes less than $0.0001I_{\text{max}}$ are not considered further, as they are too weak to be reliable candidates. Higher fluxes lead to lower **F** scores, indicating a better match.

#### **3. Multiplet Check (M)**
The **multiplet check** evaluates whether the putative ID belongs to a multiplet and, if so, whether other members of the multiplet are also detected in the *Input Line List*. Each putative ID’s multiplet members are examined, and both the number of total multiplet members ($P$) and the number of detected multiplet members ($D$) are considered. The **M** score is determined according to the following conditions:

- **M = 0** if $P = 2$ & $D = 2$, or $D > 2$
- **M = 1** if $P = 1$ & $D = 1$, or $P > 2$ & $D = 2$
- **M = 2** if $P = 0$ & $D = 0$, or $P > 1$ & $D = 1$
- **M = 3** if $P = 1$ & $D = 0$, or $P = 2$ & $D = 0$
- **M = 4** if $P > 2$ & $D = 0$

The presence of multiple detected multiplet members improves the match, resulting in a lower **M** score.

#### **4. Identification Index (IDI)**
Finally, the **Identification Index (IDI)** is used to quantify the overall likelihood of each putative ID being the plausible identification. The **IDI** is calculated as the sum of the scores from the three components:  

$$ IDI = W + F + M. $$

Lower **IDI** values correspond to better matches, meaning the candidate line satisfies the criteria more closely. PyEMILI provides users with multiple potential IDs for each line, ranked by their **IDI** scores, allowing users to select the most plausible identification based on the context of analysis.

#### Summary of Identification Criteria
| W | Condition                   | F | Condition             | M   | Condition              |
|---|-----------------------------|---|-----------------------|-----|------------------------|
| 0 | $\Delta\lambda\leq 1 \sigma$          | 0 | $I\geq 0.1 I_{max}$    | 0  | $P = 2$ & $D = 2$, or $D > 2$     |
| 1 | $1\sigma <\Delta\lambda \leq 2\sigma$ | 1 | $I\geq 0.01 I_{max}$   | 1  | $P = 1$ & $D = 1$, or $P > 2$ & $D = 2$ |
| 2 | $2\sigma <\Delta\lambda \leq 3\sigma$ | 2 | $I\geq 0.001 I_{max}$  | 2  | $P = 0$ & $D = 0$, or $P > 1$ & $D = 1$ |
| 3 | $3\sigma <\Delta\lambda \leq 4\sigma$ | 3 | $I\geq 0.0001 I_{max}$ | 3  | $P = 1$ & $D = 0$, or $P = 2$ & $D = 0$ |
| 4 | $4\sigma <\Delta\lambda \leq 5\sigma$ |   |                        | 4  | $P > 2$ & $D = 0$              |

### Definition of energy bins and ionization & velocity structure models

All ions are separated into 5 bins based on the ionization energies to produce such ions. The bounds of each bin are listed below.

| bin 1   | bin 2       | bin 3       | bin 4      | bin 5  |
|---------|-------------|-------------|------------|--------|
|0-13.6 eV| 13.6-24.7 eV| 24.7-55 eV| 55-100 eV| >100 eV|

Each bin has the parameters of ionization and velocity structure models. The ionization structure model refers to the proportion of ions in a given ionization bin to the total elemental abundance of that ion. **The modified abundances of ions are derived through the total abundances of the elements multiplied by the values of the ionization structure model for the bins in which the ions reside**. Baldwin et. al. (2000) have shown that in the rest frame of a nebula, the magnitude of the velocity difference between the observed and laboratory wavelengths is correlated with the parent ions' ionization energies. Thus, the velocity structure model aims to correct this correlation for different ions with different ionization energies.

For example, the minimum energy to produce $\mathrm{O^{2+}}$ is 35.11 eV, which also means the ionization energy of the $\mathrm{O^{1+}}$. Thus, $\mathrm{O^{2+}}$ will be classified into the bin 3. And all the [O III] collisionally excited lines and O II recombination lines, whose parent ion is $\mathrm{O^{2+}}$, will use the corresponding values of ionization and velocity structure models for bin 3 to calculate their scores.

Initial ionization structure value for each bin:

| bin 1   | bin 2       | bin 3       | bin 4      | bin 5  |
|---------|-------------|-------------|------------|--------|
|0.01| 0.5| 0.4| 0.1| 0.001|

Initial velocity structure value for each bin (in km/s):

| bin 1   | bin 2       | bin 3       | bin 4      | bin 5  |
|---------|-------------|-------------|------------|--------|
|0| 0| 0| 0| 0|

### Output files

After running `pyemili.Lines.Line_list.identify()`, two files ending with **'.dat'** and **'.out'** will be generated in the directory. The '.out' file contains complete candidate IDs of each input observed line. As an example of the results, The following presents a demonstrative output of PyEMILI's identification of an emission line observed:

>`Number 98 	Observed line:	3918.93000 	 1.07E-03	SNR:226.4	FWHM:18.7`

The '.out' file is composed of individual blocks, with each block corresponding to the line identification of a single observed line. Each block begins with header information as shown above, which means this observed line is the 98th line of the *Input Line List*. The observed wavelength is 3918.93 angstroms. Flux is 1.07e-03 relative to the $\mathrm{H}_\beta=1$.

![out file](./pic/out.png)

>**Column (1)** Wavelength of a candidate ID corrected using the velocity-structure parameter.  The '`+`' notation in front of the wavelength means that the velocity difference in Column (9) is less than $1\sigma$.  
>**Column (2)** Laboratory wavelength of a candidate ID.  An asterisk '**$\ast$**' in front of the laboratory wavelength means it is a weighted-average wavelength if all the fine-structure transitions (including this candidate ID) within a multiplet have very close wavelengths, e.g., with differences smaller than the instrumental resolution.  No asterisk in this example, because it is not applicable to this case.  
>**Columns (3)** The emitting ion.  
>**Columns (4)** The lower and upper spectral terms.  
>**Columns (5)** The predicted template flux relative to $\mathrm{H}_\beta$. The value with a tilde '**$\sim$**' in front is calculated using the effective recombination coefficient, while that with an asterisk '**$\ast$**' in front indicates the transition has no transition probability in our atomic transition database and its predicted template flux is estimated using a default transition probability (see paper section 2.4.2).  
>**Columns (6)** The number of associated multiplet members that are expected to be observed, versus the number of transitions actually detected. For example, `5/1` means 5 multiplet lines should be possibly observed in the input line list, but only 1 multiplet line has been detected by the code.  
>**Columns (7)** The identification index (IDI, IDI=W+F+M) of the candidate ID.  
>**Column (8)** Ranking of the candidate ID assigned by PyEMILI.  The wedge symbol '**$\wedge$**' in front of a ranking means at least one detected multiplet member has been ranked as 'A' for the corresponding observed line.  In this case, the candidate ID is treated as one of the best identifications regardless of whether it has been ranked as 'A'.  
>**Columns (9)** Velocity difference (in km/s) between the wavelengths in Column (1) and Column (2).  
>**Columns (10) & (11)** The associated multiplet members (i.e. fine-structure components belonging to the same multiplet) detected from the *Input Line List*, and their velocity difference, respectively.  Up to three multiplet lines are displayed in this output.  The wedge symbol '**$\wedge$**' in front of the wavelength means this multiplet member is assigned an 'A' ranking in the corresponding observed line.  
>**Columns (12)** Electron configurations of the lower and upper states.  
>**Columns (13) & (14)** Statistical weights of the lower and upper states, respectively.  
>**NOTE**: In each block for one observed line, the last candidate line has the **IDI=99 in column (7)**. It indicates that this line does not participate in the ranking. This line is the predicted strongest line within a wavelength difference range of 5-10 $\sigma$. It's just shown here to remind you of another possibility if the input wavelength uncertainty is too low, or the possible line blending.

Another file with '.dat' is the output of the line list with all the A ranking IDs.

![dat file](./pic/dat.png)

>**First column** is the observed wavelength.  
>**Second column** is the line flux relative to $\mathrm{H}_\beta$.  
>**On the right side of the vertical bar** are candidates with A ranking, or candidates whose multiplet member is ranked as A in other observed lines.
