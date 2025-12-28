# Reproducibility supplement for "Distribution-Free Signs and Ranks via Optimal Transport under Multivariate Symmetry: Application to One-Sample Location Testing"

This repository contains code to reproduce the results in the supplementary file of the paper "Distribution-Free Signs and Ranks via Optimal Transport under Multivariate Symmetry: Application to One-Sample Location Testing".


### Organization of the reproducibility folder

The reproducibility folder is organized as follows.  

- **`OT-visualization.ipynb`**  
  Contains all code used to generate the visualization plots in Figures 1–4 of the supplementary file, illustrating optimal transport–based signs and ranks under multivariate symmetry.

- **`OTsymm_simulations.ipynb`**  
  Contains simulation code for Appendix D.1 and Appendix D.2. This notebook implements the OT-Wilcox and OT-sign testing procedures using both random and Halton-transformed reference vectors, produces results for the simulated examples (C1–C10, S1–S10, Sp1–Sp10), and includes a section to generate the raw outputs (`Gaussian.csv` and `Epanechnikov.csv`) used for Figure 5.

- **`SymOT.R`**  
  Reads `Gaussian.csv` and `Epanechnikov.csv` and generates Figure 5 in the supplementary file.

- **`OT_symm_real_data.ipynb`**  
  Reproduces all real-data analyses reported in Appendix D.3.

- **`Gaussian.csv` and `Epanechnikov.csv`**  
  Store precomputed simulation results for Figure 5 to avoid long computation times.

- **`metadata.csv`**  
  Contains the NASA prognostics data used in the real-data example.

- **`aids_clinical_trials_X.csv` and `aids_clinical_trials_y.csv`**  
  Contain the healthcare data used in the real-data example.



### Reproducibility process


1. To produce the visualization plots Figures 1-4 in the supplementary file, run the Jupyter notebook `OT-visualization.ipynb`.

    Note: We use the package "Python Optimal Transport" via `import ot`. If this package is not installed, run `%pip install pot` at the beginning of the notebook; otherwise, ignore that block.

2. To produce Figure 5 in the supplementary file, first run the section "Code to produce Figure 5" in the Jupyter notebook `OTsymm_simulations.ipynb`, which produces two files `Gaussian.csv` and `Epanechnikov.csv` that store the results. Then run `SymOT.R` which generates Figure 5 based on these two files. 

    Note: Obtaining `Gaussian.csv` and `Epanechnikov.csv` could take 1-2 days on a local machine (MacBook Pro 2.4 GHz Intel Core i5). Copies of these two files are saved in this repository.

3. To reproduce the results of the real data examples in Appendix D.3, run the Jupyter notebook `OT_symm_real_data.ipynb`.

    Note: Copies of the data sets used are included in this repository (NASA prognostics data `metadata.csv` and healthcare data `aids_clinical_trials_X.csv` and `aids_clinical_trials_y.csv`).

4. To implement our methods (OT-Wilcox and OT-sign) with random and Halton transformed reference vectors for simulated examples (C1-C10, S1-S10, Sp1-Sp10) in Appendix D.2, run the section "Code to implement OT-Wilcox, OT-sign for simulated examples in Appendix D.2" in the Jupyter notebook `OTsymm_simulations.ipynb`. 

    Note: Running all these simulations could take a long time. The results have already been printed in the notebook.


