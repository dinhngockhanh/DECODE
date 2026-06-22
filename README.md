#   Deconvolution of mutations into neutral tail and clusters, with corrections for sequencing and mutation calling biases

##  Installation

[DECODE (**De**ciphering **C**ancer **O**rigin from **D**NA **E**volution)](https://doi.org/10.64898/2026.06.15.732415) is an algorithm to decompose genomic variants into the neutral tail and mutation clusters based on their Variant Allele Frequencies (VAFs), with corrections for sample-specific DNA-sequencing coverage distribution and mutation calling biases in the Site Frequency Spectrum (SFS).

The DECODE library can be installed with

```R
devtools::install_github("dinhngockhanh/DECODE")
```

Detailed descriptions of how to run DECODE and analyze its output can be viewed in the [introductory vignette](https://dinhngockhanh.github.io/DECODE/vignettes/DECODE.html).

##  Methodology

DECODE infers clonality in a tumor sample based on the site frequency spectrum (SFS), the distribution of somatic mutations by their variant allele frequency (VAF).
The shape of the SFS encodes evolutionary history.
Clonal mutations form a cluster at high VAFs, and each subclone contributes a lower-VAF cluster.
Furthermore, regardless of whether the tumor evolves neutrally or undergoes selective sweeps, neutral mutations that accumulate as cells divide generate a characteristic power-law tail in the SFS.
By decomposing this spectrum and recovering the number and sizes of subclones present, DECODE estimates intra-tumor heterogeneity (ITH) and the ongoing tumor evolution.

DECODE is based on [our mathematical framework for the SFS](https://doi.org/10.1214/19-STS7561), which corrects for sample-specific sequencing coverage and mutation calling biases.
It implements [ABC-SMC-DRF](https://doi.org/10.1007/s11222-025-10748-x), our general likelihood-free inference method available as a stand-alone [R package](https://github.com/dinhngockhanh/abcsmcrf), which incorporates random forests into the framework of sequential Monte Carlo to accurately and efficiently infer the parameter posterior distribution.

<p align="center">
  <img src="Fig_schematics.jpg" alt="DECODE methodology" width="100%">
</p>
<p align="center">
  <small><em>
    A: correspondence between clonal evolution and the SFS. 
    B: schematic overview of DECODE's algorithm.
  </em></small>
</p>

Given a DNA-sequencing sample, DECODE first selects thresholds $(L,M)$ for three different data subsamples, termed **inference A**, **inference B** and **validation** (**step 1**).
The variant read count threshold is different in each subset, so the SFS from each filtered subsample assumes a different shape (**step 2**).
For a given cluster count $H$, DECODE applies ABC-SMC-DRF to infer the parameters $\theta_H=\left(\alpha,p_1,\dots,p_H,\omega_0,\omega_1,\dots,\omega_H\right)$, which characterize the exponent of the tail, each cluster's mean VAF, and the mutation count in each component.
ABC-SMC-DRF is modified to maintain cluster ordering, where cluster 1 is truncal and cluster $H$ corresponds to the rarest subclone.
DECODE thus finds the distribution of $\theta_H$ such that the predicted SFS best matches the empirical SFS from **inference A** and **inference B** subsamples (**step 3**).

To determine the parsimonious cluster configuration, DECODE tests the capacity of the inferred parameter distribution to predict the SFS from the **validation** subsample (**step 4**).
The accuracy of the prediction, balanced against the model's complexity as determined by cluster count $H$, is quantified with the Generalized Information Criterion (GIC) (**step 5**).
DECODE compares the GIC densities based on an increasing sequence of cluster counts.
It selects the more complex result as the better model if its GIC is lower, and continues adding more clusters (**step 6**).
Otherwise, it selects the model with the lower cluster count (**step 7**) as the most parsimonious decomposition (**step 8**).
DECODE performs this pipeline to fit the data with and without a neutral tail, then select the better model by GIC.

##  References

1.  Chen Y, Jaksik R, Terranova P, El Baghdadi S, Koval A, Kurpas MK, Tavaré S, Kimmel M, Dinh KN. [Accurate detection of tumor clonality and ongoing expansion mode from genomic data](https://doi.org/10.64898/2026.06.15.732415). bioRxiv (2026). 
2.  Dinh KN, Jaksik R, Kimmel M, Lambert A, Tavaré S. [Statistical inference for the evolutionary history of cancer genomes](https://doi.org/10.1214/19-STS7561). Statistical Science 35(1):129 (2020).
3.  Dinh KN, Liu C, Xiang Z, Liu Z, Tavaré S. [Approximate Bayesian computation sequential Monte Carlo via random forests](https://doi.org/10.1007/s11222-025-10748-x). Statistics and Computing 35, 219 (2025). 
