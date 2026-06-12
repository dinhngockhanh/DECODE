#   Deconvolution of mutations into neutral tail and clusters, with corrections for sequencing and mutation calling biases

##  Installation

DECODE (Deciphering Cancer Origin from DNA Evolution) is an algorithm to decompose genomic variants into the neutral tail and mutation clusters based on their Variant Allele Frequencies (VAFs), with corrections for sample-specific DNA-sequencing coverage distribution and mutation calling biases in the Site Frequency Spectrum (SFS).

The DECODE library can be installed with

```R
devtools::install_github("dinhngockhanh/DECODE")
```

Detailed descriptions of how to run DECODE and analyze its output can be viewed in the INTRODUCTORY VIGNETTE.

##  Methodology

UP

<p align="center">
  <img src="Fig_schematics.jpg" alt="DECODE methodology" width="100%">
</p>
<p align="center">
  <small><em>
    Left: DECODE's model for clonal evolution and the corresponding Site Frequency Spectrum (SFS). 
    Clonal mutations are present in a tumor's MRCA and shared across all cancer cells, forming the truncal cluster in the SFS.
    Each subclone results from an ongoing selection sweep, where the mutations in the subclone's MRCA are present in its cells and constitute a subclonal cluster at lower Variant Allele Frequencies (VAFs).
    Mitoses of cells within all subclones induce additional mutations, resulting in a neutral tail in the lowest-VAF region of the SFS.
    Right: Schematic overview of DECODE's algorithm.
    For each genomic sample, DECODE extracts the SFS from three filtering strategies, based on variant and total read counts (\textit{\textbf{step 1-2}}).
    For each cluster count $H$, it estimates model parameters with ABC-SMC-DRF from two inference SFS subsamples (\textit{\textbf{step 3}}), then finds the Generalized Information Criterion (GIC) density based on the predicted SFS for the validation subsample (\textit{\textbf{step 4-5}}).
    The GIC distributions from successive $H$'s are compared to find the parsimonious fit (\textit{\textbf{step 6-8}}).
  </em></small>
</p>

DOWN