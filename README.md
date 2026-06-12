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
    Left: clonal evolution's connection to the SFS. 
    Right: schematic overview of DECODE's algorithm.
  </em></small>
</p>

DOWN