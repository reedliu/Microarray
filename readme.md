# Introduction to mRNA microarray

> **Creator:** Yunze Liu (Reed)  CAAS/AGIS/SDAU
>
> **Log:** created on 12/20/2018

### About terminology

- **target:** DNA hybridized to the array, mobile substrate
- **probe:** DNA spotted on the array (spot) ; one of these 25-mer oligonucleutides
- **probe set:** a collection of probes (e.g. 11) targeting the same transcript
- **print-tip-group**: collection of spots printed using the same print-tip (or pin), aka. grid.
- **reporter:** the sequence
- **feature:** a physical patch on the array with molecules intended to have the same reporter sequence (one reporter can be represented by multiple features)
- **accuracy:** how close are the estimated values to the truth
- **precision:**  how variable are the estimates
- **calibration/ normalization:**  adjust for systematic drifts associated with dye, array
- **background correction:**  adjust fot the non-linearity at the lower end of the dynamic range
- **PM:** Perfect Match; **MM:** MisMatch

### About microarray data

- Extracted from the images with artifacts (e.g. cross-talk) removed
- Final store in a matrix: row for probes, column for samples
- For each sample, each probe has one number from one-color arrays and two numbers for tow-color arrays 
- Expression levels for the same genes from different arrays can be **compared, after proper normalization**
- Only calculate relative expressions

### About history

- Evolved from Southern blotting, which is a procedure to detect and quantify a specific DNA sequence. Microarray can be thought as parallelized Southern blotting.
- First influential paper: Schena *et al*. (1995) *Science* [study the expression of 45 Arabidopsis genes]
- mRNA microarray got 50,000+ hits on pubmed for the past 20 years
- Main pros: lower costs, easier experimental procedure and more established analysis methods
- **Brief timeline**
  - Late 1980s: Lennon, Lehrach: cDNAs spotted on nylon membranes
  - 1990s: Affymetrix adapts oligonucleotide synthesis technology to build microchip with patent (So, the *Genechip* cannot be used by others)
  - 1990s: Brown lab in Stanford develops two-colour spotted array tech (open and free)
  - 1998: Yeast cell cycle expression profiling on spotted arrays (Spellmann) and Affymetrix (Cho)
  - 1999: Tumor type discrimination based on mRNA profiles (Golub)
  - 2000 - 2004: Affymetrix dominates the microarray market
  - Since 2003: Nimblegen, Illumina, Agilent
  - Since 2000: CGH, CNVs, SNPs, ChIP, tiling arrays
  - Since 2007: NGS (454, Solexa, ABI Solid, ...)

### About technology and design

![1.png](https://upload-images.jianshu.io/upload_images/9376801-d34fc2a99b777601.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

![3.png](https://upload-images.jianshu.io/upload_images/9376801-8a8553ee72af09dd.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
- Collection of DNA spot on a solid surface
- Each spot contains many copies of the same DNA sequence ("probes")
  [probes are designed to target specific genes]
- Part of genes sequence are complementary to a probe ("hybridize" or "stick to")
- The amount of mRNA for target gene is measured with the amount of hybridization

### Platforms

- Affymetrix
- Agilent
- Nimblegene
- Illumina
- ABI
- Spotted cDNA

### One-color vs. two-color arrays

> color means channel

- One-color arrays hybridize **one sample per array**
  ex. **Affymetrix, Illumina**
  Easier but need twice arrays
- Two-color arrays hybridize **two samples on the same array**
  ex. **Agilent, Nimblegen**
  Each spot produce two numbers

---

### Most famous - Affymetrix
![2.png](https://upload-images.jianshu.io/upload_images/9376801-22b32206524ffae1.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

#### U133 chip

- HG U133A genechip represents more than 22,000 full length genes and EST clusters
- Around 20 probes (11 probe pairs)  in a probe set to target the same transcript
- Not necessarily evenly spaced: sequence property matters
- Probes are randomly located on the chip to avergae out the effects of the array surface
- High signal intensity and lower noise

![6.png](https://upload-images.jianshu.io/upload_images/9376801-3d5397af68437ba1.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

#### Probe set naming

- ..._at: (anti-sense target) detects antisense strand of given gene, these are unique probes

- _a_at: (gene family probeset) recognize multiple transcripts of same gene

  _s_at: (identical probeset) recognize multiple transcripts from different genes

  _x_at: (mixed probe set) cross-hybrization with other sequences used for design

  > a, s, x derived from gene cluster + gene family info 



---

### Analysis Tasks

- IMPORT : reformating and setup/curation of the metadata
- NORMALIZATION
- QUALITY ASSESSMENT & CONTROL
- DIFFERENTIAL EXPRESSION
- Annotation
- Gene set enrichment
- Clustering and classification
- Integration of other datasets

### Analysis challenges

- **Data normalization:** remove systematic technical artifacts
  - Within array: variations of probe intensities
    - cross-hybridization: probe capture the "wrong " target
    - probe sequences: some probes bind too tight like "sticker"
    - chip factors: spot sizes, smoothness of array surface
  - Between array: variations of sample processing, image reader etc.
- **Calculation of gene expressions**: standard of summarizing multiple probes belonging to the same gene into one number
- **Differential expression detection**: find different genes between different conditions, e.g. case vs. control

---

### Data normalization

#### Why?

Artifacts are introduced:

- Sample preparation: PCR effects
- Array itself: array surface effects, printing-tip effects 
- Hybridization: non-specific binding, GC effects
- Scanning: scanner effects

> So, the normalization is mean to ensure differences in intensities due to differential expressions, not artifacts

#### Types

- Within-array: remove array-specific artifacts individually and obtain true signals 
- Between-array: put different arrays at the same baseline to make sure that numbers are comparable

#### First: Within array normalization (two-color)

Most common problem: intensity dependent effect

Most popular: **loess normalization** 

##### What is MA plot?

Widely used diagnostic plot fot microarray and sequencing data.

**M** measures relative expression; **A** measures total expression

##### What is loess normalization?

ASSUMPTIONS : 1. most genes are not DE (M=0); 2. M and A are independent (MA plot should be flat and centered at 0)

Loess (lowess): locally weighted scatterplot smoothing [to fit a smooth curve between two variables]

![4.png](https://upload-images.jianshu.io/upload_images/9376801-801df7f790b70522.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

#### Second: Within array normalization (one-color)

Background error models: **RMA( Robust Multi-array Average**) 
Model fitting: Median Polish

> Also, the MA plot and background error models (e.g. RMA) are popular in other microarrays, ChIP-seq, RNA-seq

#### Third: Between array

Artifacts reasons could be:

-  Total amount of mRNA used
- Properties of the agents used 
- Array properties
- Settings of laser scanners

##### Lineat scaling method

Affymetrix software **MAS**: use a number of "housekeeping" genes and assume their expressions are identical, then rescale all data [ based on one-step Tukey Biweight (TBW) ]

##### Non-linear smoothing based

used in dChip: based on PM-MM

PROCEDURE: find a set of genes invariant across arrays => find a "baseline" array => ohter arrays fit a smooth curve on expressions of invariant genes => normalize based on the fitted curve (Q-Q plot)

##### Quantile normalization

Force the distribution of all data from all arrays to be the same, but keep the ranks of the genes

>  Simple, Fast, Easy and can correct for quite nasty non-linearities (saturation, background); no matter how the data is good/bad

PROCEDURE: find SMALLEST value for each sample => AVERAGE them => replace each value by the average => find the NEXT SMALLEST, then AVERAGE => replace 

RESULT: 1. The values taken in each column are exactly the same; 2. The ranks of genes in each column are the same as before normalization

CAUTION: QN is too strong and often remove the true signals

![5.png](https://upload-images.jianshu.io/upload_images/9376801-0aafc8a59cb25547.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

---

### Microarray with bioconductor

> 200+ packages for microarray

#### Platform-specific data import

http://www.bioconductor.org/docs/workflows/oligoarrays/

- Affymetrix 3' IVT (e.g. Human U133 Plus 2.0, Mouse 430 2.0) => **affy**
  [ one of the earliest packages ]
- Affymetrix Exon (e.g. Human Exon 1.0 ST) => **oligo, exonmap, xps**
- Affymetrix SNP arrays => **oligo**
  [ designed to replace affy package ]
- Nimblegen tiling arrays (e.g. for ChIP-chip)=> **Ringo**
- Affymetrix tiling arrays (e.g. for ChIP-chip) => **Starr**
- Illumina bead arrays => **beadarray, lumi**



