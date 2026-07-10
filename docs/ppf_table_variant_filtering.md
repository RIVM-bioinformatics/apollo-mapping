## PPF table thresholds for variant filtering

### The Percent Point Function (PPF)

The **Percent Point Function (PPF)** is also known as the **inverse cumulative distribution function (inverse CDF)** or the **quantile function**. This function calculates the threshold value below which a given probability falls.

For a specified probability $p$, the PPF answers the following question:
> *"What value $x$ satisfies the equation: $P(X \le x) = p$?"*

### PPF thresholds for variant filtering.

In this pipeline PE Illumina fastq data is mapped to its corresponding reference species (fasta/genome). Later, variants are called using freebayes, yielding (raw) variant calling files (vcf). Both the observed genomic coverage (from BigWig) and the observed variant scores, are converted into a modelled (normal) distribution describing the fit.

Based on these modelled distributions, (raw) called variants are filtered out as being unlikely to be biologically true. The modelled fits are transformed into Percent Point Function (PPF) tables, from which the user can provide in the CLI tool of this pipeline (apollo_mapping.py) at which lower and/or upper boundary of these distributions, variants should be discarded as being unlikely true.

### Choosing PPF thresholds using CLI flags

In a most simplified example, one could decide to use ```--min-ppf-QUAL-species 0.01``` and ```--max-ppf-QUAL-species 0.99```, meaning that the lower 1% and upper 1% of the distrubution's tails are discarded.

Realize this can be adjusted by:
- lower (--min-*) and upper (--max-*) direction of the scale independently
- at the --*-species and --*-strain mapping level independently
- so make sure you realize the difference, and the meaning of the "clade" or "strain" concept in this pipeline!
- thresholds are asymmetrical! So 0.01 and 0.99 combined remove 1+1 = 2% of the total width of the distribution
- one can choose value ```NA```, meaning that at this direction the thresholding is omitted
- one can use e.g. ```--omit-ppf-QUAL-species```, which is identical to ```--min-ppf-QUAL-species NA --max-ppf-QUAL-species NA ```

### Traceability of variant filtering.

In the output files of the variant calling and filtering phase of this pipeline, 5 files per sample are generated:

1. variants/raw/species-SRR14906880.vcf
   - raw, unfiltered variants
   - rule: freebayes in call_variants.smk
2. variants/norm/species-SRR14906880.vcf
   - normalized, per SNP/indel splitted representation (splits TYPE=complex and TYPE=mnp)
   - rule: vcflib_breakmulti_wave in call_variants.smk
   - realize number of variants in norm > raw !
3. variants/filtered/SRR14906880_on_species-labeled.vcf
   - each variant has FILTER field stated with PASS or semicolon-separated (multiple) reasons-of-FAIL:
     - Blacklisted: variant disqualified by blacklist overlap (only in case of multiclade concept)
     - Complex: variant disqualified because it's a TYPE=complex one
     - Indel: variant disqualified because it's an indel
     - CovFit: variant disqualified because it's out of bounds of the PPF-thresholds from the fitted WGS coverage
     - QualFit: variant disqualified because it's out of bounds of the PPF-thresholds from the fitted QUAL score
   - realize most disqualified SNPs are disqualified by multiple reasons.
   - realize number of variants in labeled == norm !
4. variants/filtered/SRR14906880_on_species-simplified.vcf
   - all non-PASS value in the FILTER field rewritten as explicit FAIL
   - realize number of variants in labeled == simplified == norm !
5. variants/filtered/SRR14906880_on_species-passed.vcf
   - only all variants (in current pipeline SNPs only) with FILTER value PASS, thus accepted SNPs
   - realize number of variants in passed < simplified !

Realize that in 'multiclade' fashion all files will occur twice (so 2x5 files, the extra ones with label "strain")
