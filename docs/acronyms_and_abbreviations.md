## Acronyms, Abbreviations and concept names used

| Acronym         | Topic     | Full                                | Explanation                                                                                  |
|-----------------|-----------|-------------------------------------|----------------------------------------------------------------------------------------------|
| Apollo          |           |                                     |                                                                                              |
| PPF             |           | Percent Point Function              |                                                                                              |
| SNP             | vcf       | Single Nucleotide Polymporphism     | single nucleotide variant                                                                    |
| MNP             | vcf       | Multi Nucleotide Polymporphism      | directly consequential SNPs (e.g. REF/ALT == AA/TT)                                          |
| INDEL           | vcf       | (IN)sertion or DEL(eletion)         | of a (short) DNA sequence                                                                    |
| "complex"       | vcf       | Complex variant                     | short DNA segment with various SNPs and indels considered as one larger variable micro-region|
| rule            | snakemake | single step in a pipeline           |                                                                                              |
| module          | snakemake | set of rules in separate *.smk file |                                                                                              |
| harmonize       | TODO      | pipeline harmonization project      | future refactoring/improvement/simplification once RIVM pipeline harmonisation kicks in      |
| blacklist       | pipeline  | blacklist of genomic regions (.bed) | regions identified to yield unreliable SNPs in a multi-clade concept                         |


## Variable names explained

| name            | topic | Explanation                                                                               |
|-----------------|-------|-------------------------------------------------------------------------------------------|
| checkpoint_path | rules | full path to a yaml file created by apollo_reference.matchreferencedata.matchref function |
|                 |       |                                                                                           |
|                 |       |                                                                                           |


