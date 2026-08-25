## required files for multiclade central reference blacklisting.

Apollo-reference supports multi-species via a (preloaded) database concept.
A customization to this, is that different (highly distinct) clades of the same species can be defined in a "cladegroup".
Current single example is Candia auris (Canidozyma auris), which has at this moment in time 6 defined claded (I-II-III-IV-V-VI).
In order to allow across-clades SNP comparisons, regions where SNP calling is not reliable (in **any**, not per se all of the clades!), are "blacklisted".
This is explained more detailedly in xxxx and yyyyy.
- data/cauris-GCA_002759435.3-vs-fastq-blacklist.bed contains all (bedtools merge) regions with access softclipping.
- data/cauris-GCA_002759435.3-vs-WGA-blacklist.bed contains all inter-assembly comparisons, emphasising where a.o. segmental deletions are located
- data/cauris-GCA_002759435.3-blacklist.bed contains the (bedtools merge) union of these two above stated files
 
### what to do with these files?
Since the workflows that generate this blacklist aren't fully automated (yet), for now the required files are deposited here.
- Please copy them to your downloaded "reference database folder" (see apollo-reference)
- Place all these *.bed files in the literal subfolder "multiclade"
- Don't rename any of these files.
- apollo-mapping.py, when called in multiclade mode, will check (and if not fail directly) if these files are present
- At this moment, data/cauris-GCA_002759435.3-vs-fastq-blacklist.bed still lacks clade III data (since we've failed to identify a publicly available fastq sample).
- At this moment, data/cauris-GCA_002759435.3-vs-WGA-blacklist.bed is still empty (NA)

### code snippet to run when first-ever piperun fails on this scenario

- assuming your working directory is (freshly cloned) apollo_mapping
- assuming the path to your reference species database is /mnt/db/apollo/reference

```bash
mkdir -p /mnt/db/apollo/reference/multiclade/
cp data/*.bed  /mnt/db/apollo/reference/multiclade/
```
