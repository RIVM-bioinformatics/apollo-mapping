""" python code to filter variants, which is to complex for the shell: directive """

# Exact code below was first a run: directive, which does not support conda environment.
# Here, arguments known in the *.smk are extracted from global snakemake
from snakemake.shell import shell
params = snakemake.params
input = snakemake.input
output = snakemake.output
log = snakemake.log

class AttrDict(dict):
    """ mixin between dictionary and object, with both options for key lookup """

    # TODO: refactor / move elsewhere!
    # https://stackoverflow.com/questions/4984647/accessing-dict-keys-like-an-attribute
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

params.ths = AttrDict(params.ths)

# -----------------------------------------------------------------
# start of originally written run: directive
# -----------------------------------------------------------------

# gather extra header lines
extra_header_lines = []

# Always filter indels and complex.
steps = [
    "bcftools filter -e 'TYPE=\"indel\"' -s Indel -m+ {input.vcf}",
    "bcftools filter -e 'INFO/TYPE=\"complex\"' -s Complex -m+"
]

# In 'species" mode AND multiclade, filter on the blacklist
# TODO: !!important!! at this point the pipeline will crash in some conditions
#       e.g. non C.auris ;-) e.g. not-multiclade.
#       Need to get fixed ASAP ... but can't be done without extra input
#       Probably it's needed to add to input: ref=MultiReferenceProvider.get_ref_path
#       Based on this, we can check if it's linked to a blacklist
if params.trigger_multiclade_masking_workflow == "True" and wildcards.ref_type == "species":
    header_line = '##FILTER=<ID=Blacklisted,Description="Variant overlaps with softclipped and/or segmental variable regions">'
    steps.append(
        """bcftools annotate -a <(cat {input.bed} | awk '{{ OFS=\"\\t\"; print $1,$2+1,$3,\"Blacklisted\" }}') """
        f"-c CHROM,FROM,TO,FILTER "
        f'--header-lines <(echo \'{header_line}\')'
    )

# Conditionally (default=yes) filter on quality
if not params.omit_qual_thresholding:
    header_qual = '##FILTER=<ID=QualFit,Description="Variant quality (QUAL) is outside the fitted model thresholds">'
    extra_header_lines.append(header_qual)
    steps.append("bcftools filter -e 'QUAL < {params.ths.min_QUAL} || QUAL > {params.ths.max_QUAL}' -s QualFit -m+ "
                 )

# Conditionally (default=yes) filter on coverage
if not params.omit_cov_thresholding:
    header_cov = '##FILTER=<ID=CovFit,Description="Variant coverage (INFO/DP) is outside the fitted model thresholds">'
    extra_header_lines.append(header_cov)
    steps.append("bcftools filter -e 'INFO/DP < {params.ths.min_DP} || INFO/DP > {params.ths.max_DP}' -s CovFit -m+ "
                 )

if extra_header_lines:
    # add extra header lines; better done using (blanco) annotate, not bcftools reheader
    header_lines = "\\n".join(extra_header_lines)
    steps.append(
        f"bcftools annotate --header-lines <(echo -e '{header_lines}')"
    )

## Add write-to-output-file to the final command, merge CLI commands and add stderr logfile writing
steps[-1] += " -o {output.vcf}"
full_command = " | ".join(steps)
full_command += " > {log} 2>&1"

# Execute the full command!
shell(full_command)