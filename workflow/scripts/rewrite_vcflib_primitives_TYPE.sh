#!/usr/bin/env bash
# script that explicitly converts remaining TYPE=complex and TYPE=mnp in vcf INFO field;
# use-case is vcfallelicprimitives, that splits these but leaves the INFO tag unmodified
set -euo pipefail

awk 'BEGIN{FS=OFS="\t"}
/^#/ {print; next}
{
    len_ref = length($4)
    len_alt = length($5)
    new_type = ""

    if (len_ref == 1 && len_alt == 1) {
        new_type = "TYPE=snp"
    } else if (len_ref < len_alt && len_ref == 1) {
        new_type = "TYPE=ins"
    } else if (len_ref > len_alt && len_alt == 1) {
        new_type = "TYPE=del"
    }

    if (new_type != "") {
        gsub(/TYPE=[^;\t]+/, new_type, $8)
    }
    print
}'