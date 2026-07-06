#!/usr/bin/env python3
""" output a yaml containing basic statistics to stdout from one-value-per-line numerical input"""

# Python Imports
import pandas as pd
import sys

def main() -> None:
    """ """
    # read stdin
    df = pd.read_csv(sys.stdin,sep="\t",header=None)
    C = 'vals'
    df.columns = [C]

    # feature: "gracefully" solve number of decimals
    if df.dtypes[C] == int:
        num_decimals = 3
        def summed(v):
            return v
    else:
        df['decimals']=df[C].astype('str').str.split('.', expand=True)[1].apply(lambda x: len(x))
        num_decimals = int(df['decimals'].median()) + 1
        del(df['decimals'])
        def summed(v):
            return round(v,num_decimals-1)

    # generate stats
    stats = {
        "model":    "norm",
        "n":        df.shape[0],
        "sum":      summed(df[C].sum()),
        "min":      df[C].min(),
        "max":      df[C].max(),
        "median":   df[C].median(),
        "mean":     round(df[C].mean(),num_decimals),
        "stdv":     round(df[C].std(),num_decimals),
    }

    indentation=""
    if "--indented" in sys.argv:
        print("fit:")
        indentation="  "
    # print stats in yaml format
    for k,v in stats.items():
        print("%s%s: %s" % (indentation,k,v))
    sys.exit(0)

if __name__ == "__main__":
    main()