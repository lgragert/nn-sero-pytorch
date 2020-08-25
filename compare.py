import pandas as pd

loci = ["A", "B", "C", "DPB1", "DQB1", "DRB1"]
for locus in loci:
    oldset = pd.read_csv("old_sets/test/" + locus + "_test.csv", low_memory=False)
    oldset = oldset.drop(columns=['serology'])
    newset = pd.read_csv("RSNNS_fixed/testing/" + locus + "_test.csv", low_memory=False)
    newset = newset.drop(columns=['serology'])

    oldset = oldset.set_index('allele')
    newset = newset.set_index('allele')

    newset = newset[newset.index.isin(oldset.index)]
    oldset = oldset[oldset.index.isin(newset.index)]

    col_drop = []
    for col in newset.columns:
        if newset[col].equals(oldset[col]):
            col_drop.append(col)

    oldlen = len(oldset)
    newlen = len(newset)

    merged = pd.concat([oldset,newset]).drop_duplicates(keep=False)
    merged = merged.drop(columns=col_drop)
    merged.to_csv(locus + "_diff.csv", index=True)
