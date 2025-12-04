import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt

def get_diffy_expressed_genes(filepath,n):
    """
    Transforms the cuffdiff dataframe specified at file_path to a pd:Dataframe with columns. Removes all rows with significance == no
    -  gene_id
    - log2(fold_change)
    - z_norm_log2FC: Normalized FC based on rows from the input dataframe whose log2(fold_change) != -inf. I 
    - p_value
    - q_value

    """
    D = pd.read_csv(filepath, delimiter="\t")
    #remove all rows whose log2FC is inf
    D = D.replace([float("inf"), float("-inf")], pd.NA).dropna(subset=["log2(fold_change)"])
    D["z_norm_log2FC"] = (D["log2(fold_change)"] - D["log2(fold_change)"].mean() ) / D["log2(fold_change)"].std()    
    D = D[D["significant"] == "yes"]

    D = D[ ["gene_id", "log2(fold_change)", "z_norm_log2FC", "p_value", "q_value"] ]
    D.sort_values("z_norm_log2FC",inplace=True,ascending=False)
    
    # print("upregulated:")
    # print(D.head(n))
    # print("downregulated:")
    # print(D.tail(n))
    # print("differentially regulated")
    
    D.sort_values(by="z_norm_log2FC", inplace=True, key = lambda x: abs(x))
    return D


def main():
    file = "diff/5s_vs_4s.gene_exp.diff"
    n = 20
    assert(len(sys.argv) in [1,2,3])
    if len(sys.argv) == 1:
        print(f"using default params:\n\tn ={n}\n\tinput file = {file}")
    else:
        file = sys.argv[1]
        n = int(sys.argv[2])

    D = get_diffy_expressed_genes(file,n)

    print(D.head(n).to_string())

if __name__ == "__main__":
    main()


