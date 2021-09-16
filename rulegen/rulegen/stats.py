import rulegen.xyz
import pandas as pd
import seaborn as sns
import os
import matplotlib.pyplot as plt
from typing import List


def stringfiles2table(path: str):
    column_names = ["name", "node", "energy"]
    rows = []
    _, _, filenames = next(os.walk(path))

    for fn in filenames:
        if not fn.startswith("stringfile"):
            continue

        name = fn[-4:]
        reactionpath = rulegen.xyz.read_stringfile(os.path.join(path, fn))
        for node, m in enumerate(reactionpath.nodes):
            rows.append([name, node, m.energy])

    return pd.DataFrame(rows, columns=column_names)


def aggregate_stringfiles(path: str):
    tbl = stringfiles2table(path)
    print(tbl)
    tbl_barrier = tbl.groupby("name").energy.agg(["max", "last"])
    tbl_barrier = tbl_barrier.nsmallest(10, "max")
    print(tbl_barrier.index)
    sns.lineplot(data=tbl[tbl.name.isin(tbl_barrier.index)], x="node", y="energy", hue="name")
    plt.savefig("stat.pdf")
