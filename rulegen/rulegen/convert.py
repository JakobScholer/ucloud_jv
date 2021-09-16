import rulegen.xyz
import os


def reactionpath2rxn(reactionpath):
    obr = reactionpath.nodes[0].obmol.write("mdl")
    obp = reactionpath.nodes[-1].obmol.write("mdl")
    obrxnstr = f"""
$RXN

OpenBabel

1  1
$MOL
{obr}
$MOL
{obp}

    """
    return obrxnstr


def stringfiles2rxn(path, outdir=None):
    if outdir is not None and not os.path.exists(outdir):
        raise Exception("Output directory does not exists.")

    rps = rulegen.xyz.read_stringfiles(path)
    for rp in rps:
        res = reactionpath2rxn(rp)
        if outdir is None:
            print(res)
        else:
            with open(f"{outdir}/{rp.name}.rxn", "w") as f:
                f.write(res)


