from jax import numpy as jnp

def get_grid():
    import pyscf
    from pyscf import gto, scf, dft
    from pyscf.data.nist import BOHR
    mol = pyscf.M(atom='''
                H  0. 0. 0.
                H  0. 0. 0.7408481486
            ''', basis='sto-3g')

    mf = scf.RKS(mol)
    # S = mf.get_ovlp(mol)  # overlap matrix

    mf.kernel()
    mf.grids.build(with_non0tab=True)
    Rmol = jnp.array(mf.grids.coords)
    W = jnp.array(mf.grids.weights)
    return Rmol, W

x,w = get_grid()

# print to file
with open("gridpoints.txt", "w") as f:
    for R in x:
        RR = str(R).strip("]").strip("[")+"\n"
        f.write(RR)
