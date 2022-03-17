import numpy as np

import autotst.species

import zmatrix  # https://github.com/wutobias/r2z
from simtk import unit


def write_scan_file(fname, conformer, torsion_index, degree_delta=10.0):
    """Function to write a Gaussian rotor scan
    Takes an autoTST conformer and a rotor index
    """

    # header
    scan_job_lines = [
        "%mem=5GB",
        "%nprocshared=24",
        "#P m062x/cc-pVTZ",
        "Opt=(CalcFC,ModRedun)",
        "",
        "Gaussian input for a rotor scan",
        "",
        f"0 {conformer.rmg_molecule.multiplicity}",
    ]
    rdmol = conformer._rdkit_molecule
    cart_crds = np.array(rdmol.GetConformers()[0].GetPositions()) * unit.angstrom
    zm = zmatrix.ZMatrix(conformer._rdkit_molecule)

    zm_text = zm.build_pretty_zcrds(cart_crds)
    zm_lines = zm_text.splitlines()
    bonds = []
    angles = []
    dihedrals = []
    b = 1  # indices for bonds
    a = 1
    d = 1
    for line in zm_lines:
        tokens = line.split()
        # print(line)
        if len(tokens) == 1:
            pass
        elif len(tokens) == 3:
            bonds.append(f'B{b}        {tokens[2]}')
            tokens[2] = f'B{b}'
            b += 1
        elif len(tokens) == 5:
            bonds.append(f'B{b}        {tokens[2]}')
            tokens[2] = f'B{b}'
            b += 1
            angles.append(f'A{a}        {tokens[4]}')
            tokens[4] = f'A{a}'
            a += 1
        elif len(tokens) == 7:
            bonds.append(f'B{b}        {tokens[2]}')
            tokens[2] = f'B{b}'
            b += 1
            angles.append(f'A{a}        {tokens[4]}')
            tokens[4] = f'A{a}'
            a += 1
            dihedrals.append(f'D{d}        {tokens[6]}')
            tokens[6] = f'D{d}'
            d += 1
            tokens.append('0')
        else:
            raise NotImplementedError
        
        scan_job_lines.append(' '.join(tokens))
    
    scan_job_lines.append("")
    for bond in bonds:
        scan_job_lines.append(bond)
    for angle in angles:
        scan_job_lines.append(angle)
    for dihedral in dihedrals:
        scan_job_lines.append(dihedral)
    scan_job_lines.append("")

    # bond order and connectivity - not sure if this is needed

    # dihedral to scan
    indices = conformer.torsions[torsion_index].atom_indices
    N_increments = int(360.0 / degree_delta)
    scan_job_lines.append(f"D {indices[0] + 1} {indices[1] + 1} {indices[2] + 1} {indices[3] + 1} S {N_increments} {float(degree_delta)}")

    scan_job_lines.append("")
    with open(fname, 'w') as f:
        for line in scan_job_lines:
            f.write(line + '\n')
