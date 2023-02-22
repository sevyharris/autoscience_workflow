import numpy as np

import autotst.species

import zmatrix  # https://github.com/wutobias/r2z
from simtk import unit


def write_scan_file(fname, conformer, torsion_index, degree_delta=20.0):
    """Function to write a Gaussian rotor scan
    Takes an autoTST conformer and a rotor index
    """

    # header
    scan_job_lines = [
        "%mem=5GB",
        "%nprocshared=48",
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

    # # bond order and connectivity - not sure if this is needed
    # rdkit_bonds = zm.rdmol.GetBonds()
    # bond_list = [f'{i + 1}' for i in range(0, len(rdkit_bonds))]
    # for bond in rdkit_bonds:
    #     bond_start = zm.a2z(bond.GetBeginAtomIdx())
    #     bond_end = zm.a2z(bond.GetEndAtomIdx())
    #     bond_order = bond.GetBondTypeAsDouble()
    #     bond_list[bond_start] += f' {bond_end + 1} {bond_order}'
    # for bond in bond_list:
    #     scan_job_lines.append(bond)
    # scan_job_lines.append("")

    # dihedral to scan
    indices = conformer.torsions[torsion_index].atom_indices

    # # don't convert to z-matrix index??
    # first = indices[0] + 1
    # second = indices[1] + 1
    # third = indices[2] + 1
    # fourth = indices[3] + 1

    # convert to z-matrix index
    first = zm.a2z(indices[0]) + 1
    second = zm.a2z(indices[1]) + 1
    third = zm.a2z(indices[2]) + 1
    fourth = zm.a2z(indices[3]) + 1

    N_increments = int(360.0 / degree_delta)
    scan_job_lines.append(f"D {first} {second} {third} {fourth} S {N_increments} {float(degree_delta)}")

    scan_job_lines.append("")
    with open(fname, 'w') as f:
        for line in scan_job_lines:
            f.write(line + '\n')
