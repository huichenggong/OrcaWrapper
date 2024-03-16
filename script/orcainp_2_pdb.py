#!/usr/bin/env python3
import argparse
import MDAnalysis as mda
from orcawrapper import orca_input

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert ORCA input to PDB format. If pdb file is given, residue '
                                                 'information will be read from pdb file.')

    # Add the arguments
    parser.add_argument('-pdb',
                        help='Template PDB file.')
    parser.add_argument('-inp',
                        required=True,
                        help='ORCA input file. Coordinate will be read from this file.')
    parser.add_argument('-o',
                        required=True,
                        help='Output pdb file.')

    # Parse the arguments
    args = parser.parse_args()

    # Your script can now use args.pdb, args.inp, and args.o for further processing
    print(f"pdb file: {args.pdb}")

    if args.pdb:
        print(f"pdb file: {args.pdb}")
        print(f"Orca input file: {args.inp}")
        print(f"Output file: {args.o}")
        u = mda.Universe(args.pdb)
        orca_inp = orca_input(args.inp)
        u.atoms.positions = orca_inp.xyz
        u.atoms.write(args.o)
    else:
        print(f"Orca input file: {args.inp}")
        print(f"Output file: {args.o}")
        print("Pdb file is not given. Writing coordinate with element name has not been implemented yet.")
