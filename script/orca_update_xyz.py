#!/usr/bin/env python3
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='update the coordinates in a ORCA input file with xyz/PDB file.')

    # Add the arguments
    parser.add_argument('-inp',
                        required=True,
                        help='Template ORCA input file.')
    parser.add_argument('-pdb',
                        help='Coordinate input.')
    parser.add_argument('-xyz',
                        help='Coordinate input.')
    parser.add_argument('-o',
                        required=True,
                        help='Output ORCA input file.')

    # Parse the arguments
    # only one of pdb and xyz can be used
    args = parser.parse_args()
    if args.pdb and args.xyz:
        raise ValueError("Only one of pdb and xyz can be used.")
