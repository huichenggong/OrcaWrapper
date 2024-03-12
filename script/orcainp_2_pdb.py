#!/usr/bin/env python3
import argparse



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert ORCA input to PDB format.')
    
    # Add the arguments
    parser.add_argument('-pdb',
                        required=True,
                        help='Template PDB file.')
    parser.add_argument('-inp',
                        required=True,
                        help='Coordinate will be read from this ORCA input file.')
    parser.add_argument('-o',
                        required=True,
                        help='Output pdb file.')
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Your script can now use args.pdb, args.inp, and args.o for further processing
    print(f"PDB file: {args.pdb}")
    print(f"INP file: {args.inp}")
    print(f"Output file: {args.o}")
