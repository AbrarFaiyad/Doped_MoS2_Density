#!/usr/bin/env python3
"""
Quick test of the MoS2 structure building function using ASE's mx2
"""

import sys
import os

# Add the src directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

try:
    from mos2_npt_simulation import build_mos2_structure
    
    print("Testing MoS2 structure builder using ASE's mx2 function...")
    
    # Create a small test structure
    atoms = build_mos2_structure(supercell_size=(2, 2), layers=2, vacuum=10)
    
    print(f"Success! Created structure with {len(atoms)} atoms")
    print(f"Chemical symbols: {set(atoms.get_chemical_symbols())}")
    print(f"Cell: {atoms.get_cell()}")
    print(f"Volume: {atoms.get_volume():.2f} Å³")
    
    # Calculate expected number of atoms
    expected_atoms = 2 * 2 * 2 * 3  # supercell_x * supercell_y * layers * (1 Mo + 2 S per unit cell)
    print(f"Expected atoms: {expected_atoms}, Actual: {len(atoms)}")
    
    # Save test structure
    from ase.io import write
    write("test_structure.xyz", atoms)
    print("Test structure saved as test_structure.xyz")
    
    # Print first few atomic positions
    print("\nFirst few atomic positions:")
    for i in range(min(6, len(atoms))):
        pos = atoms.positions[i]
        symbol = atoms.symbols[i]
        print(f"  {symbol}: {pos[0]:8.3f} {pos[1]:8.3f} {pos[2]:8.3f}")
    
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
