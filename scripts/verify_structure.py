#!/usr/bin/env python3
"""
Verify the initial MoS2 structure before running NPT simulation
"""

import sys
import os
import numpy as np

# Add the src directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

try:
    from mos2_npt_simulation import build_mos2_structure
    from ase.io import write
    from ase.neighborlist import NeighborList, natural_cutoffs
    
    print("Building and verifying MoS2 structure...")
    
    # Create test structure 
    atoms = build_mos2_structure(supercell_size=(10, 10), layers=8, vacuum=0)
    
    print(f"\nStructure verification:")
    print(f"  - Total atoms: {len(atoms)}")
    print(f"  - Chemical symbols: {set(atoms.get_chemical_symbols())}")
    
    # Check atomic positions
    positions = atoms.get_positions()
    print(f"  - Position ranges:")
    print(f"    X: {positions[:,0].min():.2f} to {positions[:,0].max():.2f} Å")
    print(f"    Y: {positions[:,1].min():.2f} to {positions[:,1].max():.2f} Å") 
    print(f"    Z: {positions[:,2].min():.2f} to {positions[:,2].max():.2f} Å")
    
    # Check cell parameters
    cell = atoms.get_cell()
    print(f"  - Cell vectors:")
    print(f"    a: [{cell[0,0]:.2f}, {cell[0,1]:.2f}, {cell[0,2]:.2f}] Å")
    print(f"    b: [{cell[1,0]:.2f}, {cell[1,1]:.2f}, {cell[1,2]:.2f}] Å")
    print(f"    c: [{cell[2,0]:.2f}, {cell[2,1]:.2f}, {cell[2,2]:.2f}] Å")
    print(f"  - Volume: {atoms.get_volume():.2f} Å³")
    
    # Check nearest neighbor distances
    try:
        cutoffs = natural_cutoffs(atoms)
        nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
        nl.update(atoms)
        
        min_distance = float('inf')
        for i in range(len(atoms)):
            indices, offsets = nl.get_neighbors(i)
            if len(indices) > 0:
                distances = atoms.get_distances(i, indices, mic=True)
                min_distance = min(min_distance, distances.min())
        
        print(f"  - Minimum interatomic distance: {min_distance:.3f} Å")
        
        if min_distance < 1.0:
            print("  ⚠️  WARNING: Very short bond detected!")
        elif min_distance > 4.0:
            print("  ⚠️  WARNING: No close neighbors found!")
        else:
            print("  ✓ Interatomic distances look reasonable")
            
    except Exception as e:
        print(f"  Could not compute neighbor distances: {e}")
    
    # Save structure for inspection
    write("verified_structure.xyz", atoms)
    print(f"\n✓ Structure saved as 'verified_structure.xyz'")
    
    # Check stoichiometry
    symbols = atoms.get_chemical_symbols()
    mo_count = symbols.count('Mo')
    s_count = symbols.count('S')
    ratio = s_count / mo_count if mo_count > 0 else 0
    
    print(f"\nStoichiometry check:")
    print(f"  - Mo atoms: {mo_count}")
    print(f"  - S atoms: {s_count}")
    print(f"  - S:Mo ratio: {ratio:.2f} (should be 2.00)")
    
    if abs(ratio - 2.0) < 0.01:
        print("  ✓ Correct MoS2 stoichiometry")
    else:
        print("  ⚠️  WARNING: Incorrect stoichiometry!")
    
    print("\n" + "="*50)
    print("STRUCTURE VERIFICATION COMPLETE")
    print("="*50)
    
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
