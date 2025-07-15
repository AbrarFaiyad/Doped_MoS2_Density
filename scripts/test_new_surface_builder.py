#!/usr/bin/env python3
"""
Test script to verify the new MoS2 surface building using SurfaceBuilder.

This script will:
1. Test the new SurfaceBuilder-based MoS2 construction
2. Compare with the old mx2-based approach
3. Validate structure properties and stoichiometry
4. Export structures for visualization
"""

import sys
import os
sys.path.append('../src')

from mos2_npt_simulation import build_mos2_structure
from ase.io import write
from ase.build import mx2
import numpy as np

def test_new_surface_builder():
    """Test the new SurfaceBuilder-based MoS2 construction."""
    print("="*60)
    print("TESTING NEW SURFACEBUILDER-BASED MoS2 CONSTRUCTION")
    print("="*60)
    
    try:
        # Test parameters
        supercell_size = (3, 3)  # Smaller for testing
        layers = 2
        
        print(f"\nBuilding MoS2 structure with SurfaceBuilder:")
        print(f"  - Supercell: {supercell_size[0]}x{supercell_size[1]}")
        print(f"  - Layers: {layers}")
        
        # Build structure using new method
        atoms_new = build_mos2_structure(
            supercell_size=supercell_size,
            layers=layers,
            vacuum=10.0
        )
        
        # Analyze structure
        symbols = atoms_new.get_chemical_symbols()
        cell = atoms_new.get_cell()
        positions = atoms_new.get_positions()
        
        print(f"\n✓ Structure Analysis:")
        print(f"  - Total atoms: {len(atoms_new)}")
        print(f"  - Mo atoms: {symbols.count('Mo')}")
        print(f"  - S atoms: {symbols.count('S')}")
        print(f"  - S/Mo ratio: {symbols.count('S')/symbols.count('Mo'):.3f}")
        print(f"  - Cell vectors:")
        for i, vec in enumerate(cell):
            print(f"    {chr(97+i)}: [{vec[0]:8.3f}, {vec[1]:8.3f}, {vec[2]:8.3f}] Å")
        print(f"  - Volume: {atoms_new.get_volume():.2f} Å³")
        
        # Calculate density
        masses = atoms_new.get_masses()
        total_mass = np.sum(masses) * 1.66054e-24  # amu to grams
        volume_cm3 = atoms_new.get_volume() * 1e-24  # Å³ to cm³
        density = total_mass / volume_cm3
        print(f"  - Density: {density:.3f} g/cm³")
        
        # Check layer structure
        z_coords = positions[:, 2]
        print(f"  - Z-coordinate range: {z_coords.min():.3f} to {z_coords.max():.3f} Å")
        print(f"  - Structure height: {z_coords.max() - z_coords.min():.3f} Å")
        
        # Identify layers by grouping atoms with similar z-coordinates
        sorted_z = np.sort(np.unique(np.round(z_coords, 1)))
        print(f"  - Detected {len(sorted_z)} distinct z-levels")
        
        # Analyze Mo-S distances
        mo_indices = [i for i, sym in enumerate(symbols) if sym == 'Mo']
        s_indices = [i for i, sym in enumerate(symbols) if sym == 'S']
        
        if mo_indices and s_indices:
            min_mo_s_dist = float('inf')
            for mo_idx in mo_indices[:5]:  # Check first 5 Mo atoms
                mo_pos = positions[mo_idx]
                for s_idx in s_indices:
                    s_pos = positions[s_idx]
                    dist = np.linalg.norm(mo_pos - s_pos)
                    min_mo_s_dist = min(min_mo_s_dist, dist)
            print(f"  - Minimum Mo-S distance: {min_mo_s_dist:.3f} Å")
        
        # Save structure
        write("test_mos2_surfacebuilder.xyz", atoms_new)
        print(f"\n✓ Structure saved as 'test_mos2_surfacebuilder.xyz'")
        
        return atoms_new
        
    except Exception as e:
        print(f"\n❌ Error in SurfaceBuilder test: {e}")
        raise

def compare_with_mx2():
    """Compare new method with old mx2 approach."""
    print("\n" + "="*60)
    print("COMPARING WITH OLD MX2 METHOD")
    print("="*60)
    
    try:
        # Test parameters
        supercell_size = (3, 3)
        layers = 2
        
        print(f"\nBuilding MoS2 structure with mx2:")
        
        # Build using old mx2 method
        atoms_old = mx2(
            formula='MoS2',
            kind='2H',
            a=3.18,
            thickness=3.19,
            size=(supercell_size[0], supercell_size[1], layers),
            vacuum=10.0
        )
        atoms_old.set_pbc([True, True, True])
        
        # Analyze old structure
        symbols_old = atoms_old.get_chemical_symbols()
        
        print(f"✓ MX2 Structure Analysis:")
        print(f"  - Total atoms: {len(atoms_old)}")
        print(f"  - Mo atoms: {symbols_old.count('Mo')}")
        print(f"  - S atoms: {symbols_old.count('S')}")
        print(f"  - S/Mo ratio: {symbols_old.count('S')/symbols_old.count('Mo'):.3f}")
        print(f"  - Volume: {atoms_old.get_volume():.2f} Å³")
        
        # Calculate density for comparison
        masses_old = atoms_old.get_masses()
        total_mass_old = np.sum(masses_old) * 1.66054e-24
        volume_cm3_old = atoms_old.get_volume() * 1e-24
        density_old = total_mass_old / volume_cm3_old
        print(f"  - Density: {density_old:.3f} g/cm³")
        
        # Save old structure
        write("test_mos2_mx2.xyz", atoms_old)
        print(f"✓ MX2 structure saved as 'test_mos2_mx2.xyz'")
        
        return atoms_old
        
    except Exception as e:
        print(f"\n❌ Error in MX2 comparison: {e}")
        print("Note: mx2 function may not be available in current ASE version")
        return None

def main():
    """Main test function."""
    print("Testing MoS2 Structure Building Methods")
    print("="*60)
    
    # Test new SurfaceBuilder method
    atoms_new = test_new_surface_builder()
    
    # Compare with old mx2 method
    atoms_old = compare_with_mx2()
    
    # Summary comparison
    if atoms_old is not None:
        print("\n" + "="*60)
        print("COMPARISON SUMMARY")
        print("="*60)
        
        print(f"SurfaceBuilder method:")
        print(f"  - Atoms: {len(atoms_new)}")
        print(f"  - Volume: {atoms_new.get_volume():.2f} Å³")
        
        print(f"MX2 method:")
        print(f"  - Atoms: {len(atoms_old)}")
        print(f"  - Volume: {atoms_old.get_volume():.2f} Å³")
        
        vol_diff = abs(atoms_new.get_volume() - atoms_old.get_volume())
        print(f"Volume difference: {vol_diff:.2f} Å³")
    
    print("\n✓ Structure building test completed!")
    print("\nOutput files:")
    print("  - test_mos2_surfacebuilder.xyz (new method)")
    if atoms_old is not None:
        print("  - test_mos2_mx2.xyz (old method)")

if __name__ == "__main__":
    main()
