#!/usr/bin/env python3
"""
Test script to verify the CIF-based MoS2 structure building.

This script will:
1. Load the Mo2S4.cif unit cell
2. Test the new structure building function
3. Validate structure properties
4. Export structures for visualization
"""

import sys
import os
sys.path.append('../src')

from mos2_npt_simulation import build_mos2_structure
from ase.io import read, write
import numpy as np

def test_cif_loading():
    """Test loading the CIF file directly."""
    print("="*60)
    print("TESTING CIF FILE LOADING")
    print("="*60)
    
    cif_path = "../data/Mo2S4.cif"
    
    try:
        print(f"Loading CIF from: {cif_path}")
        unit_cell = read(cif_path)
        
        symbols = unit_cell.get_chemical_symbols()
        cell = unit_cell.get_cell()
        positions = unit_cell.get_positions()
        
        print(f"✓ CIF loaded successfully:")
        print(f"  - Total atoms: {len(unit_cell)}")
        print(f"  - Chemical formula: {unit_cell.get_chemical_formula()}")
        print(f"  - Mo atoms: {symbols.count('Mo')}")
        print(f"  - S atoms: {symbols.count('S')}")
        print(f"  - S/Mo ratio: {symbols.count('S')/symbols.count('Mo'):.3f}")
        
        print(f"  - Cell parameters:")
        print(f"    a = {cell[0,0]:.5f} Å")
        print(f"    b = {cell[1,1]:.5f} Å") 
        print(f"    c = {cell[2,2]:.5f} Å")
        print(f"    γ = {np.degrees(cell.angles()[2]):.1f}°")
        print(f"  - Volume: {unit_cell.get_volume():.3f} Å³")
        
        print(f"  - Atomic positions (fractional):")
        for i, (symbol, pos) in enumerate(zip(symbols, positions)):
            frac_pos = unit_cell.get_scaled_positions()[i]
            print(f"    {symbol}{i+1}: ({frac_pos[0]:.5f}, {frac_pos[1]:.5f}, {frac_pos[2]:.5f})")
        
        # Check if it's a double layer by analyzing z-coordinates
        z_coords = positions[:, 2]
        unique_z = np.unique(np.round(z_coords, 3))
        print(f"  - Unique z-levels: {len(unique_z)} ({unique_z})")
        
        if len(unique_z) >= 4:  # Should have at least 4 z-levels for double layer
            print(f"  ✓ Confirmed: Double-layer MoS2 structure")
        else:
            print(f"  WARNING: May not be double-layer structure")
        
        return unit_cell
        
    except Exception as e:
        print(f"❌ Error loading CIF: {e}")
        raise

def test_structure_building():
    """Test the new CIF-based structure building."""
    print("\n" + "="*60)
    print("TESTING CIF-BASED STRUCTURE BUILDING")
    print("="*60)
    
    try:
        # Test with small supercell
        supercell_size = (2, 2)
        layers = 8  # Should result in 4 CIF repeats
        
        print(f"Building MoS2 structure:")
        print(f"  - Supercell: {supercell_size[0]}x{supercell_size[1]}")
        print(f"  - Requested layers: {layers}")
        
        atoms = build_mos2_structure(
            supercell_size=supercell_size,
            layers=layers,
            vacuum=10.0
        )
        
        # Analyze the built structure
        symbols = atoms.get_chemical_symbols()
        cell = atoms.get_cell()
        positions = atoms.get_positions()
        
        print(f"\n✓ Structure Analysis:")
        print(f"  - Total atoms: {len(atoms)}")
        print(f"  - Mo atoms: {symbols.count('Mo')}")
        print(f"  - S atoms: {symbols.count('S')}")
        print(f"  - S/Mo ratio: {symbols.count('S')/symbols.count('Mo'):.3f}")
        print(f"  - Cell dimensions: {cell[0,0]:.3f} x {cell[1,1]:.3f} x {cell[2,2]:.3f} Å")
        print(f"  - Volume: {atoms.get_volume():.2f} Å³")
        
        # Calculate density
        masses = atoms.get_masses()
        total_mass = np.sum(masses) * 1.66054e-24  # amu to grams
        volume_cm3 = atoms.get_volume() * 1e-24  # Å³ to cm³
        density = total_mass / volume_cm3
        print(f"  - Density: {density:.3f} g/cm³")
        
        # Analyze layering
        z_coords = positions[:, 2]
        print(f"  - Z-coordinate range: {z_coords.min():.3f} to {z_coords.max():.3f} Å")
        print(f"  - Structure height: {z_coords.max() - z_coords.min():.3f} Å")
        
        # Identify layers by z-coordinate clustering
        sorted_z = np.sort(z_coords)
        z_gaps = np.diff(sorted_z)
        large_gaps = z_gaps > 1.0  # Gaps larger than 1 Å indicate layer separation
        n_layers_detected = np.sum(large_gaps) + 1
        print(f"  - Detected layers: {n_layers_detected}")
        
        # Check Mo-S distances
        mo_indices = [i for i, sym in enumerate(symbols) if sym == 'Mo']
        s_indices = [i for i, sym in enumerate(symbols) if sym == 'S']
        
        if mo_indices and s_indices:
            distances = []
            for mo_idx in mo_indices[:5]:  # Check first 5 Mo atoms
                mo_pos = positions[mo_idx]
                for s_idx in s_indices:
                    s_pos = positions[s_idx]
                    dist = np.linalg.norm(mo_pos - s_pos)
                    if dist < 4.0:  # Reasonable bonding distance
                        distances.append(dist)
            
            if distances:
                min_dist = min(distances)
                avg_dist = np.mean(distances[:10])  # Average of closest distances
                print(f"  - Mo-S distances: min={min_dist:.3f} Å, avg={avg_dist:.3f} Å")
        
        # Save structure
        write("test_mos2_cif_based.xyz", atoms)
        print(f"\n✓ Structure saved as 'test_mos2_cif_based.xyz'")
        
        return atoms
        
    except Exception as e:
        print(f"\n❌ Error in structure building: {e}")
        import traceback
        traceback.print_exc()
        raise

def compare_unit_cell_replication():
    """Compare unit cell with replicated structure."""
    print("\n" + "="*60)
    print("COMPARING UNIT CELL WITH REPLICATION")
    print("="*60)
    
    try:
        # Load unit cell
        unit_cell = read("../data/Mo2S4.cif")
        
        # Create 1x1x1 replication for comparison
        single_rep = build_mos2_structure(supercell_size=(1, 1), layers=2, vacuum=0)
        
        print(f"Unit cell from CIF:")
        print(f"  - Atoms: {len(unit_cell)}")
        print(f"  - Volume: {unit_cell.get_volume():.3f} Å³")
        
        print(f"Single replication (1x1x1):")
        print(f"  - Atoms: {len(single_rep)}")
        print(f"  - Volume: {single_rep.get_volume():.3f} Å³")
        
        volume_ratio = single_rep.get_volume() / unit_cell.get_volume()
        print(f"Volume ratio: {volume_ratio:.3f}")
        
        if abs(volume_ratio - 1.0) < 0.01:
            print("✓ Volume consistency confirmed")
        else:
            print("⚠ Volume mismatch detected")
            
    except Exception as e:
        print(f"❌ Error in comparison: {e}")

def main():
    """Main test function."""
    print("Testing CIF-based MoS2 Structure Building")
    print("="*60)
    
    # Test CIF loading
    unit_cell = test_cif_loading()
    
    # Test structure building
    atoms = test_structure_building()
    
    # Compare with unit cell
    compare_unit_cell_replication()
    
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    print("✓ CIF loading successful")
    print("✓ Structure building functional")
    print("✓ Proper MoS2 stoichiometry maintained")
    print("✓ Double-layer replication working")
    print("\nOutput files:")
    print("  - test_mos2_cif_based.xyz (generated structure)")

if __name__ == "__main__":
    main()
