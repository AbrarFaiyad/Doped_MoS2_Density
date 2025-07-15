#!/usr/bin/env python3
"""
Test script for MoS2 NPT simulation

This script tests the NPT simulation setup without running the full simulation.
It validates imports, structure building, and calculator setup.
"""

import sys
import traceback

def test_imports():
    """Test all required imports."""
    print("Testing imports...")
    try:
        from ase import units
        from ase.io import read, write, Trajectory
        from ase.md.nptberendsen import NPTBerendsen
        from ase.md import MDLogger
        from fairchem.core import pretrained_mlip, FAIRChemCalculator
        import numpy as np
        print("✓ All ASE and Fairchem imports successful")
        return True
    except ImportError as e:
        print(f"✗ Import error: {e}")
        return False

def test_surface_builder():
    """Test surface builder import."""
    print("\nTesting surface builder import...")
    try:
        sys.path.append('/home/afaiyad/borgstore/Energy_Profile_ML_Pot_vs_DFT/energy_profile_calculator')
        from energy_profile_calculator.surfaces import SurfaceBuilder
        print("✓ Surface builder import successful")
        return True
    except ImportError as e:
        print(f"✗ Surface builder import error: {e}")
        return False

def test_simulation_class():
    """Test MoS2NPTSimulation class initialization."""
    print("\nTesting simulation class...")
    try:
        from mos2_npt_simulation import MoS2NPTSimulation
        
        # Test with small parameters for quick testing
        sim = MoS2NPTSimulation(
            layers=2,
            supercell_size=(2, 2),
            temperature=300,
            pressure=1.0,
            timestep=1.0
        )
        print("✓ Simulation class initialization successful")
        return True
    except Exception as e:
        print(f"✗ Simulation class error: {e}")
        traceback.print_exc()
        return False

def main():
    """Run all tests."""
    print("="*50)
    print("MoS2 NPT SIMULATION TESTS")
    print("="*50)
    
    tests = [
        test_imports,
        test_surface_builder,
        test_simulation_class
    ]
    
    results = []
    for test in tests:
        results.append(test())
    
    print("\n" + "="*50)
    print("TEST SUMMARY")
    print("="*50)
    
    passed = sum(results)
    total = len(results)
    
    print(f"Tests passed: {passed}/{total}")
    
    if passed == total:
        print("✓ All tests passed! Ready to run simulation.")
        return True
    else:
        print("✗ Some tests failed. Check errors above.")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
