#!/usr/bin/env python3
"""
Check ASE units module to see available pressure units
"""

from ase import units
import inspect

print("Available units in ase.units:")
for attr in sorted(dir(units)):
    if not attr.startswith('_'):
        value = getattr(units, attr)
        if isinstance(value, (int, float)):
            print(f"  {attr}: {value}")

print("\nPressure-related units:")
pressure_units = [attr for attr in dir(units) if 'pa' in attr.lower() or 'bar' in attr.lower() or 'atm' in attr.lower()]
for unit in pressure_units:
    value = getattr(units, unit)
    print(f"  {unit}: {value}")
