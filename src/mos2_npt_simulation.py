#!/usr/bin/env python3
"""
MoS2 NPT Simulation Script

This script creates a layered MoS2 system (8 layers of 10x10 supercells) 
and runs an NPT simulation to study density and structural properties.

Features:
- 8-layer MoS2 system with 10x10 supercells per layer
- Fairchem UMA-S-1.1 calculator for accurate ML potential
- NPTBerendsen ensemble simulation for pressure and volume control
- Density monitoring and analysis
- Comprehensive logging and trajectory saving

Author: A. Faiyad
Date: July 2025
"""

from ase import units, Atoms
from ase.io import read, write, Trajectory
from ase.md.nptberendsen import NPTBerendsen
from ase.md import MDLogger
from fairchem.core import pretrained_mlip, FAIRChemCalculator
from fairchem.core.units.mlip_unit.api.inference import InferenceSettings
import numpy as np
import time
import os
import sys

# Import SurfaceBuilder for proper MoS2 construction
from surfaces import SurfaceBuilder

def build_mos2_structure(supercell_size=(10, 10), layers=8, vacuum=0):
    """
    Build a layered MoS2 structure using the comprehensive SurfaceBuilder class.
    
    Args:
        supercell_size (tuple): Size of supercell in x and y directions
        layers (int): Number of MoS2 layers
        vacuum (float): Vacuum space above and below structure in Angstroms
    
    Returns:
        Atoms: ASE Atoms object containing the MoS2 structure
    """
    print(f"Building MoS2 structure using SurfaceBuilder:")
    print(f"  - Supercell: {supercell_size[0]}x{supercell_size[1]}")
    print(f"  - Layers: {layers}")
    print(f"  - Vacuum: {vacuum:.2f} Å")
    
    # Initialize SurfaceBuilder
    surface_builder = SurfaceBuilder()
    
    # Build MoS2 structure using the comprehensive 2D material builder
    # This uses proper TMD lattice parameters and structure
    atoms = surface_builder.build_2d_material(
        material='MoS2',
        size=supercell_size,
        vacuum=vacuum if vacuum > 0 else 0,  # Only add vacuum if requested
        layers=layers
    )
    
    # Set periodic boundary conditions
    atoms.set_pbc([True, True, True])
    
    # Get structure information
    symbols = atoms.get_chemical_symbols()
    cell = atoms.get_cell()
    
    print(f"✓ Structure created with {len(atoms)} atoms")
    print(f"  - Mo atoms: {symbols.count('Mo')}")
    print(f"  - S atoms: {symbols.count('S')}")
    print(f"  - Cell: {cell[0,0]:.2f} x {cell[1,1]:.2f} x {cell[2,2]:.2f} Å")
    print(f"  - Volume: {atoms.get_volume():.2f} Å³")
    
    # Validate proper MoS2 stoichiometry
    mo_count = symbols.count('Mo')
    s_count = symbols.count('S')
    expected_ratio = 2.0  # S:Mo ratio should be 2:1
    actual_ratio = s_count / mo_count if mo_count > 0 else 0
    
    print(f"  - Stoichiometry check: S/Mo ratio = {actual_ratio:.3f} (expected: {expected_ratio:.3f})")
    
    if abs(actual_ratio - expected_ratio) > 0.01:
        print(f"  WARNING: Stoichiometry deviation detected!")
    else:
        print(f"  ✓ Proper MoS2 stoichiometry confirmed")
    
    return atoms

class MoS2NPTSimulation:
    """Class to handle MoS2 NPT simulation setup and execution."""
    
    def __init__(self, layers=8, supercell_size=(10, 10), 
                 temperature=300, pressure=1.0, timestep=1.0):
        """
        Initialize MoS2 NPT simulation parameters.
        
        Args:
            layers (int): Number of MoS2 layers
            supercell_size (tuple): Size of each layer (nx, ny)
            temperature (float): Target temperature in K
            pressure (float): Target pressure in atm
            timestep (float): MD timestep in fs
        """
        self.layers = layers
        self.supercell_size = supercell_size
        self.temperature = temperature
        self.pressure = pressure * 101325  # Convert atm to Pa
        self.timestep = timestep
        self.atoms = None
        self.calc = None
        
        print(f"Initializing MoS2 NPT simulation:")
        print(f"  - Layers: {self.layers}")
        print(f"  - Supercell size: {self.supercell_size[0]}x{self.supercell_size[1]}")
        print(f"  - Temperature: {self.temperature} K")
        print(f"  - Pressure: {pressure} atm ({self.pressure:.0f} Pa)")
        print(f"  - Timestep: {self.timestep} fs")
    
    def setup_calculator(self):
        """Initialize the Fairchem calculator."""
        print("\n" + "="*50)
        print("SETTING UP FAIRCHEM CALCULATOR")
        print("="*50)
        
        try:
            print("Loading UMA-S-1.1 model...")
            predictor = pretrained_mlip.get_predict_unit(
                "uma-s-1p1", 
                device="cuda", 
            )
            self.calc = FAIRChemCalculator(predictor, task_name="omc")
            print("✓ Fairchem calculator initialized successfully")
        except Exception as e:
            print(f"Error initializing calculator: {e}")
            raise
    
    def build_structure(self):
        """Build the layered MoS2 structure."""
        print("\n" + "="*50)
        print("BUILDING MoS2 STRUCTURE")
        print("="*50)
        
        try:
            # Build MoS2 supercell with specified parameters using SurfaceBuilder
            print(f"Building {self.layers}-layer MoS2 with {self.supercell_size[0]}x{self.supercell_size[1]} supercells...")
            self.atoms = build_mos2_structure(
                supercell_size=self.supercell_size,
                layers=self.layers,
                vacuum=0  # No vacuum for bulk simulation
            )
            
            # Assign calculator
            self.atoms.calc = self.calc
            
            # Print structure information
            positions = self.atoms.get_positions()
            cell = self.atoms.get_cell()
            symbols = self.atoms.get_chemical_symbols()
            
            print(f"✓ Structure built successfully:")
            print(f"  - Total atoms: {len(self.atoms)}")
            print(f"  - Mo atoms: {symbols.count('Mo')}")
            print(f"  - S atoms: {symbols.count('S')}")
            print(f"  - Cell dimensions: {cell[0,0]:.2f} x {cell[1,1]:.2f} x {cell[2,2]:.2f} Å")
            print(f"  - Cell volume: {self.atoms.get_volume():.2f} Å³")
            
            # Calculate initial density
            initial_density = self.calculate_density()
            print(f"  - Initial density: {initial_density:.3f} g/cm³")
            
            # Save initial structure
            write("initial_structure.xyz", self.atoms)
            print("✓ Initial structure saved as 'initial_structure.xyz'")
            
        except Exception as e:
            print(f"Error building structure: {e}")
            raise
    
    def calculate_density(self):
        """Calculate the density of the system in g/cm³."""
        if self.atoms is None:
            return 0.0
        
        # Get atomic masses
        masses = self.atoms.get_masses()
        total_mass = np.sum(masses)  # in amu
        
        # Convert to grams
        total_mass_g = total_mass * 1.66054e-24  # amu to grams
        
        # Get volume in cm³
        volume_ang3 = self.atoms.get_volume()
        volume_cm3 = volume_ang3 * 1e-24  # Å³ to cm³
        
        # Calculate density
        density = total_mass_g / volume_cm3
        return density
    
    def estimate_temperature(self):
        """Estimate temperature from kinetic energy."""
        try:
            ekin = self.atoms.get_kinetic_energy()
            dof = 3 * len(self.atoms)
            return ekin / (0.5 * dof * units.kB)
        except:
            return 0.0
    
    def run_npt_simulation(self, steps=50000, log_interval=100, traj_interval=100):
        """
        Run NPTBerendsen simulation.
        
        Args:
            steps (int): Number of simulation steps
            log_interval (int): Logging interval
            traj_interval (int): Trajectory saving interval
        """
        print("\n" + "="*50)
        print("RUNNING NPTBerendsen SIMULATION")
        print("="*50)
        
        print(f"Simulation parameters:")
        print(f"  - Steps: {steps}")
        print(f"  - Temperature: {self.temperature} K")
        print(f"  - Pressure: {self.pressure:.0f} Pa")
        print(f"  - Timestep: {self.timestep} fs")
        print(f"  - Log interval: {log_interval}")
        print(f"  - Trajectory interval: {traj_interval}")
        
        try:
            # Convert pressure to atomic units (eV/Å³)
            # Using the conversion: 1 Pa = 6.242e-12 eV/Å³
            pressure_au = self.pressure * 6.242e-12  # Convert Pa to eV/Å³
            
            # Compressibility for MoS2 (realistic value in atomic units)
            # Literature value for MoS2: ~1-5 × 10^-5 bar^-1
            # Conversion: 1 bar^-1 = 1.602e-6 Å³/eV
            compressibility_bar = 2.0e-5  # bar^-1 (conservative estimate for MoS2)
            compressibility_au = compressibility_bar * 1.602e-6  # Convert to Å³/eV
            
            # Additional safety checks for stability
            if compressibility_au > 1e-4:
                print(f"Warning: Compressibility {compressibility_au:.2e} may be too high for stability")
                compressibility_au = 1e-4  # Cap at safe value
            
            print(f"Pressure conversion:")
            print(f"  - Pressure in Pa: {self.pressure:.0f}")
            print(f"  - Pressure in atomic units: {pressure_au:.2e} eV/Å³")
            print(f"  - Compressibility (bar^-1): {compressibility_bar:.2e}")
            print(f"  - Compressibility (Å³/eV): {compressibility_au:.2e}")
            
            # Use more conservative NPT parameters for stability
            temperature_coupling = 500 * units.fs  # Longer coupling time
            pressure_coupling = 2000 * units.fs    # Longer coupling time
            
            print(f"NPT coupling parameters:")
            print(f"  - Temperature coupling: {temperature_coupling/units.fs:.0f} fs")
            print(f"  - Pressure coupling: {pressure_coupling/units.fs:.0f} fs")
            
            # Initialize NPTBerendsen dynamics
            dyn = NPTBerendsen(
                self.atoms,
                timestep=self.timestep * units.fs,
                temperature_K=self.temperature,
                pressure_au=pressure_au,
                taut=temperature_coupling,  # Conservative temperature coupling
                taup=pressure_coupling,     # Conservative pressure coupling  
                compressibility_au=compressibility_au,
                fixcm=True
            )
            
            # Set up logging
            logger = MDLogger(
                dyn, 
                self.atoms, 
                "npt_simulation.log", 
                header=True, 
                stress=True, 
                peratom=False
            )
            dyn.attach(logger, interval=log_interval)
            
            # Set up trajectory saving
            trajectory = Trajectory("npt_trajectory.traj", "w", self.atoms)
            dyn.attach(trajectory.write, interval=traj_interval)
            
            # Custom function to log density and check for instability
            def log_density():
                current_density = self.calculate_density()
                current_temp = self.estimate_temperature()
                current_volume = self.atoms.get_volume()
                step = dyn.nsteps
                
                # Check for system instability
                positions = self.atoms.get_positions()
                max_position = np.max(np.abs(positions))
                cell_volume = self.atoms.get_volume()
                
                # Safety checks
                if max_position > 1000:  # Positions > 1000 Å indicate instability
                    print(f"WARNING: Large atomic displacement detected at step {step}")
                    print(f"  Max position: {max_position:.2e} Å")
                
                if cell_volume > 1e6:  # Volume > 1M Å³ indicates runaway expansion
                    print(f"WARNING: Large volume expansion at step {step}")
                    print(f"  Volume: {cell_volume:.2e} Å³")
                
                with open("density_log.txt", "a") as f:
                    f.write(f"{step:8d} {current_temp:8.2f} {current_volume:12.3f} {current_density:8.4f} {max_position:12.3f}\n")
            
            # Initialize density log file with additional column
            with open("density_log.txt", "w") as f:
                f.write("# Step    Temp(K)   Volume(Å³)   Density(g/cm³)  MaxPos(Å)\n")
            
            dyn.attach(log_density, interval=log_interval)
            
            print("\nStarting NPT simulation...")
            start_time = time.time()
            
            # Run simulation
            dyn.run(steps=steps)
            
            simulation_time = time.time() - start_time
            trajectory.close()
            
            # Final statistics
            final_density = self.calculate_density()
            final_temp = self.estimate_temperature()
            final_volume = self.atoms.get_volume()
            
            print(f"\n✓ NPT simulation completed successfully!")
            print(f"  - Simulation time: {simulation_time:.1f} seconds")
            print(f"  - Final temperature: {final_temp:.1f} K")
            print(f"  - Final volume: {final_volume:.2f} Å³")
            print(f"  - Final density: {final_density:.4f} g/cm³")
            
            # Save final structure
            write("final_structure.xyz", self.atoms)
            print(f"  - Final structure saved as 'final_structure.xyz'")
            
            return {
                'simulation_time': simulation_time,
                'final_temperature': final_temp,
                'final_volume': final_volume,
                'final_density': final_density
            }
            
        except Exception as e:
            print(f"Error during NPT simulation: {e}")
            raise
    
    def analyze_results(self):
        """Analyze simulation results and generate summary."""
        print("\n" + "="*50)
        print("ANALYSIS SUMMARY")
        print("="*50)
        
        try:
            # Read density log
            data = np.loadtxt("density_log.txt", skiprows=1)
            
            if len(data) > 0:
                steps = data[:, 0]
                temperatures = data[:, 1]
                volumes = data[:, 2]
                densities = data[:, 3]
                
                print(f"Simulation statistics:")
                print(f"  - Data points: {len(data)}")
                print(f"  - Temperature range: {temperatures.min():.1f} - {temperatures.max():.1f} K")
                print(f"  - Temperature average: {temperatures.mean():.1f} ± {temperatures.std():.1f} K")
                print(f"  - Volume range: {volumes.min():.2f} - {volumes.max():.2f} Å³")
                print(f"  - Volume average: {volumes.mean():.2f} ± {volumes.std():.2f} Å³")
                print(f"  - Density range: {densities.min():.4f} - {densities.max():.4f} g/cm³")
                print(f"  - Density average: {densities.mean():.4f} ± {densities.std():.4f} g/cm³")
                
                # Save analysis
                with open("analysis_summary.txt", "w") as f:
                    f.write("MoS2 NPT Simulation Analysis\n")
                    f.write("="*30 + "\n")
                    f.write(f"System: {self.layers} layers, {self.supercell_size[0]}x{self.supercell_size[1]} supercells\n")
                    f.write(f"Total atoms: {len(self.atoms)}\n")
                    f.write(f"Target temperature: {self.temperature} K\n")
                    f.write(f"Target pressure: {self.pressure/101325:.1f} atm\n\n")
                    f.write("Results:\n")
                    f.write(f"Temperature: {temperatures.mean():.1f} ± {temperatures.std():.1f} K\n")
                    f.write(f"Volume: {volumes.mean():.2f} ± {volumes.std():.2f} Å³\n")
                    f.write(f"Density: {densities.mean():.4f} ± {densities.std():.4f} g/cm³\n")
                
                print("✓ Analysis saved to 'analysis_summary.txt'")
            else:
                print("Warning: No data found in density log")
                
        except Exception as e:
            print(f"Error during analysis: {e}")

def main():
    """Main function to run the MoS2 NPT simulation."""
    print("="*60)
    print("MoS2 NPT SIMULATION")
    print("="*60)
    
    # Create simulation instance with conservative parameters
    sim = MoS2NPTSimulation(
        layers=8,
        supercell_size=(10, 10),
        temperature=300,  # K
        pressure=1.0,     # atm
        timestep=0.5      # Reduced timestep for stability
    )
    
    try:
        # Setup and run simulation
        sim.setup_calculator()
        sim.build_structure()
        
        # Run NPT simulation with shorter initial run for testing
        results = sim.run_npt_simulation(
            steps=10000,      # Reduced steps for testing stability
            log_interval=50,  # More frequent logging
            traj_interval=100
        )
        
        # Analyze results
        sim.analyze_results()
        
        print("\n" + "="*60)
        print("SIMULATION COMPLETED SUCCESSFULLY!")
        print("="*60)
        print("\nOutput files generated:")
        print("  - initial_structure.xyz (initial structure)")
        print("  - final_structure.xyz (final structure)")
        print("  - npt_trajectory.traj (full trajectory)")
        print("  - npt_simulation.log (MD log)")
        print("  - density_log.txt (density evolution)")
        print("  - analysis_summary.txt (statistical analysis)")
        
    except Exception as e:
        print(f"\nSimulation failed with error: {e}")
        raise

if __name__ == "__main__":
    main()
