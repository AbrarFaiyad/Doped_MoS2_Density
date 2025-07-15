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
from ase.build import mx2
from fairchem.core import pretrained_mlip, FAIRChemCalculator
from fairchem.core.units.mlip_unit.api.inference import InferenceSettings
import numpy as np
import time
import os
import sys

def build_mos2_structure(supercell_size=(10, 10), layers=8, a=3.18, vacuum=0):
    """
    Build a layered MoS2 structure using ASE's mx2 function.
    
    Args:
        supercell_size (tuple): Size of supercell in x and y directions
        layers (int): Number of MoS2 layers
        a (float): In-plane lattice parameter in Angstroms
        vacuum (float): Vacuum space above structure in Angstroms
    
    Returns:
        Atoms: ASE Atoms object containing the MoS2 structure
    """
    print(f"Building MoS2 structure using ASE's mx2 function:")
    print(f"  - Supercell: {supercell_size[0]}x{supercell_size[1]}")
    print(f"  - Layers: {layers}")
    print(f"  - Lattice parameter a: {a:.2f} Å")
    print(f"  - Vacuum: {vacuum:.2f} Å")
    
    # Build MoS2 structure using ASE's mx2 function
    # The thickness parameter controls the layer separation
    atoms = mx2(
        formula='MoS2',
        kind='2H',  # 2H polytype (most common)
        a=a,
        thickness=3.19,  # Standard MoS2 layer thickness
        size=(supercell_size[0], supercell_size[1], layers),
        vacuum=vacuum
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
            # Build MoS2 supercell with specified parameters using ASE's mx2 function
            print(f"Building {self.layers}-layer MoS2 with {self.supercell_size[0]}x{self.supercell_size[1]} supercells...")
            self.atoms = build_mos2_structure(
                supercell_size=self.supercell_size,
                layers=self.layers,
                a=3.18,  # MoS2 lattice parameter
                vacuum=0
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
            pressure_au = self.pressure * units.Pa  # Convert Pa to eV/Å³
            
            # Compressibility for MoS2 (approximate value in atomic units)
            # Using a typical value for layered materials
            compressibility_au = 4.57e-5 / units.bar  # Å³/eV
            
            # Initialize NPTBerendsen dynamics
            dyn = NPTBerendsen(
                self.atoms,
                timestep=self.timestep * units.fs,
                temperature_K=self.temperature,
                pressure_au=pressure_au,
                taut=100 * units.fs,  # Temperature coupling time (100 fs)
                taup=1000 * units.fs,  # Pressure coupling time (1 ps)
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
            
            # Custom function to log density
            def log_density():
                current_density = self.calculate_density()
                current_temp = self.estimate_temperature()
                current_volume = self.atoms.get_volume()
                step = dyn.nsteps
                
                with open("density_log.txt", "a") as f:
                    f.write(f"{step:8d} {current_temp:8.2f} {current_volume:12.3f} {current_density:8.4f}\n")
            
            # Initialize density log file
            with open("density_log.txt", "w") as f:
                f.write("# Step    Temp(K)   Volume(Å³)   Density(g/cm³)\n")
            
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
    
    # Create simulation instance
    sim = MoS2NPTSimulation(
        layers=8,
        supercell_size=(10, 10),
        temperature=300,  # K
        pressure=1.0,     # atm
        timestep=1.0      # fs
    )
    
    try:
        # Setup and run simulation
        sim.setup_calculator()
        sim.build_structure()
        
        # Run NPT simulation
        results = sim.run_npt_simulation(
            steps=50000,      # 50 ps simulation
            log_interval=100,
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
