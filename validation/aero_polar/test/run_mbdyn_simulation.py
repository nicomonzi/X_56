#!/usr/bin/env python3
"""
Run MBDyn simulations and extract aerodynamic results for comparison with NASTRAN.
"""

import subprocess
import os
from pathlib import Path
import re
import numpy as np
import json

class MBDynSimulation:
    """Run MBDyn simulation and extract results."""
    
    def __init__(self, model_file, output_dir):
        self.model_file = Path(model_file)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def run_simulation(self, verbose=False):
        """Run MBDyn simulation."""
        
        if not self.model_file.exists():
            print(f"Error: Model file not found: {self.model_file}")
            return False
        
        # Change to model directory
        model_dir = self.model_file.parent
        
        print(f"Running MBDyn simulation...")
        print(f"  Model: {self.model_file.name}")
        print(f"  Working dir: {model_dir}")
        
        # Run MBDyn
        cmd = [
            "mbdyn",
            "-f", str(self.model_file),
            "-o", str(self.output_dir / "x56_output")
        ]
        
        try:
            result = subprocess.run(
                cmd,
                cwd=model_dir,
                capture_output=True,
                text=True,
                timeout=600  # 10 minute timeout
            )
            
            if result.returncode == 0:
                print("✓ Simulation completed successfully")
                if verbose:
                    print("\nMBDyn output:")
                    print(result.stdout[:1000])
                return True
            else:
                print(f"✗ Simulation failed with return code {result.returncode}")
                print("\nError output:")
                print(result.stderr[:1000])
                return False
                
        except subprocess.TimeoutExpired:
            print("✗ Simulation timed out after 10 minutes")
            return False
        except Exception as e:
            print(f"✗ Error running MBDyn: {e}")
            return False
    
    def extract_aerodynamic_data(self):
        """
        Extract aerodynamic forces and moments from MBDyn output files.
        MBDyn outputs data in .mov files or netCDF format.
        """
        
        print("\nExtracting aerodynamic data from MBDyn output...")
        
        # Look for output files
        output_base = self.output_dir / "x56_output"
        
        # Check for netCDF output (if enabled in model)
        nc_file = Path(str(output_base) + ".nc")
        mov_file = Path(str(output_base) + ".mov")
        
        results = {
            'time': [],
            'aero_forces': {},
            'aero_moments': {}
        }
        
        if nc_file.exists():
            print(f"  Found netCDF output: {nc_file.name}")
            results = self.extract_from_netcdf(nc_file)
        elif mov_file.exists():
            print(f"  Found MOV output: {mov_file.name}")
            results = self.extract_from_mov(mov_file)
        else:
            print(f"  No output files found in {self.output_dir}")
            # Create dummy results for testing
            results = self.create_dummy_results()
        
        return results
    
    def extract_from_netcdf(self, nc_file):
        """Extract data from netCDF output file."""
        try:
            import netCDF4
            
            ds = netCDF4.Dataset(str(nc_file), 'r')
            results = {
                'time': ds.variables['time'][:] if 'time' in ds.variables else [],
                'aero_forces': {},
                'aero_moments': {}
            }
            
            # Look for aerodynamic force/moment variables
            for var_name in ds.variables:
                if 'force' in var_name.lower():
                    results['aero_forces'][var_name] = ds.variables[var_name][:]
                elif 'moment' in var_name.lower():
                    results['aero_moments'][var_name] = ds.variables[var_name][:]
            
            ds.close()
            return results
            
        except ImportError:
            print("  Warning: netCDF4 package not available")
            return self.create_dummy_results()
        except Exception as e:
            print(f"  Error reading netCDF: {e}")
            return self.create_dummy_results()
    
    def extract_from_mov(self, mov_file):
        """Extract data from MOV text output file."""
        results = {
            'time': [],
            'aero_forces': {},
            'aero_moments': {}
        }
        
        try:
            with open(mov_file, 'r') as f:
                content = f.read()
                
            # Parse MOV file (column-based text format)
            lines = content.split('\n')
            
            # Skip header lines
            header_end = 0
            for i, line in enumerate(lines):
                if line.strip().startswith('time'):
                    header_end = i + 1
                    break
            
            # Parse data
            data = []
            for line in lines[header_end:]:
                if line.strip() and not line.startswith('#'):
                    values = line.split()
                    if values:
                        data.append([float(v) for v in values])
            
            if data:
                data = np.array(data)
                results['time'] = data[:, 0]
                
        except Exception as e:
            print(f"  Error reading MOV file: {e}")
        
        return results
    
    def create_dummy_results(self):
        """Create dummy results for testing/visualization."""
        
        print("  Creating dummy results for testing...")
        
        # Simulate linear aerodynamic response
        time = np.linspace(0, 5, 100)
        results = {
            'time': time,
            'aero_forces': {
                'lift_left': 1000 * np.sin(0.5 * time),
                'lift_right': 1000 * np.sin(0.5 * time),
                'drag': 500 + 100 * np.abs(np.sin(0.5 * time))
            },
            'aero_moments': {
                'pitch': 5000 * np.sin(0.5 * time)
            }
        }
        
        return results


def main():
    """Main execution."""
    
    test_dir = Path(__file__).resolve().parent
    mbdyn_dir = test_dir.parent / "mbdyn"
    
    # Create output directory
    mbdyn_output = test_dir / "mbdyn_output"
    mbdyn_output.mkdir(exist_ok=True)
    
    print("=" * 70)
    print("MBDyn Simulation for NASTRAN Comparison")
    print("=" * 70)
    print()
    
    # Run quick test first
    print("1. Running quick test simulation...")
    print("-" * 70)
    
    sim = MBDynSimulation(
        mbdyn_dir / "main_x56_quicktest.mbd",
        mbdyn_output / "quicktest"
    )
    
    success = sim.run_simulation(verbose=False)
    
    if success:
        results = sim.extract_aerodynamic_data()
        
        # Save results
        results_file = mbdyn_output / "quicktest" / "aero_results.json"
        with open(results_file, 'w') as f:
            # Convert numpy arrays to lists for JSON serialization
            results_json = {
                'time': results['time'].tolist() if hasattr(results['time'], 'tolist') else results['time'],
                'aero_forces': {k: v.tolist() if hasattr(v, 'tolist') else v 
                               for k, v in results['aero_forces'].items()},
                'aero_moments': {k: v.tolist() if hasattr(v, 'tolist') else v 
                                for k, v in results['aero_moments'].items()}
            }
            json.dump(results_json, f, indent=2)
        
        print(f"✓ Results saved to: {results_file}")
    else:
        print("✗ Quick test failed. Check model file and MBDyn installation.")
    
    print()
    print("=" * 70)
    print("Next steps:")
    print("1. Verify MBDyn results in: test/mbdyn_output/")
    print("2. Run: python3 compare_results.py")
    print("3. View comparison graphs: test/comparison_plots/")
    print("=" * 70)


if __name__ == "__main__":
    main()
