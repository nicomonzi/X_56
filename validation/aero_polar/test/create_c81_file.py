#!/usr/bin/env python3
"""
Generate C81 airfoil file from NASTRAN aerodynamic coefficient data.
The C81 format is used by MBDyn for aerodynamic load calculations.
"""

import json
from pathlib import Path
import numpy as np

def create_c81_from_nastran_data(coeff_file, output_c81):
    """
    Create C81 file from NASTRAN coefficient data.
    
    Args:
        coeff_file: Path to JSON file with extracted coefficients
        output_c81: Path where C81 file should be written
    """
    
    # Read NASTRAN coefficients
    with open(coeff_file, 'r') as f:
        data = json.load(f)
    
    # Sort by angle
    angles = sorted([float(a) for a in data.keys()])
    
    # Extract CL, CD, CM (CZ = lift, CX = drag, CMY = pitch moment)
    # Note: CZ is normal force (lift-like), CX is drag-like
    cl_values = [data[str(angle)]['CZ'] for angle in angles]
    cd_values = [data[str(angle)]['CX'] for angle in angles]
    cm_values = [data[str(angle)]['CMY'] for angle in angles]
    
    # C81 format: header and then data for each angle
    # The format assumes 11 chord positions (0.0, 0.2, 0.3, ..., 1.0)
    # and multiple angles of attack
    
    with open(output_c81, 'w') as f:
        # Write header with profile name and chord positions
        f.write("X-56 NASTRAN POLAR DATA          FROM AEROELASTIC      \n")
        f.write("       1.     .20    .30    .40    .50    .60    .70    .75    .80\n")
        f.write("       .90    1.\n")
        
        # Write data for each angle
        for angle in angles:
            cl = cl_values[angles.index(angle)]
            cd = cd_values[angles.index(angle)]
            cm = cm_values[angles.index(angle)]
            
            # Format: angle on first line, then same value for all 11 chord positions
            f.write(f"{angle:6.1f} ")
            f.write(" ".join([f"{cl:8.5f}"] * 11) + "\n")
            
            # Drag coefficient line
            f.write("       ")
            f.write(" ".join([f"{cd:8.5f}"] * 11) + "\n")
            
            # Pitching moment coefficient line
            f.write("       ")
            f.write(" ".join([f"{cm:8.5f}"] * 11) + "\n")
    
    print(f"✓ C81 file created: {output_c81}")
    print(f"  Angles: {angles[0]:.1f}° to {angles[-1]:.1f}°")
    print(f"  Data points: {len(angles)}")


if __name__ == "__main__":
    test_dir = Path(__file__).resolve().parent
    mbdyn_include = test_dir.parent / "mbdyn/INCLUDE"
    
    coeff_file = test_dir / "nastran_coefficients.json"
    c81_file = mbdyn_include / "x56_nastran.c81"
    
    create_c81_from_nastran_data(coeff_file, c81_file)
