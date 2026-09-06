#!/usr/bin/env python3
"""
Generate 3 C81 files with span-wise variations from NASTRAN aerodynamic data.
- c81_root.c81: Root section (inner wing) - higher lift coefficient
- c81_wing.c81: Mid-wing section (main wing) - nominal coefficients
- c81_winglet.c81: Winglet section (outer wing/winglet) - reduced lift coefficient
"""

import json
from pathlib import Path
import numpy as np

def create_spanwise_c81_files():
    """
    Create 3 C81 files with span-wise variations based on NASTRAN data.
    The variations represent typical span-wise aerodynamic changes.
    """
    
    test_dir = Path(__file__).resolve().parent
    mbdyn_include = test_dir.parent / "mbdyn/INCLUDE"
    
    # Read NASTRAN coefficients
    with open(test_dir / "nastran_coefficients.json", 'r') as f:
        data = json.load(f)
    
    # Sort by angle
    angles = sorted([float(a) for a in data.keys()])
    
    # Extract base coefficients
    base_cl = {angle: data[str(angle)]['CZ'] for angle in angles}
    base_cd = {angle: data[str(angle)]['CX'] for angle in angles}
    base_cm = {angle: data[str(angle)]['CMY'] for angle in angles}
    
    print("Generating 3 C81 files with span-wise variations...")
    print("=" * 70)
    
    # Configuration for three sections
    sections = {
        'root': {
            'file': mbdyn_include / 'x56_root.c81',
            'description': 'Root section (inner wing)',
            'cl_factor': 1.15,    # Higher lift at root
            'cd_factor': 1.10,    # Slightly higher drag
            'cm_factor': 1.05,    # Slightly higher moment
        },
        'wing': {
            'file': mbdyn_include / 'x56_wing.c81',
            'description': 'Mid-wing section (main wing)',
            'cl_factor': 1.00,    # Nominal coefficients
            'cd_factor': 1.00,
            'cm_factor': 1.00,
        },
        'winglet': {
            'file': mbdyn_include / 'x56_winglet.c81',
            'description': 'Winglet section (outer wing)',
            'cl_factor': 0.85,    # Reduced lift at winglet
            'cd_factor': 0.90,    # Slightly reduced drag
            'cm_factor': 0.95,    # Slightly reduced moment
        }
    }
    
    # Create each C81 file
    for section_name, config in sections.items():
        create_c81_file(
            filename=config['file'],
            angles=angles,
            base_cl=base_cl,
            base_cd=base_cd,
            base_cm=base_cm,
            cl_factor=config['cl_factor'],
            cd_factor=config['cd_factor'],
            cm_factor=config['cm_factor'],
            description=f"X-56 {section_name.upper()}"
        )
        
        print(f"✓ {section_name.upper():8s} - {config['file'].name}")
        print(f"           {config['description']}")
        print(f"           CL factor: {config['cl_factor']:.2f}, CD factor: {config['cd_factor']:.2f}, CM factor: {config['cm_factor']:.2f}")
        print()
    
    print("=" * 70)
    print("Created 3 C81 files for span-wise aerodynamic analysis")
    
    return sections

def create_c81_file(filename, angles, base_cl, base_cd, base_cm, 
                    cl_factor=1.0, cd_factor=1.0, cm_factor=1.0, description="AIRFOIL"):
    """
    Create a single C81 file with scaled coefficients.
    
    Args:
        filename: Output C81 filename
        angles: List of angle of attack values
        base_cl, base_cd, base_cm: Dictionaries of base coefficients by angle
        cl_factor, cd_factor, cm_factor: Scaling factors for each section
        description: Header description for the airfoil
    """
    
    with open(filename, 'w') as f:
        # Write header
        f.write(f"{description:40s} FROM NASTRAN\n")
        f.write("       1.     .20    .30    .40    .50    .60    .70    .75    .80\n")
        f.write("       .90    1.\n")
        
        # Write data for each angle with span-wise scaling
        for angle in angles:
            cl = base_cl[angle] * cl_factor
            cd = base_cd[angle] * cd_factor
            cm = base_cm[angle] * cm_factor
            
            # Format: angle on first line with 11 chord positions
            f.write(f"{angle:6.1f} ")
            f.write(" ".join([f"{cl:8.5f}"] * 11) + "\n")
            
            # Drag coefficient line
            f.write("       ")
            f.write(" ".join([f"{cd:8.5f}"] * 11) + "\n")
            
            # Pitching moment coefficient line
            f.write("       ")
            f.write(" ".join([f"{cm:8.5f}"] * 11) + "\n")
    
    # Print sample data
    print(f"  Sample data for {filename.name}:")
    for angle in sorted(angles)[::10]:  # Print every 10th angle
        cl = base_cl[angle] * cl_factor
        cd = base_cd[angle] * cd_factor
        cm = base_cm[angle] * cm_factor
        print(f"    α = {angle:6.1f}°: CL = {cl:8.5f}, CD = {cd:8.5f}, CM = {cm:8.5f}")


if __name__ == "__main__":
    sections = create_spanwise_c81_files()
