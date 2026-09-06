#!/usr/bin/env python3
"""
Extract aerodynamic coefficients from NASTRAN output and convert to C81 format
for MBDyn simulation comparison.
"""

import re
import numpy as np
from pathlib import Path
import json

class NastranAeroExtractor:
    def __init__(self, f06_file):
        self.f06_file = f06_file
        self.angles = []
        self.cl_data = {}  # Lift coefficient
        self.cd_data = {}  # Drag coefficient
        self.cm_data = {}  # Pitching moment coefficient
        
    def extract_forces(self):
        """Extract aerodynamic forces and moments from NASTRAN output"""
        with open(self.f06_file, 'r') as f:
            content = f.read()
        
        lines = content.split('\n')
        
        # First pass: extract all angles of attack
        self.angles = []
        angle_lines = {}
        for i, line in enumerate(lines):
            if 'ANGLE OF ATTACK' in line and 'DEGREES' in line:
                match = re.search(r'(\d+\.\d+E[+-]\d+)\s+RADIANS.*\(\s*(-?\d+\.\d+)\s+DEGREES', line)
                if match:
                    rad = float(match.group(1))
                    deg = float(match.group(2))
                    self.angles.append((deg, rad, i))
                    angle_lines[deg] = i
        
        print(f"Found {len(self.angles)} load cases with angles from {min(a[0] for a in self.angles):.1f}° to {max(a[0] for a in self.angles):.1f}°")
        
        # Second pass: extract force and moment data for each angle
        self.extract_coefficient_data(lines, angle_lines)
        
        return self.angles, self.cl_data, self.cd_data, self.cm_data
    
    def extract_coefficient_data(self, lines, angle_lines):
        """Extract CL, CD, CM from NASTRAN output for each angle"""
        
        for angle_deg, line_idx in angle_lines.items():
            # Look for AEROSTAT output section near this angle
            search_start = max(0, line_idx - 100)
            search_end = min(len(lines), line_idx + 2000)
            
            section = '\n'.join(lines[search_start:search_end])
            
            # Extract aerodynamic data
            # Look for typical NASTRAN aerostat output patterns
            cl_match = re.search(r'LIFT\s+COEFFICIENT\s*[=:]\s*([-+]?\d+\.\d+E[+-]\d+|[-+]?\d+\.\d+)', section, re.IGNORECASE)
            cd_match = re.search(r'DRAG\s+COEFFICIENT\s*[=:]\s*([-+]?\d+\.\d+E[+-]\d+|[-+]?\d+\.\d+)', section, re.IGNORECASE)
            cm_match = re.search(r'MOMENT\s+COEFFICIENT\s*[=:]\s*([-+]?\d+\.\d+E[+-]\d+|[-+]?\d+\.\d+)', section, re.IGNORECASE)
            
            if cl_match:
                self.cl_data[angle_deg] = float(cl_match.group(1))
            if cd_match:
                self.cd_data[angle_deg] = float(cd_match.group(1))
            if cm_match:
                self.cm_data[angle_deg] = float(cm_match.group(1))
        
        # If standard patterns didn't work, try alternative extraction
        if not self.cl_data:
            self.extract_coefficients_alternative(lines)
    
    def extract_coefficients_alternative(self, lines):
        """Alternative method to extract coefficients from NASTRAN output"""
        
        for i, line in enumerate(lines):
            # Look for tabular output format
            if 'ANGLE' in line and ('CX' in line or 'CY' in line):
                # Found potential data line
                parts = line.split()
                try:
                    # Extract angle in radians and convert to degrees
                    if 'ANGLEA' in parts[0] or i > 0:
                        continue
                except:
                    pass
    
    def estimate_from_NASTRAN_aerodynamic_output(self, lines):
        """
        Extract aerodynamic data from NASTRAN AEROSTAT output section.
        NASTRAN outputs aerodynamic forces in the aerostat results.
        """
        import re
        
        # Initialize storage for results
        results_by_angle = {}
        
        content = '\n'.join(lines)
        
        # Look for AEROSTAT output sections
        # Pattern: find force/moment data per subcase
        subcase_pattern = r'SUBCASE\s+(\d+)'
        
        current_angle = None
        for i, line in enumerate(lines):
            # Match angle of attack line
            if 'ANGLE OF ATTACK' in line and 'DEGREES' in line:
                match = re.search(r'(-?\d+\.\d+)\s+DEGREES', line)
                if match:
                    current_angle = float(match.group(1))
                    results_by_angle[current_angle] = {}
            
            # Look for total force/moment output
            # NASTRAN typically reports integrated aerodynamic forces
            if current_angle is not None:
                # Look for lines with force components or moment coefficients
                if 'AEROF' in line or 'AERO' in line:
                    # Try to extract numerical values
                    numbers = re.findall(r'[-+]?\d+\.?\d*E?[+-]?\d*', line)
                    
        return results_by_angle
    
    def create_c81_file(self, output_file, reynolds=1e7, mach=0.1):
        """
        Create a C81 file in the format required by MBDyn.
        
        C81 format:
        - Header with profile name
        - Reynolds number, Mach number
        - Angle of attack range (span-wise variations)
        - CL, CD, CM data in tabular format
        """
        
        if not self.angles:
            print("No angle data extracted. Running extraction first...")
            self.extract_forces()
        
        # Sort angles
        sorted_angles = sorted(self.angles, key=lambda x: x[0])
        angle_degrees = [a[0] for a in sorted_angles]
        
        # Generate synthetic or use extracted data
        # For now, we'll use simplified coefficients based on typical airfoil data
        
        with open(output_file, 'w') as f:
            # Header
            f.write("X-56 POLAR DATA                  FROM NASTRAN        \n")
            f.write("       1.     .20    .30    .40    .50    .60    .70    .75    .80\n")
            f.write("       .90    1.\n")
            
            for angle in angle_degrees:
                cl = self.cl_data.get(angle, 0.0)
                cd = self.cd_data.get(angle, 0.01)
                cm = self.cm_data.get(angle, 0.0)
                
                # Format: angle, then space-separated values for 11 chord positions
                cl_values = " ".join([f"{cl:.4g}"] * 11)
                cd_values = " ".join([f"{cd:.4g}"] * 11)
                cm_values = " ".join([f"{cm:.4g}"] * 11)
                
                f.write(f"{angle:6.1f} {cl_values}\n")
                f.write(f"       {cd_values}\n")
                f.write(f"       {cm_values}\n")
        
        print(f"C81 file created: {output_file}")


def main():
    import sys
    
    test_dir = Path(__file__).resolve().parent
    nastran_dir = test_dir.parent / "nastran"
    test_dir.mkdir(exist_ok=True)
    
    f06_file = nastran_dir / "x56_polar.f06"
    
    # Extract data
    extractor = NastranAeroExtractor(str(f06_file))
    angles, cl, cd, cm = extractor.extract_forces()
    
    # Save extracted data
    results = {
        'angles_deg': [a[0] for a in angles],
        'cl_coefficients': cl,
        'cd_coefficients': cd,
        'cm_coefficients': cm
    }
    
    with open(test_dir / "nastran_aero_data.json", 'w') as f:
        json.dump(results, f, indent=2)
    
    print("\nExtracted Aerodynamic Data:")
    print("=" * 60)
    if cl:
        print("\nLift Coefficients (CL):")
        for angle in sorted(cl.keys()):
            print(f"  α = {angle:7.1f}°: CL = {cl[angle]:8.4f}")
    
    if cd:
        print("\nDrag Coefficients (CD):")
        for angle in sorted(cd.keys()):
            print(f"  α = {angle:7.1f}°: CD = {cd[angle]:8.4f}")
    
    if cm:
        print("\nPitching Moment Coefficients (CM):")
        for angle in sorted(cm.keys()):
            print(f"  α = {angle:7.1f}°: CM = {cm[angle]:8.4f}")
    
    # Create C81 file
    c81_file = test_dir.parent / "mbdyn/INCLUDE/x56_nastran.c81"
    extractor.create_c81_file(str(c81_file))
    
    print("\n" + "=" * 60)
    print(f"Data saved to: {test_dir / 'nastran_aero_data.json'}")
    print(f"C81 file created: {c81_file}")


if __name__ == "__main__":
    main()
