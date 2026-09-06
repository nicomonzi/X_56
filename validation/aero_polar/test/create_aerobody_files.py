#!/usr/bin/env python3
"""
Create modified aerobody files that use span-wise C81 data:
- aerobody_root.mbd: Uses x56_root.c81 for root sections (bodies 1-5)
- aerobody_wing.mbd: Uses x56_wing.c81 for mid-wing sections (bodies 6-19)
- aerobody_winglet.mbd: Uses x56_winglet.c81 for winglet sections (bodies 20, 40, 39)
- aerobody_control.mbd: Uses x56_wing.c81 for control surfaces (bodies 41-50)
"""

from pathlib import Path
import re

def create_spanwise_aerobody_files():
    """
    Create multiple aerobody files with span-wise C81 references.
    """
    
    x56_aero_dir = Path(__file__).resolve().parent.parent / "mbdyn/INCLUDE"
    
    # Read original aerobody
    with open(x56_aero_dir / "aerobody.mbd", 'r') as f:
        content = f.read()
    
    print("Creating aerobody files with span-wise C81 distribution...")
    print("=" * 70)
    
    # 1. Root section aerobody (bodies 1-5)
    root_aero = content
    root_aero = root_aero.replace(
        'c81 data: 12, "./naca0012.c81";',
        'c81 data: 1, "./x56_root.c81";'
    )
    
    # Remove wing and winglet sections, keep only bodies 1-5
    root_aero = keep_bodies_range(root_aero, [1, 2, 3, 4, 5])
    root_aero = root_aero.replace('3, c81, 12, jacobian, yes;', '3, c81, 1, jacobian, yes;')
    
    with open(x56_aero_dir / "aerobody_root.mbd", 'w') as f:
        f.write(root_aero)
    print("✓ aerobody_root.mbd (bodies 1-5, using x56_root.c81)")
    
    # 2. Wing section aerobody (bodies 6-19)
    wing_aero = content
    wing_aero = wing_aero.replace(
        'c81 data: 12, "./naca0012.c81";',
        'c81 data: 2, "./x56_wing.c81";'
    )
    wing_aero = keep_bodies_range(wing_aero, list(range(6, 20)))
    wing_aero = wing_aero.replace('3, c81, 12, jacobian, yes;', '3, c81, 2, jacobian, yes;')
    
    with open(x56_aero_dir / "aerobody_wing.mbd", 'w') as f:
        f.write(wing_aero)
    print("✓ aerobody_wing.mbd (bodies 6-19, using x56_wing.c81)")
    
    # 3. Winglet section aerobody (bodies 20, 40, 39)
    winglet_aero = content
    winglet_aero = winglet_aero.replace(
        'c81 data: 12, "./naca0012.c81";',
        'c81 data: 3, "./x56_winglet.c81";'
    )
    winglet_aero = keep_bodies_range(winglet_aero, [20, 40, 39])
    winglet_aero = winglet_aero.replace('3, c81, 12, jacobian, yes;', '3, c81, 3, jacobian, yes;')
    
    with open(x56_aero_dir / "aerobody_winglet.mbd", 'w') as f:
        f.write(winglet_aero)
    print("✓ aerobody_winglet.mbd (bodies 20, 40, 39, using x56_winglet.c81)")
    
    # 4. Control surfaces aerobody (bodies 41-50)
    control_aero = content
    control_aero = control_aero.replace(
        'c81 data: 12, "./naca0012.c81";',
        'c81 data: 2, "./x56_wing.c81";'
    )
    control_aero = keep_bodies_range(control_aero, list(range(41, 51)))
    control_aero = control_aero.replace('3, c81, 12, jacobian, yes;', '3, c81, 2, jacobian, yes;')
    
    with open(x56_aero_dir / "aerobody_control.mbd", 'w') as f:
        f.write(control_aero)
    print("✓ aerobody_control.mbd (bodies 41-50, using x56_wing.c81)")
    
    # 5. Combined aerobody file with all sections
    combined_aero = create_combined_aerobody(content, x56_aero_dir)
    with open(x56_aero_dir / "aerobody_combined.mbd", 'w') as f:
        f.write(combined_aero)
    print("✓ aerobody_combined.mbd (all bodies 1-50, with span-wise C81)")
    
    print("=" * 70)
    print("All aerobody files created successfully!")


def keep_bodies_range(content, body_ids):
    """
    Keep only specified aerodynamic bodies, remove others.
    """
    lines = content.split('\n')
    result = []
    skip_body = False
    
    for line in lines:
        # Check if this is an aerodynamic body definition
        match = re.match(r'\s*aerodynamic body:\s*(\d+),', line)
        
        if match:
            body_id = int(match.group(1))
            skip_body = body_id not in body_ids
        
        if not skip_body:
            result.append(line)
        elif skip_body and re.match(r'\s*aerodynamic body:', line):
            # Skip until next definition or end
            continue
        elif skip_body and (line.strip() == '' or line.strip().startswith('#')):
            # Keep comments and blank lines between definitions
            result.append(line)
    
    return '\n'.join(result)


def create_combined_aerobody(content, output_dir):
    """
    Create a combined aerobody file with all bodies using span-wise C81.
    Replaces the single c81 data definition with multiple definitions.
    """
    
    # Replace single c81 definition with three definitions
    header = content[:content.find('c81 data:')]
    rest = content[content.find('c81 data:'):]
    
    # Build new header with three c81 definitions
    new_header = header
    if 'c81 data:' in rest:
        # Find and replace the first c81 data line
        rest = rest.split('\n', 1)
        rest = rest[1] if len(rest) > 1 else ''
    
    combined = new_header
    combined += '    c81 data: 1, "./x56_root.c81";\n'
    combined += '    c81 data: 2, "./x56_wing.c81";\n'
    combined += '    c81 data: 3, "./x56_winglet.c81";\n\n'
    combined += rest
    
    # Replace all c81 references with appropriate ones based on body ID
    lines = combined.split('\n')
    result = []
    
    for i, line in enumerate(lines):
        match = re.match(r'\s*aerodynamic body:\s*(\d+),', line)
        if match:
            body_id = int(match.group(1))
            
            # Determine which C81 to use based on body ID
            if body_id in range(1, 6):  # Root sections
                c81_id = 1
            elif body_id in range(6, 20):  # Wing sections
                c81_id = 2
            elif body_id in [20, 40, 39]:  # Winglet sections
                c81_id = 3
            else:  # Control surfaces and others
                c81_id = 2
            
            result.append(line)
            
            # Find and replace the c81 reference in this body definition
            for j in range(i+1, min(i+15, len(lines))):
                if '3, c81, 12, jacobian, yes;' in lines[j]:
                    lines[j] = lines[j].replace('12', str(c81_id))
                    break
        else:
            result.append(line)
    
    return '\n'.join(result)


if __name__ == "__main__":
    create_spanwise_aerobody_files()
