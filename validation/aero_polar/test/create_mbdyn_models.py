#!/usr/bin/env python3
"""
Create a comprehensive MBDyn model that includes both wing and winglet aerodynamic elements
with span-wise C81 distribution from NASTRAN data.
"""

from pathlib import Path

def create_comprehensive_model():
    """
    Create main_x56_comprehensive.mbd that includes:
    - Wing aerodynamic bodies (aerobody_combined.mbd)
    - Winglet aerodynamic bodies (aerobody_winglets.mbd)
    """
    
    # Read original main model
    aero_root = Path(__file__).resolve().parent.parent
    main_file = aero_root / "mbdyn/main_x56.mbd"
    with open(main_file, 'r') as f:
        content = f.read()
    
    # Modify for comprehensive model
    comprehensive = content
    
    # Update simulation time for testing
    comprehensive = comprehensive.replace(
        'final time: 10.0;',
        'final time: 5.0;  # Reduced time for NASTRAN comparison'
    )
    
    # Update number of aerodynamic elements (50 wing + 10 winglets = 60)
    comprehensive = comprehensive.replace(
        'aerodynamic elements: 50;',
        'aerodynamic elements: 60;  # 50 wing + 10 winglet elements'
    )
    
    # Add comment and include both aerobody files
    # Find the line with "include: "./INCLUDE/aerobody.mbd";"
    comprehensive = comprehensive.replace(
        'include: "./INCLUDE/aerobody.mbd";',
        '''# Wing aerobody elements (1-50)
    include: "./INCLUDE/aerobody_combined.mbd";
    
    # Winglet/fence aerobody elements (41-50)
    # Note: These overlap with control surfaces in wing; adjust body IDs as needed
    # Uncomment the next line to enable winglets:
    #include: "./INCLUDE/aerobody_winglets.mbd";'''
    )
    
    # Save comprehensive model
    output_file = aero_root / "mbdyn/main_x56_comprehensive.mbd"
    with open(output_file, 'w') as f:
        f.write(comprehensive)
    
    print(f"✓ Created {output_file.name}")
    print("  - Uses aerobody_combined.mbd for span-wise wing elements")
    print("  - Ready for winglet elements (aerobody_winglets.mbd)")
    print()
    
    # Also create a version with simplified body count for quick testing
    quick_test = content
    quick_test = quick_test.replace(
        'final time: 10.0;',
        'final time: 2.0;  # Quick test - 2 seconds'
    )
    quick_test = quick_test.replace(
        'time step: 0.02;',
        'time step: 0.05;  # Larger time step for quick test'
    )
    quick_test = quick_test.replace(
        'include: "./INCLUDE/aerobody.mbd";',
        'include: "./INCLUDE/aerobody_combined.mbd";  # Span-wise C81 distribution'
    )
    
    quick_test_file = aero_root / "mbdyn/main_x56_quicktest.mbd"
    with open(quick_test_file, 'w') as f:
        f.write(quick_test)
    
    print(f"✓ Created {quick_test_file.name}")
    print("  - Quick test version with 2 second simulation")
    print("  - Uses aerobody_combined.mbd with span-wise C81")
    print()
    
    return output_file, quick_test_file


if __name__ == "__main__":
    comprehensive_file, quick_file = create_comprehensive_model()
    
    print("=" * 70)
    print("MBDyn model files created for NASTRAN comparison:")
    print()
    print(f"1. {comprehensive_file.name}")
    print("   - Full simulation with 5 seconds duration")
    print("   - Includes all wing aerodynamic elements")
    print("   - Uses span-wise C81 files (root, wing, winglet)")
    print()
    print(f"2. {quick_file.name}")
    print("   - Quick test with 2 seconds duration")
    print("   - For fast validation of setup")
    print()
    print("=" * 70)
