#!/usr/bin/env python3
"""
FEM File Parser for MBDyn Modal Data Files

This parser reads MBDyn modal FEM files and allows querying node coordinates by ID.
"""

from typing import Dict, List, Tuple, Optional
import numpy as np

FEM_FILE = "blade-modes-reference.fem" 
GRID_POINT = np.array(
        [674, 785, 776, 767, 758, 749, 740, 1581, 1563, 1545, 1405]
    )

class FEMParser:
    """Parser for MBDyn FEM modal data files."""
    
    def __init__(self, filename: str):
        """
        Initialize the parser and load the FEM file.
        
        Args:
            filename: Path to the .fem file
        """
        self.filename = filename
        self.node_ids: List[int] = []
        self.coordinates: Dict[int, Tuple[float, float, float]] = {}
        self._parse_file()
    
    def _parse_file(self):
        """Parse the FEM file and extract node IDs and coordinates."""
        with open(self.filename, 'r') as f:
            lines = f.readlines()
        
        # Find the line numbers for each section
        node_list_start = None
        x_coords_start = None
        y_coords_start = None
        z_coords_start = None
        
        for i, line in enumerate(lines):
            if 'RECORD GROUP 2' in line and 'FINITE ELEMENT NODE LIST' in line:
                node_list_start = i + 1
            elif 'RECORD GROUP 5' in line and 'NODAL X COORDINATES' in line:
                x_coords_start = i + 1
            elif 'RECORD GROUP 6' in line and 'NODAL Y COORDINATES' in line:
                y_coords_start = i + 1
            elif 'RECORD GROUP 7' in line and 'NODAL Z COORDINATES' in line:
                z_coords_start = i + 1
        
        # Parse node IDs
        if node_list_start is not None:
            self.node_ids = self._parse_node_list(lines, node_list_start)
        
        # Parse coordinates
        if x_coords_start is not None and y_coords_start is not None and z_coords_start is not None:
            x_coords = self._parse_coordinate_section(lines, x_coords_start, y_coords_start)
            y_coords = self._parse_coordinate_section(lines, y_coords_start, z_coords_start)
            z_coords = self._parse_coordinate_section(lines, z_coords_start, len(lines))
            
            # Combine into coordinate dictionary
            for i, node_id in enumerate(self.node_ids):
                self.coordinates[node_id] = (x_coords[i], y_coords[i], z_coords[i])
    
    def _parse_node_list(self, lines: List[str], start_line: int) -> List[int]:
        """Parse the finite element node list."""
        node_ids = []
        for i in range(start_line, len(lines)):
            line = lines[i].strip()
            if line.startswith('**'):
                break
            # Parse integers from the line
            ids = [int(x) for x in line.split()]
            node_ids.extend(ids)
        return node_ids
    
    def _parse_coordinate_section(self, lines: List[str], start_line: int, end_marker: int) -> List[float]:
        """Parse a coordinate section (X, Y, or Z)."""
        coords = []
        for i in range(start_line, end_marker):
            line = lines[i].strip()
            if line.startswith('**'):
                break
            if line:
                # Parse floating point numbers in scientific notation
                try:
                    value = float(line)
                    coords.append(value)
                except ValueError:
                    continue
        return coords
    
    def get_coordinates(self, node_id: int) -> Optional[Tuple[float, float, float]]:
        return self.coordinates.get(node_id)
    
    def get_node_ids(self) -> List[int]:
        return self.node_ids.copy()
    
    def get_all_coordinates(self) -> Dict[int, Tuple[float, float, float]]:
        return self.coordinates.copy()
    
    def find_nodes_in_range(self, x_range: Optional[Tuple[float, float]] = None,
                           y_range: Optional[Tuple[float, float]] = None,
                           z_range: Optional[Tuple[float, float]] = None) -> List[int]:
        """
        Find nodes within specified coordinate ranges.
        """
        matching_nodes = []
        for node_id, (x, y, z) in self.coordinates.items():
            if x_range and not (x_range[0] <= x <= x_range[1]):
                continue
            if y_range and not (y_range[0] <= y <= y_range[1]):
                continue
            if z_range and not (z_range[0] <= z <= z_range[1]):
                continue
            matching_nodes.append(node_id)
        return matching_nodes
    
    def __repr__(self):
        return f"FEMParser('{self.filename}', nodes={len(self.node_ids)})"


def main():
    """Main function to demonstrate FEMParser usage.""" 
    fem = FEMParser(FEM_FILE)
    
    interface_nodes = GRID_POINT
    for node_id in interface_nodes:
        coords = fem.get_coordinates(node_id)
        if coords:
            print(f"{node_id:4d}: {coords[0]:10.8f}, {coords[1]:10.8f}, {coords[2]:10.8f}")

if __name__ == "__main__":
    main()
