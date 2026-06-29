"""
This script generates the Aetherion dependency graph from CMakeLists.txt

Parses CMakeLists.txt to extract:
  - Executables (add_executable)
  - Libraries (add_library)
  - Dependencies (target_link_libraries)
  - Include directories (target_include_directories) [TODO: actually implement this]

Outputs:
  - DEPENDENCY_GRAPH.csv (edges)
  - DEPENDENCY_GRAPH_NODES.csv (node metadata)
  - Updates Mermaid diagram in DEPENDENCY_GRAPH.md 

Usage:
  python3 scripts/generate_dependency_graph.py

"""

# ============================================================================
# imports
# ============================================================================
import re
import csv
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional
from dataclasses import dataclass
import sys 

# ============================================================================
# Data Structures
# ============================================================================

@dataclass
class Target:
    """Represents the CMake target."""
    name: str
    target_type: str  # 'executable', 'library', 'interface'
    category: str    # 'MainApp', 'SimulationApp', 'Module', etc.
    description: str
    file_path: str
    deps: List[str]  # Raw dependency names
    is_resolved: bool = False  # unused but i'll probably need it later

@dataclass
class Edge:
    """Represents a dependency relationship."""
    source: str
    target: str
    dep_type: str      # 'external', 'internal'
    category: str      # 'Framework', 'Library', 'Module', etc.
    required: bool

# ============================================================================
# CMake Parser
# ============================================================================

class CMakeParser:
    """Parses CMakeLists.txt and extracts targets and dependencies."""

    def __init__(self, cmake_file: Path):
        self.cmake_file = cmake_file
        self.content = cmake_file.read_text()
        self.content_no_comments = self._strip_comments(self.content)
        self.targets: Dict[str, Target] = {}
        self.external_libs: Set[str] = set()

    @staticmethod
    def _strip_comments(content_raw: str) -> str:
        """Remove CMake comments (#...) from content."""
        lines = content_raw.split('\n')
        result = []
        for line_text in lines:
            # note: This may break strings with a # in them, but we won't be using that anyways.
            hash_pos = line_text.find('#')
            if hash_pos != -1:
                line_text = line_text[:hash_pos]
            result.append(line_text)
        return '\n'.join(result)

    @staticmethod
    def _find_matching_paren(content: str, start_idx: int) -> int:
        """Find the index of the closing parenthesis matching the one at start."""
        paren_depth = 0
        for idx in range(start_idx, len(content)):
            char = content[idx]
            if char == '(':
                paren_depth += 1
            elif char == ')':
                paren_depth -= 1
                if paren_depth == 0:
                    return idx
        return -1  # sad trombone

    def parse(self) -> Dict[str, Target]:
        """Parse the CMakeLists.txt file."""
        self._extract_targets()
        self._classify_dependencies()
        return self.targets

    def _extract_targets(self):
        """Extract all add_executable and add_library calls."""
        # Extract executables using parenthesis matching
        searchPos = 0
        while True:
            searchPos = self.content_no_comments.find('add_executable(', searchPos)
            if searchPos == -1:
                break
            paren_start = self.content_no_comments.find('(', searchPos)
            paren_end = self._find_matching_paren(self.content_no_comments, paren_start)
            if paren_end == -1:
                searchPos += 1
                continue
            
            block = self.content_no_comments[paren_start + 1:paren_end]
            tokens = re.split(r'\s+', block.strip())
            if tokens:
                targetName = tokens[0]
                self.targets[targetName] = Target(
                    name=targetName,
                    target_type='executable',
                    category='MainApp' if targetName == 'blackhole-sim' else 'SimulationApp',
                    description='',
                    file_path='src/',
                    deps=[]
                )
            searchPos = paren_end + 1
    
        # Extract libraries
        lib_searchPos = 0
        while True:
            lib_searchPos = self.content_no_comments.find('add_library(', lib_searchPos)
            if lib_searchPos == -1:
                break
            paren_start = self.content_no_comments.find('(', lib_searchPos)
            paren_end = self._find_matching_paren(self.content_no_comments, paren_start)
            if paren_end == -1:
                lib_searchPos += 1
                continue
            
            block = self.content_no_comments[paren_start + 1:paren_end]
            tokens = re.split(r'\s+', block.strip())
            if len(tokens) >= 2:
                libName = tokens[0]
                libType = tokens[1]  # STATIC, SHARED, INTERFACE
                self.targets[libName] = Target(
                    name=libName,
                    target_type=libType.lower(),
                    category='Library' if libType == 'STATIC' else 'Interface',
                    description='',
                    file_path='external/' if libName == 'aetherion_imgui' else 'src/',
                    deps=[]
                )
            lib_searchPos = paren_end + 1

    def _classify_dependencies(self):
        """Extract target_link_libraries for each target."""
        link_searchPos = 0
        while True:
            link_searchPos = self.content_no_comments.find('target_link_libraries(', link_searchPos)
            if link_searchPos == -1:
                break
            paren_start = self.content_no_comments.find('(', link_searchPos)
            paren_end = self._find_matching_paren(self.content_no_comments, paren_start)
            if paren_end == -1:
                link_searchPos += 1
                continue
            
            block = self.content_no_comments[paren_start + 1:paren_end]
            tokens = [t for t in re.split(r'\s+', block.strip()) if t]
            
            if tokens:
                targetName = tokens[0]
                if targetName in self.targets:
                    # Filter out keywords (keywords list could be a constant but whatever)
                    deps = [t for t in tokens[1:] if t not in ('PRIVATE', 'PUBLIC', 'INTERFACE')]
                    self.targets[targetName].deps = deps
                    
                    # Classify as external
                    for dep_name in deps:
                        # TODO: maybe extract this to a helper function??
                        if '::' in dep_name or dep_name.startswith('Qt') or dep_name.startswith('GLEW') or dep_name.startswith('OpenGL'):
                            self.external_libs.add(dep_name)
            
            link_searchPos = paren_end + 1

# ============================================================================
# Graph Generation
# ==========================================================================
def generate_edges(targets: Dict[str, Target], external_libs: Set[str]) -> List[Edge]:
    """Generate edge list from targets and dependencies."""
    edges: List[Edge] = []
    
    for targetName, targetObj in targets.items():
        for depName in targetObj.deps:
            # Classify dependency
            isInternal = depName in targets
            depType = 'internal' if isInternal else 'external'
            
            # Categorize (this could probably be extracted to a function)
            if '::' in depName:
                category = 'Library'
            elif depName.startswith('Qt'):
                category = 'Framework'
            elif depName.startswith('GLEW') or depName.startswith('OpenGL'):
                category = 'Library'
            elif isInternal:
                category = 'Module'
            else:
                category = 'Library'  # default to library ig
            
            edges.append(Edge(
                source=targetName,
                target=depName,
                dep_type=depType,
                category=category,
                required=True
            ))
    
    return edges

def generate_nodes_metadata() -> Dict[str, Tuple[str, str, str, str]]:
    """
    Return manually curated node metadata.
    (name, type, category, description, filepath)
    """
    return {
        'blackhole-sim': ('Executable', 'MainApp', 'Qt-based launcher and unified application container', 'src/QT-LAUNCHER/'),
        'blackhole-2D': ('Executable', 'SimulationApp', 'Standalone 2D black hole ray-tracing simulator', 'src/2D/BlackHole2D.cpp'),
        'blackhole-3D': ('Executable', 'SimulationApp', 'Standalone 3D GPU-accelerated black hole renderer', 'src/3D/BlackHole3D.cpp'),
        'physics-regression-tests': ('Executable', 'Testing', 'Physics validation test suite', 'tests/physics_regression_tests.cpp'),
        'aetherion_imgui': ('Library', 'ImGUI', 'Static library bundling Dear ImGui + ImGui-SFML + OpenGL3 backend', 'external/imgui/'),
        '2D-core': ('Module', 'Physics', 'Core types and body visual structures', 'src/2D/2D-core/'),
        '2D-physics': ('Module', 'Physics', 'Physics solvers (Schwarzschild/Kerr geodesics; RK4/RKF45 integrators)', 'src/2D/2D-physics/'),
        '2D-simulation': ('Module', 'Physics', 'Central simulation manager orchestrating photon/body physics', 'src/2D/2D-simulation/'),
        '2D-rendering': ('Module', 'Rendering', 'SFML-based 2D camera and rendering layer', 'src/2D/2D-rendering/'),
        '2D-visualization': ('Module', 'Rendering', 'Ray/orbit visualization and visual data preprocessing', 'src/2D/2D-visualization/'),
        '2D-ui': ('Module', 'UI', '2D simulation input controls and keybind handling', 'src/2D/2D-ui/'),
        '2D-utils': ('Module', 'Utilities', 'Configuration presets and utilities', 'src/2D/2D-utils/'),
        'simulation_2d_widget': ('Module', 'UI', 'Qt widget embedding SFML 2D simulation canvas', 'src/QT-LAUNCHER/simulation_2d_widget.cpp'),
        'simulation_3d_widget': ('Module', 'UI', 'Qt widget embedding SFML 3D simulation canvas', 'src/QT-LAUNCHER/simulation_3d_widget.cpp'),
        'mainwindow': ('Module', 'UI', 'Qt main window orchestrating simulation widgets and dialogs', 'src/QT-LAUNCHER/mainwindow.cpp'),
        'sfml_canvas': ('Module', 'UI', 'Abstract Qt+SFML integration bridge for embedding SFML in Qt', 'src/QT-LAUNCHER/sfml_canvas.cpp'),
        'keybindbuttonwidget': ('Module', 'UI', 'Interactive key remapping widget for action bindings', 'src/QT-LAUNCHER/keybindbuttonwidget.cpp'),
        'custom_bh_dialog': ('Module', 'UI', 'Custom black hole parameter input dialog', 'src/QT-LAUNCHER/custom_bh_dialog.cpp'),
        'updater': ('Module', 'Tools', 'GitHub Releases-based auto-update checker', 'src/QT-LAUNCHER/updater.cpp'),
        'bh_preview_widget': ('Module', 'UI', 'Object library preview and visual asset browser', 'src/QT-LAUNCHER/bh_preview_widget.cpp'),
    }

# ============================================================================
# CSV Export                                                                 =
# ============================================================================

def write_edges_csv(edges: List[Edge], output_path: Path):
    """Write edges to DEPENDENCY_GRAPH.csv."""
    with open(output_path, 'w', newline='') as outfile:
        fieldz = ['Source', 'Target', 'Type', 'Category', 'Required']
        writer = csv.DictWriter(outfile, fieldnames=fieldz)
        writer.writeheader()
        for edge in edges:
            writer.writerow({
                'Source': edge.source,
                'Target': edge.target,
                'Type': edge.dep_type,
                'Category': edge.category,
                'Required': str(edge.required).upper()
            })

def write_nodes_csv(targets: Dict[str, Target], external_libs: Set[str], output_path: Path):
    """Write nodes to DEPENDENCY_GRAPH_NODES.csv."""
    metadata = generate_nodes_metadata()
    
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['Node', 'Type', 'Category', 'Description', 'FilePath'])
        writer.writeheader()
        
        # Internal targets
        for name, target in targets.items():
            meta = metadata.get(name, (target.target_type.capitalize(), target.category, '', target.file_path))
            writer.writerow({
                'Node': name,
                'Type': meta[0],
                'Category': meta[1],
                'Description': meta[2],
                'FilePath': meta[3]
            })
        
        # External libraries
        for lib in sorted(external_libs):
            writer.writerow({
                'Node': lib,
                'Type': 'External',
                'Category': 'Library' if '::' in lib else 'Framework',
                'Description': '',
                'FilePath': 'N/A'
            })

# ============================================================================
# Main
# ============================================================================

def main():
    local_dir = Path(__file__).resolve().parent
    cmakeFile = local_dir / 'CMakeLists.txt'
    
    if not cmakeFile.exists():
        print(f"CMakeLists.txt not found at {cmakeFile}")
        return 1
    
    print(f"📖 Parsing {cmakeFile}...")
    parser = CMakeParser(cmakeFile)
    targets = parser.parse()
    
    targetCount = len(targets)
    libCount = len(parser.external_libs)
    print(f"Found {targetCount} targets")
    print(f"Found {libCount} external libraries")
    
    # Generate edges
    edges = generate_edges(targets, parser.external_libs)
    print(f"Generated {len(edges)} dependency edges")
    
    # Write CSVs
    edgesCSV = local_dir / 'DEPENDENCY_GRAPH.csv'
    nodesCSV = local_dir / 'DEPENDENCY_GRAPH_NODES.csv'
    
    write_edges_csv(edges, edgesCSV)
    print(f"📝 Wrote {edgesCSV}")
    
    write_nodes_csv(targets, parser.external_libs, nodesCSV)
    print(f"📝 Wrote {nodesCSV}")
    
    print("\n✨ Dependency graph regenerated from CMakeLists.txt")
    return 0

if __name__ == '__main__':
    exit(main())

# Side note for future: this file could be refactored to use a more robust CMake parser library if one exists, or to handle more complex CMake constructs. 
# For now, it works for the specific structure of the Aetherion project.
# Eventually i'm going to get rid of this or just keep it in a local clone for myself, since the dependency graphs won't truly matter unless multiple contributors work on this overtime.
# For now, it's just a convenient way to visualize the project structure and dependencies.