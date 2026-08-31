"""RootCTrait: 3D root system architecture trait extraction from segmented CT volumes."""
__version__ = "1.0.0"

from .io_volume import load_volume
from .graph_extraction import prune_skeleton
from .root_decomposition import decompose_root_system
from .decontamination import decontaminate, keep_base_component
from .detection_hypocotyle import collar_and_hypocotyl
from .root_traits_full import compute_all_traits
