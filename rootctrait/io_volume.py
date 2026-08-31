"""Loading 3D volumes for the root phenotyping pipeline.

Reads an already segmented binary mask from various formats and returns a 3D
numpy array. The axis convention expected by the rest of the pipeline is
col0 = Y (depth), col1 = X, col2 = Z.

Recognized formats (by extension):
  .mat            MATLAB v7.3 (HDF5, via h5py) or older 5/6/7 (via scipy.io)
  .tif / .tiff    3D TIFF stack (via tifffile)
  .npy / .npz     NumPy array
  .nii / .nii.gz  NIfTI (via nibabel)

tifffile and nibabel are imported only if the corresponding format is actually
used; they are therefore not required to process .mat files only.

Orientation: each format stores axes its own way. If the loaded volume is not in
the order col0 = depth, provide axis_order to permute the axes (for example
axis_order=(2, 1, 0)). Always check the result on the HTML figure produced by the
pipeline: the depth axis must point downward.

Expected input: a binary (or near-binary) foreground mask. The pipeline then
binarizes with V > 0.5 * max(V). A multi-label segmentation must therefore be
reduced to the root label beforehand.
"""
import os
import numpy as np
import h5py
from scipy.io import loadmat


def _load_mat(path, transpose_v73=True):
    """Load a .mat file.

    v7.3 (HDF5): h5py returns the axes reversed (MATLAB is column-major, h5py
    reads row-major), hence the transpose to recover col0 = depth.
    Older format (scipy.io.loadmat): direct order, no transpose.
    """
    try:
        with h5py.File(path, 'r') as f:
            # pick the largest 3D array among the variables (robust when the file
            # holds several arrays; generic, no naming assumption)
            cand = [(k, f[k].shape) for k in f.keys()
                    if not k.startswith('#') and getattr(f[k], 'ndim', 0) == 3]
            if not cand:
                raise KeyError("No 3D variable in %s" % path)
            key = max(cand, key=lambda kv: int(np.prod(kv[1])))[0]
            arr = np.array(f[key])
            return arr.T if transpose_v73 else arr
    except OSError:
        d = loadmat(path)
        cand = [(k, np.asarray(v)) for k, v in d.items()
                if not k.startswith('__') and np.asarray(v).ndim == 3]
        if not cand:
            raise KeyError("No 3D variable in %s" % path)
        return max(cand, key=lambda kv: kv[1].size)[1]


def _load_npz(path):
    with np.load(path) as d:
        keys = list(d.keys())
        if not keys:
            raise KeyError("Empty .npz archive: %s" % path)
        return np.asarray(d[keys[0]])


def load_volume(path, axis_order=None, transpose_v73=True):
    """Load a 3D volume from .mat, .tif/.tiff, .npy/.npz or .nii/.nii.gz.

    Parameters
    ----------
    path : str
        File path. The format is inferred from the extension.
    axis_order : tuple or None
        Axis permutation applied after loading, to obtain col0 = depth
        (e.g. (2, 1, 0)). None = no permutation.
        For a v7.3 .mat, the internal transpose is handled separately;
        axis_order applies on top if provided (leave None for .mat files
        already processed with this pipeline).
    transpose_v73 : bool
        Internal transpose of v7.3 .mat files. Leave True unless you have a
        verified special case.

    Returns
    -------
    numpy.ndarray, 3D
    """
    low = path.lower()
    ext = '.nii.gz' if low.endswith('.nii.gz') else os.path.splitext(path)[1].lower()

    if ext == '.mat':
        V = _load_mat(path, transpose_v73=transpose_v73)
    elif ext in ('.tif', '.tiff'):
        import tifffile
        V = np.asarray(tifffile.imread(path))
    elif ext == '.npy':
        V = np.asarray(np.load(path))
    elif ext == '.npz':
        V = _load_npz(path)
    elif ext in ('.nii', '.nii.gz'):
        import nibabel as nib
        V = np.asarray(nib.load(path).get_fdata())
    else:
        raise ValueError(
            "Unrecognized format: '%s'. Supported formats: "
            ".mat, .tif/.tiff, .npy/.npz, .nii/.nii.gz" % ext)

    V = np.asarray(V)

    if axis_order is not None:
        if len(axis_order) != V.ndim:
            raise ValueError(
                "axis_order %s incompatible with a %dD volume (shape %s)"
                % (axis_order, V.ndim, V.shape))
        V = np.transpose(V, axes=axis_order)

    if V.ndim != 3:
        raise ValueError(
            "Expected a 3D volume, got %dD (shape %s) from %s. "
            "Use axis_order or check the file contents."
            % (V.ndim, V.shape, path))

    return V
