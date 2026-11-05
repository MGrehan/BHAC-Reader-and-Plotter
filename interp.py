import numpy as np
import time
from pykdtree.kdtree import KDTree
from numba import njit, prange, types
from typing import Union, List, Dict, Tuple

import numba


DEFAULT_KD_TREE_LEAFSIZE = 16
DEFAULT_KD_NEIGHBOURS = 4

# Define Numba signatures for better performance
sig_linear_interp_single = types.float32[:](
    types.int32[:,:],      # neighbor_indices
    types.float32[:,:],    # neighbor_distances
    types.float32[:],      # values
    types.int32            # k
)

sig_linear_interp_multi = types.float32[:,:](
    types.int32[:,:],      # neighbor_indices
    types.float32[:,:],    # neighbor_distances  
    types.float32[:,:],    # all_values
    types.int32,           # k
    types.int32            # n_vars
)


@njit(sig_linear_interp_single, parallel=True, cache=True, fastmath=True)
def _fast_linear_interpolate(neighbor_indices, neighbor_distances, values, k):
    """
    Perform fast linear interpolation using pre-computed nearest neighbors.
    Uses inverse distance weighting with Numba JIT compilation.
    """
    n_points = neighbor_indices.shape[0]
    result = np.empty(n_points, dtype=np.float32)
    
    for i in prange(n_points):
        if neighbor_distances[i, 0] < 1e-10:
            # Exact match with a source point
            result[i] = values[neighbor_indices[i, 0]]
        else:
            # Inverse distance weighting
            weights = 1.0 / (neighbor_distances[i, :k] + 1e-10)
            weight_sum = np.sum(weights)
            weighted_values = 0.0
            for j in range(k):
                weighted_values += values[neighbor_indices[i, j]] * weights[j]
            result[i] = weighted_values / weight_sum
    
    return result


@njit(sig_linear_interp_multi, parallel=True, cache=True, fastmath=True)
def _fast_interpolate_multiple_vars(neighbor_indices, neighbor_distances, all_values, k, n_vars):
    """
    Interpolate multiple variables at once to maximize cache reuse.
    This is the key optimization - process all variables in a single pass.
    """
    n_points = neighbor_indices.shape[0]
    results = np.empty((n_vars, n_points), dtype=np.float32)
    
    for i in prange(n_points):
        # Calculate weights once for all variables
        if neighbor_distances[i, 0] < 1e-10:
            # Exact match
            idx = neighbor_indices[i, 0]
            for v in range(n_vars):
                results[v, i] = all_values[v, idx]
        else:
            # Calculate inverse distance weights
            inv_distances = 1.0 / (neighbor_distances[i, :k] + 1e-10)
            weight_sum = np.sum(inv_distances)
            weights = inv_distances / weight_sum
            # Apply weights to all variables
            for v in range(n_vars):
                weighted_sum = 0.0
                for j in range(k):
                    idx = neighbor_indices[i, j]
                    weighted_sum += all_values[v, idx] * weights[j]
                results[v, i] = weighted_sum
                
    return results


def fast_interpolate_var_to_grid(data, var, Ngrid_x=2048, Ngrid_y=2048, 
                                x_range=None, y_range=None, k_neighbors=DEFAULT_KD_NEIGHBOURS,
                                leafsize=DEFAULT_KD_TREE_LEAFSIZE, verbose=False):
    """
    Fast interpolation function that replaces interpolate_var_to_grid with KD-tree optimization.
    
    Parameters:
    - data (dict): Dictionary containing 'center_x', 'center_y', and the variable data
    - var (str): Variable name to interpolate
    - Ngrid_x (int): Number of grid points in x direction
    - Ngrid_y (int): Number of grid points in y direction  
    - method (str): Interpolation method ('fast_kd', 'nearest', 'linear', 'cubic')
    - x_range (tuple): Optional x-axis bounds (xmin, xmax)
    - y_range (tuple): Optional y-axis bounds (ymin, ymax)
    - k_neighbors (int): Number of neighbors for KD-tree interpolation
    - leafsize (int): Leaf size for KD-tree construction
    
    Returns:
    - tuple: (grid_x, grid_y, interpolated_var)
    """
    if verbose:
        print('===============================')
        print(f"Started fast interpolating variable: {var}")
    start_time = time.time()
    
    center_x, center_y = data["center_x"], data["center_y"]
    
    # Create mask for spatial filtering
    mask = np.ones(center_x.shape, dtype=bool)
    
    if x_range is not None:
        x_mask = (center_x >= x_range[0]) & (center_x <= x_range[1])
        mask &= x_mask
        
    if y_range is not None:
        y_mask = (center_y >= y_range[0]) & (center_y <= y_range[1])
        mask &= y_mask
    
    # Filter data
    filtered_center_x = center_x[mask]
    filtered_center_y = center_y[mask]
    filtered_var_data = var[mask]
    
    # Create uniform grid
    grid_x_1d = np.linspace(filtered_center_x.min(), filtered_center_x.max(), Ngrid_x, dtype=np.float32)
    grid_y_1d = np.linspace(filtered_center_y.min(), filtered_center_y.max(), Ngrid_y, dtype=np.float32)
    grid_x, grid_y = np.meshgrid(grid_x_1d, grid_y_1d, indexing='ij')
    
    # Use fast KD-tree interpolation
    # Ensure data is contiguous and float32
    source_points = np.column_stack((
        np.ascontiguousarray(filtered_center_x, dtype=np.float32),
        np.ascontiguousarray(filtered_center_y, dtype=np.float32)
    ))
    
    # Build KD-tree
    tree = KDTree(source_points, leafsize=leafsize)
    
    # Query neighbors for grid points
    grid_points = np.column_stack((grid_x.ravel(), grid_y.ravel()))
    distances, indices = tree.query(grid_points, k=k_neighbors)
    
    # Prepare data for Numba function
    values = np.ascontiguousarray(filtered_var_data, dtype=np.float32)
    distances = distances.astype(np.float32)
    indices = indices.astype(np.int32)
    
    # Interpolate using fast Numba function
    interpolated_flat = _fast_linear_interpolate(indices, distances, values, k_neighbors)
    interpolated_var = interpolated_flat.reshape(grid_x.shape)
    
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    if verbose:
        print(f"Finished fast interpolating variable: {var}")
        print(f"Time taken to interpolate: {elapsed_time:.4f} seconds")
        print('===============================')
    
    return interpolated_var, (grid_x, grid_y)


def fast_interpolate_multiple_vars_to_grid(data, var_list, Ngrid_x=2048, Ngrid_y=2048,
                                         x_range=None, y_range=None, k_neighbors=DEFAULT_KD_NEIGHBOURS,
                                         leafsize=DEFAULT_KD_TREE_LEAFSIZE, verbose=False):
    """
    Highly optimized function to interpolate multiple variables simultaneously.
    This is much faster than calling fast_interpolate_var_to_grid multiple times.
    
    Parameters:
    - data (dict): Dictionary containing 'center_x', 'center_y', and variable data
    - var_list (list): List of variable names to interpolate
    - Ngrid_x (int): Number of grid points in x direction
    - Ngrid_y (int): Number of grid points in y direction
    - x_range (tuple): Optional x-axis bounds (xmin, xmax)
    - y_range (tuple): Optional y-axis bounds (ymin, ymax)
    - k_neighbors (int): Number of neighbors for KD-tree interpolation
    - leafsize (int): Leaf size for KD-tree construction
    
    Returns:
    - dict: Dictionary with keys as variable names and values as interpolated 2D arrays
    - tuple: (grid_x, grid_y) coordinate arrays
    """
    if verbose:
        print('===============================')
        print(f"Started fast interpolating {len(var_list)} variables: {var_list}")
    start_time = time.time()
        
    center_x, center_y = data["center_x"], data["center_y"]
    

    # Create mask for spatial filtering
    mask = np.ones(center_x.shape, dtype=bool)
    
    if x_range is not None:
        x_mask = (center_x >= x_range[0]) & (center_x <= x_range[1])
        mask &= x_mask
        
    if y_range is not None:
        y_mask = (center_y >= y_range[0]) & (center_y <= y_range[1])
        mask &= y_mask
    

    # Filter data
    filtered_center_x = center_x[mask]
    filtered_center_y = center_y[mask]
    
    # Create uniform grid
    grid_x_1d = np.linspace(filtered_center_x.min(), filtered_center_x.max(), Ngrid_x, dtype=np.float32)
    grid_y_1d = np.linspace(filtered_center_y.min(), filtered_center_y.max(), Ngrid_y, dtype=np.float32)
    grid_x, grid_y = np.meshgrid(grid_x_1d, grid_y_1d, indexing='ij')
    

    
    # Prepare source points for KD-tree
    source_points = np.column_stack((
        np.ascontiguousarray(filtered_center_x, dtype=np.float32),
        np.ascontiguousarray(filtered_center_y, dtype=np.float32)
    ))
    
    # Build KD-tree
    tree_start = time.time()
    
    tree = KDTree(source_points, leafsize=leafsize)
    

    
    # Query neighbors for all grid points
    grid_points = np.column_stack((grid_x.ravel(), grid_y.ravel()))
    

    
    distances, indices = tree.query(grid_points, k=k_neighbors)
    

    tree_time = time.time() - tree_start
    
    # Prepare all variable data in a single array for efficient processing
    prep_start = time.time()
    all_values = np.empty((len(var_list), len(filtered_center_x)), dtype=np.float32)
    for i, var in enumerate(var_list):
        all_values[i] = np.ascontiguousarray(var[mask], dtype=np.float32)
    

    distances = distances.astype(np.float32)
    indices = indices.astype(np.int32)
    prep_time = time.time() - prep_start    
    # Interpolate all variables at once using the optimized Numba function
    interp_start = time.time()
    interpolated_values = _fast_interpolate_multiple_vars(
        indices, distances, all_values, k_neighbors, len(var_list)
    )
    

    
    # Reshape results into dictionary
    # results = {}
    results = []
    for i, var in enumerate(var_list):
        # results[var] = interpolated_values[i].reshape(grid_x.shape)
        results.append(interpolated_values[i].reshape(grid_x.shape))
        

    
    interp_time = time.time() - interp_start
    total_time = time.time() - start_time
    

    
    if verbose:
        print(f"Interpolation completed in {interp_time:.4f} seconds")
        print(f"Total time for {len(var_list)} variables: {total_time:.4f} seconds")
        print('===============================')
    
    return results, (grid_x, grid_y)


# Import the fill_nan function from your existing library or define it here
def fill_nan(grid):
    """
    Fills NaN values in a 2D array using nearest neighbor interpolation.
    This is the same function from your existing library.
    """
    from scipy.interpolate import griddata
    
    nan_mask = np.isnan(grid)
    if np.any(nan_mask):
        # Find the indices of non-NaN values
        x_non_nan, y_non_nan = np.where(~nan_mask)
        non_nan_values = grid[~nan_mask]
        
        # Fill NaN values with nearest neighbor interpolation
        grid[nan_mask] = griddata(
            (x_non_nan, y_non_nan),
            non_nan_values,
            (np.where(nan_mask)[0], np.where(nan_mask)[1]),
            method='nearest'
        )
    return grid






# Spherical Integration 2D
sig_sph_int_2d_32 = types.float32[:](
    types.float32[:,:],
    types.float32[:,:],
    types.float32[:],
    types.int64
    )

sig_sph_int_2d_64 = types.float64[:](
    types.float64[:,:],
    types.float64[:,:],
    types.float64[:],
    types.int64
    )

@njit([sig_sph_int_2d_32, sig_sph_int_2d_64], parallel=True, fastmath=True, cache=True)
def spherical_integrate_2D_core(
    data: np.ndarray,
    r: np.ndarray,
    bin_edges: np.ndarray,
    bins: int) -> np.ndarray:
    """
    2D spherical integration using thread-local accumulation pattern.
    """
    data_flat = data.ravel()
    r_flat = r.ravel()
    n_elements = len(data_flat)
    n_edges = len(bin_edges)
    
    # Get number of threads
    n_threads = numba.config.NUMBA_NUM_THREADS
    
    # Create thread-local accumulators
    local_sums = np.zeros((n_threads, bins), dtype=data.dtype)
    
    # Parallel binning
    for i in prange(n_elements):
        r_val = r_flat[i]
        
        # Binary search for bin
        if r_val < bin_edges[0]:
            bin_idx = 0
        elif r_val >= bin_edges[-1]:
            bin_idx = n_edges
        else:
            left = 0
            right = n_edges - 1
            
            while left < right:
                mid = (left + right) // 2
                if r_val < bin_edges[mid]:
                    right = mid
                else:
                    left = mid + 1
            
            bin_idx = left
        
        # Accumulate
        if 1 <= bin_idx <= bins:
            thread_id = numba.get_thread_id()
            local_sums[thread_id, bin_idx - 1] += data_flat[i]
    
    # Merge results
    radial_sum = np.zeros(bins, dtype=data.dtype)
    for i in range(bins):
        for t in range(n_threads):
            radial_sum[i] += local_sums[t, i]
    
    return radial_sum

sig_rad_dist_2D = types.float32[:,:](
    types.UniTuple(types.int64, 2)
    )

@njit(sig_rad_dist_2D, parallel=True, fastmath=True, cache=True)
def compute_radial_distances_2D_core(
    shape : tuple) -> np.ndarray:
    """
    Parallel computation of radial distances from center
    Shape: (nx, ny)
    """
    nx, ny = shape
    center_x = (nx - 1) / 2.0
    center_y = (ny - 1) / 2.0
    
    r = np.empty(shape, dtype=np.float32)
    
    for x in prange(nx):
        for y in range(ny):
            dx = x - center_x
            dy = y - center_y
            r[x, y] = np.sqrt(dx*dx + dy*dy)
    return r




import numpy as np
from numpy.fft import fft2, rfft2, fftshift, fftfreq, rfftfreq
import time
from interp import * 


def ensure_float32(field, field_name="field"):
    """
    Ensures the input field is float32 for memory efficiency.
    """
    if field.dtype != np.float32:
        # print(f"Converting {field_name} from {field.dtype} to float32")
        field = field.astype(np.float32)
    return field


def compute_power_spectrum_2D_interpolated(
    interpolated_field: np.ndarray,
    grid_x: np.ndarray = None,
    grid_y: np.ndarray = None,
    field_name: str = "field",
    real_input: bool = True,
    norm: str = 'forward',
    verbose=False) -> np.ndarray:
    """
    Computes 2D power spectrum from uniformly interpolated data.

    Parameters:
    -----------
    interpolated_field : np.ndarray
        2D interpolated field data (output from interpolation functions)
    grid_x : np.ndarray, optional
        x-coordinate grid (for calculating physical scales). If None, assumes unit spacing.
    grid_y : np.ndarray, optional  
        y-coordinate grid (for calculating physical scales). If None, assumes unit spacing.
    field_name : str
        Name of the field for logging purposes
    real_input : bool
        Whether the input field is real-valued (affects FFT method choice)
    norm : str
        FFT normalization ('forward', 'backward', or 'ortho')
        
    Returns:
    --------
    np.ndarray
        2D power spectrum, fft-shifted so k=0 is at center
    """
    if verbose:
        print('===============================')
        print(f"Computing 2D power spectrum for field: {field_name}")
    start_time = time.time()
    
    # Ensure input is 2D
    if len(interpolated_field.shape) != 2:
        raise ValueError(f"Field should be 2D, got shape {interpolated_field.shape}")
    
    # Ensure data is float32 for memory efficiency
    field = ensure_float32(interpolated_field, field_name=field_name)
    
    # Get field dimensions
    Nx, Ny = field.shape
    
    # Calculate physical spacing if grids are provided
    if grid_x is not None and grid_y is not None:
        # Assume uniform spacing (which should be true for interpolated data)
        dx = (grid_x.max() - grid_x.min()) / (Nx - 1)
        dy = (grid_y.max() - grid_y.min()) / (Ny - 1)
    else:
        # Default to unit spacing
        dx = dy = 1.0
    
    # Compute FFT
    fft_start = time.time()
    if real_input and np.isrealobj(field):
        # Use real FFT for better efficiency with real-valued data
        field_fft = rfft2(field, norm=norm)
    else:
        # Use complex FFT
        field_fft = fft2(field, norm=norm)
    
    fft_time = time.time() - fft_start
    
    # Compute power spectrum
    power_spectrum = np.abs(field_fft)**2
    
    # Handle real FFT case - need to reconstruct full spectrum for proper shifting
    if real_input and np.isrealobj(field) and field_fft.shape[-1] != field.shape[-1]:
        # Create full power spectrum for real FFT
        full_power = np.zeros((Nx, Ny), dtype=power_spectrum.dtype)
        full_power[:, :power_spectrum.shape[-1]] = power_spectrum
        
        # Mirror the positive frequencies to negative frequencies
        if Ny % 2 == 0:
            # Even number of points
            full_power[:, -Ny//2+1:] = power_spectrum[:, 1:Ny//2][:, ::-1]
        else:
            # Odd number of points  
            full_power[:, -Ny//2:] = power_spectrum[:, 1:Ny//2+1][:, ::-1]
        
        power_spectrum = full_power
    
    # Shift so that k=0 is at the center
    power_spectrum_shifted = fftshift(power_spectrum)
    
    # Create wavenumber arrays
    # Note: fftfreq gives frequencies, multiply by 2π to get wavenumbers
    k_x_1d = fftshift(fftfreq(Nx, d=dx)) * 2 * np.pi
    k_y_1d = fftshift(fftfreq(Ny, d=dy)) * 2 * np.pi
    
    # Create 2D wavenumber grids
    k_x, k_y = np.meshgrid(k_x_1d, k_y_1d, indexing='ij')
    
    end_time = time.time()
    total_time = end_time - start_time
    
    if verbose:
        print(f"Power spectrum shape: {power_spectrum_shifted.shape}")
        print(f"Total computation time: {total_time:.4f} seconds")
        print('===============================')
        
    return power_spectrum_shifted



DEFAULT_BINS_RATIO = 2  # N // 2
DEFAULT_BIN_MIN = 0.5  # the smallest k mode

def spherical_integrate_2D(
    field: np.ndarray, 
    bins: int = None,
    field_name : str = "field") -> tuple:
    """
    2D spherical integration using JIT-compiled functions.
    """
    N = field.shape[0]
    if not bins:
        bins = N // DEFAULT_BINS_RATIO
        
    # Ensure data is float32 for memory efficiency
    field = ensure_float32(
        field, 
        field_name=field_name)
    
    # Use JIT function for distances
    r = compute_radial_distances_2D_core(
        field.shape)
    bin_edges = np.linspace(DEFAULT_BIN_MIN, bins, bins + 1).astype(field.dtype)
    radial_sum = spherical_integrate_2D_core(
        field, 
        r, 
        bin_edges, 
        bins)    
    k_modes = np.ceil((bin_edges[:-1] + bin_edges[1:]) / 2)

    return k_modes, radial_sum

