

# filename = '/Users/michaelgrehan/Downloads/data0001.hdf5'
# variable_names = ['rho', 'u1', 'u2', 'u3', 'p', 'b1', 'b2', 'b3', 'e1', 'e2', 'e3', 's', 'lfac', 'xi',
#                   'jz', 'v1', 'v2', 'v3', 'divB', 'curlbz', 'curlvz', 'divE']
# data = read_HDF5(filename=filename, variable_names=variable_names)
# plot_exact_cells_from_data(data, lambda vars: vars['b1']**2+vars['b2']**2,
#                            x_range=(-0.5,0.5), y_range=(-0.5,0.5),
#                            vmin=0, vmax=1, cmap=cmr.ember, pad=0.0, label='$B_\\perp^2/B^2_0$')


# plot_cells(data, lambda vars: np.sqrt(vars['b1']**2+vars['b2']**2) ,
#                            fig=fig, ax=ax,
#                            x_range=(-0.5,0.0), y_range=(-0.5,0.0),
#                            vmin=0, vmax=1, cmap=cmr.ember, pad=0.0, label='$B_\\perp/B_0$')

# plot_blocks(data, 
#             fig=fig, ax=ax,
#             x_range=(-0.5,0.0), y_range=(-0.5,0.0),
#             color='w', linewidth=0.25)



import numpy as np
import h5py
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import cmasher as cmr
import time

from scipy.integrate import cumulative_trapezoid as cumtrapz
from scipy.interpolate import NearestNDInterpolator
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

from interp import *

from matplotlib.patches import Wedge
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
from matplotlib.ticker import FuncFormatter

from scipy.optimize import root_scalar



def read_HDF5(filename, variable_names=None, blocks=None, coord_system='cartesian'):
    """
    Extract cell-centered simulation data from HDF5 file.
    Supports both 1D and 2D data, auto-detected from file contents.

    Parameters:
    -----------
    filename : str
        Path to HDF5 file
    variable_names : list, optional
        Names for variables. If None, uses default names.
    blocks : array-like, optional
        Specific block indices to read. If None, reads all blocks.
    coord_system : str, optional
        'cartesian' or 'polar'. Determines how cell centers are calculated.
        Only relevant for 2D data.

    Returns:
    --------
    dict : Data dictionary with flattened arrays and metadata.
             Always contains 'center_x'. For 2D data also contains 'center_y'.
             Also contains 'ndim' key indicating dimensionality.
    """

    with h5py.File(filename, 'r') as f:
        cc_vars = f['Cell-centered Variables'][:]   # (nblocks, nvars, nx) or (nblocks, nvars, nx, ny)
        grid    = f['Grid'][:]                       # (nblocks, ngrid_pts) or (nblocks, ngrid_pts, 2)
        levels  = f['Levels'][:]                     # (nblocks,)

        GZ        = int(f.attrs.get('GhostZones', 2))
        time_val  = f.attrs.get('Time', 0.0)
        iteration = f.attrs.get('Iteration', 0)

    # ------------------------------------------------------------------ #
    # Auto-detect dimensionality from cc_vars shape
    # ------------------------------------------------------------------ #
    if cc_vars.ndim == 3:
        ndim = 1
        nblocks, nvars, block_size_x = cc_vars.shape
        block_size_y = 1
    elif cc_vars.ndim == 4:
        ndim = 2
        nblocks, nvars, block_size_x, block_size_y = cc_vars.shape
    else:
        raise ValueError(f"Unexpected cc_vars shape: {cc_vars.shape}")

    print(f"Detected {ndim}D data")

    # ------------------------------------------------------------------ #
    # Work out grid dimensions
    # ------------------------------------------------------------------ #
    if ndim == 1:
        # grid shape: (nblocks, ngrid_pts) or (nblocks, ngrid_pts, 1)
        if grid.ndim == 3:
            grid = grid[:, :, 0]          # drop trailing size-1 dimension
        ngrid_pts = grid.shape[1]         # number of grid points per block
        grid_nx   = ngrid_pts             # 1-D: points along x
        grid_ny   = 1

        active_nx = block_size_x - 2 * GZ
        active_ny = 1

    else:  # 2D
        ngrid = grid.shape[1]
        s = int(np.sqrt(ngrid))
        if s * s == ngrid:
            grid_nx, grid_ny = s, s
        else:
            factors = [(i, ngrid // i) for i in range(1, int(np.sqrt(ngrid)) + 1)
                       if ngrid % i == 0]
            grid_nx, grid_ny = min(factors, key=lambda x: abs(x[0] - x[1]))

        active_nx = block_size_x - 2 * GZ
        active_ny = block_size_y - 2 * GZ

    print(f"Found {nblocks} blocks, {nvars} variables")
    print(f"Ghost zones: {GZ}, Active cells per block: {active_nx}" +
          (f"x{active_ny}" if ndim == 2 else ""))
    print(f"Coordinate system: {coord_system if ndim == 2 else 'N/A (1D)'}")

    # ------------------------------------------------------------------ #
    # Block selection and variable names
    # ------------------------------------------------------------------ #
    block_indices = np.arange(nblocks) if blocks is None else np.array(blocks, dtype=int)

    if variable_names is None:
        known_names = ['rho', 'u1', 'u2', 'u3', 'p',
                       'b1', 'b2', 'b3', 'e1', 'e2', 'e3', 's', 'lfac', 'xi']
        variable_names = [
            known_names[i] if i < len(known_names) else f'var{i}'
            for i in range(nvars)
        ]
    elif len(variable_names) != nvars:
        raise ValueError(f"Number of variable names ({len(variable_names)}) "
                         f"must match number of variables in data ({nvars})")

    # ------------------------------------------------------------------ #
    # Allocate output arrays
    # ------------------------------------------------------------------ #
    cells_per_block = active_nx * active_ny
    total_cells     = len(block_indices) * cells_per_block

    center_x    = np.zeros(total_cells)
    center_y    = np.zeros(total_cells) if ndim == 2 else None
    data_arrays = {name: np.zeros(total_cells) for name in variable_names}

    xmin_global = np.inf;  xmax_global = -np.inf
    ymin_global = np.inf;  ymax_global = -np.inf

    cell_idx = 0

    # ------------------------------------------------------------------ #
    # Main loop over blocks
    # ------------------------------------------------------------------ #
    for b in block_indices:

        if ndim == 1:
            # grid[b] shape: (ngrid_pts,) — coordinates of cell edges
            coord1 = grid[b]                              # shape (ngrid_pts,)

            # Cell centres: midpoints of adjacent edge pairs.
            # The 1D grid only stores edges of *active* cells (no ghost edges),
            # so this already gives exactly active_nx centres — no ghost stripping needed.
            center_x_block = 0.5 * (coord1[:-1] + coord1[1:])  # shape (active_nx,)

            xmin_global = min(xmin_global, coord1.min())
            xmax_global = max(xmax_global, coord1.max())
            
            # Unpack variable data — shape (nvars, block_size_x)
            for var_idx, var_name in enumerate(variable_names):
                var_data = cc_vars[b, var_idx, GZ: GZ + active_nx]
                data_arrays[var_name][cell_idx: cell_idx + active_nx] = var_data

            center_x[cell_idx: cell_idx + active_nx] = center_x_block

        if ndim == 2:  # 2D (original logic, unchanged)
            g      = grid[b].reshape(grid_nx, grid_ny, 2)
            coord1 = g[:, :, 0]
            coord2 = g[:, :, 1]

            if coord_system.lower() == 'cartesian':
                center_x_block = 0.5 * (coord1[:-1, :-1] + coord1[1:, 1:])
                center_y_block = 0.5 * (coord2[:-1, :-1] + coord2[1:, 1:])

            elif coord_system.lower() == 'polar':
                x1, x2 = coord1[:-1, :-1], coord1[1:, :-1]
                x3, x4 = coord1[:-1,  1:], coord1[1:,  1:]
                y1, y2 = coord2[:-1, :-1], coord2[1:, :-1]
                y3, y4 = coord2[:-1,  1:], coord2[1:,  1:]

                r1 = np.sqrt(x1**2 + y1**2);  r2 = np.sqrt(x2**2 + y2**2)
                r3 = np.sqrt(x3**2 + y3**2);  r4 = np.sqrt(x4**2 + y4**2)

                theta1 = np.arctan2(y1, x1);  theta2 = np.arctan2(y2, x2)
                theta3 = np.arctan2(y3, x3);  theta4 = np.arctan2(y4, x4)

                rmin = np.minimum(np.minimum(r1, r2), np.minimum(r3, r4))
                rmax = np.maximum(np.maximum(r1, r2), np.maximum(r3, r4))
                theta_min = np.minimum(np.minimum(theta1, theta2), np.minimum(theta3, theta4))
                theta_max = np.maximum(np.maximum(theta1, theta2), np.maximum(theta3, theta4))

                r_center     = 0.5 * (rmin + rmax)
                theta_center = 0.5 * (theta_min + theta_max)

                center_x_block = r_center * np.cos(theta_center)
                center_y_block = r_center * np.sin(theta_center)

            else:
                raise ValueError(f"Unknown coordinate system: {coord_system}")

            xmin_global = min(xmin_global, coord1.min())
            xmax_global = max(xmax_global, coord1.max())
            ymin_global = min(ymin_global, coord2.min())
            ymax_global = max(ymax_global, coord2.max())

            for var_idx, var_name in enumerate(variable_names):
                var_data = cc_vars[b, var_idx, GZ:-GZ, GZ:-GZ]
                data_arrays[var_name][cell_idx: cell_idx + cells_per_block] = var_data.ravel()

            center_x[cell_idx: cell_idx + cells_per_block] = center_x_block.ravel()
            center_y[cell_idx: cell_idx + cells_per_block] = center_y_block.ravel()

        cell_idx += cells_per_block

    # ------------------------------------------------------------------ #
    # Build result dict
    # ------------------------------------------------------------------ #
    result = {
        'time': time_val,
        'ndim': ndim,
        'metadata': {
            'iteration':              iteration,
            'nblocks_processed':      len(block_indices),
            'total_cells':            total_cells,
            'active_cells_per_block': cells_per_block,
            'ghost_zones':            GZ,
            'levels':                 levels if blocks is None else levels[block_indices],
            'variable_names':         variable_names,
            'filename':               filename,
            'grid_dimensions':        (grid_nx, grid_ny) if ndim == 2 else (grid_nx,),
            'coord_system':           coord_system if ndim == 2 else 'cartesian',
            'xmin': xmin_global,
            'xmax': xmax_global,
            **({'ymin': ymin_global, 'ymax': ymax_global} if ndim == 2 else {}),
        }
    }

    result['center_x'] = center_x
    if ndim == 2:
        result['center_y'] = center_y

    result['coord_system'] = coord_system if ndim == 2 else 'cartesian'
    result.update(data_arrays)

    return result



from matplotlib.collections import PatchCollection


def get_block_structure(data, variables=None):
    """
    Reconstruct block structure from flattened data.
    Supports both cartesian and polar coordinate systems.
    """
    
    coord_system = data['coord_system']
    
    if 'block_data' in data:
        print("Warning: Data already has block structure")
        return data

    filename = data['metadata']['filename']
    nblocks_processed = int(data['metadata']['nblocks_processed'])
    ghost_zones = data['metadata']['ghost_zones']
    variable_names = data['metadata']['variable_names']

    if variables is None:
        variables = variable_names

    # Read only grid info
    with h5py.File(filename, 'r') as f:
        grid = f['Grid'][:]

    # Determine grid dimensions (handle non-square grids)
    total_blocks_on_file = grid.shape[0]
    if nblocks_processed > total_blocks_on_file:
        raise ValueError(f"metadata nblocks_processed ({nblocks_processed}) "
                         f"> grid blocks in file ({total_blocks_on_file})")

    ngrid = grid.shape[1]
    s = int(np.sqrt(ngrid))
    if s * s == ngrid:
        grid_nx, grid_ny = s, s
    else:
        factors = [(i, ngrid // i) for i in range(1, int(np.sqrt(ngrid)) + 1) if ngrid % i == 0]
        grid_nx, grid_ny = min(factors, key=lambda x: abs(x[0] - x[1]))

    active_size_x = int(data['metadata']['grid_dimensions'][0]) - 1
    active_size_y = int(data['metadata']['grid_dimensions'][1]) - 1
    cells_per_block = active_size_x * active_size_y

    # Reshape grid into (nblocks, grid_nx, grid_ny, 2)
    grid_nodes = grid.reshape(grid.shape[0], grid_nx, grid_ny, 2)
    grid_nodes = grid_nodes[:nblocks_processed]
    n = nblocks_processed

    # Compute block boundaries
    if coord_system.lower() == 'cartesian':
        # For cartesian: corners are simply the four corners of the grid
        block_boundaries = np.stack([
            grid_nodes[:, 0, 0, :],     # Bottom-left
            grid_nodes[:, -1, 0, :],    # Bottom-right
            grid_nodes[:, -1, -1, :],   # Top-right
            grid_nodes[:, 0, -1, :]     # Top-left
        ], axis=1)
        
    elif coord_system.lower() == 'polar':
        # For polar: need to compute wedge bounds from all corner points
        # Extract all four corners
        bl = grid_nodes[:, 0, 0, :]      # (n, 2)
        br = grid_nodes[:, -1, 0, :]     # (n, 2)
        tr = grid_nodes[:, -1, -1, :]    # (n, 2)
        tl = grid_nodes[:, 0, -1, :]     # (n, 2)
        
        # Convert all corners to polar coordinates
        r_bl = np.sqrt(bl[:, 0]**2 + bl[:, 1]**2)
        r_br = np.sqrt(br[:, 0]**2 + br[:, 1]**2)
        r_tr = np.sqrt(tr[:, 0]**2 + tr[:, 1]**2)
        r_tl = np.sqrt(tl[:, 0]**2 + tl[:, 1]**2)
        
        theta_bl = np.arctan2(bl[:, 1], bl[:, 0])
        theta_br = np.arctan2(br[:, 1], br[:, 0])
        theta_tr = np.arctan2(tr[:, 1], tr[:, 0])
        theta_tl = np.arctan2(tl[:, 1], tl[:, 0])
        
        # Stack all angles: shape (n, 4)
        all_angles = np.stack([theta_bl, theta_br, theta_tr, theta_tl], axis=1)
        # Unwrap along axis=1 (across the 4 corners for each block)
        all_angles_unwrapped = np.unwrap(all_angles, axis=1)
        # Extract back
        theta_bl = all_angles_unwrapped[:, 0]
        theta_br = all_angles_unwrapped[:, 1]
        theta_tr = all_angles_unwrapped[:, 2]
        theta_tl = all_angles_unwrapped[:, 3]
        
        # Find wedge bounds
        rmin = np.minimum(np.minimum(r_bl, r_br), np.minimum(r_tr, r_tl))
        rmax = np.maximum(np.maximum(r_bl, r_br), np.maximum(r_tr, r_tl))
        theta_min = np.minimum(np.minimum(theta_bl, theta_br), np.minimum(theta_tr, theta_tl))
        theta_max = np.maximum(np.maximum(theta_bl, theta_br), np.maximum(theta_tr, theta_tl))
        
            
        # Store wedge parameters instead of corner points
        # Format: (rmin, rmax, theta_min, theta_max) for each block
        block_boundaries = np.stack([rmin, rmax, theta_min, theta_max], axis=1)  # (n, 4)
        
    else:
        raise ValueError(f"Unknown coordinate system: {coord_system}")

    # Prepare block grid nodes as a list of arrays
    block_grid_nodes = [grid_nodes[i].copy() for i in range(n)]

    # Vectorized unflatten of variables
    available_vars = [v for v in variables if v in data]
    reshaped_vars = {}
    for var in available_vars:
        flat = np.asarray(data[var])
        expected_len = n * cells_per_block
        if flat.size != expected_len:
            raise ValueError(f"Variable '{var}' has length {flat.size}, expected {expected_len}")
        reshaped = flat.reshape(n, active_size_x, active_size_y)
        reshaped_vars[var] = reshaped

    block_cell_data = [
        {var: reshaped_vars[var][i].copy() for var in reshaped_vars}
        for i in range(n)
    ]

    result = data.copy()
    result['block_data'] = {
        'nblocks': n,
        'active_size_x': active_size_x,
        'active_size_y': active_size_y,
        'ghost_zones': ghost_zones,
        'coord_system': coord_system,
        'block_indices': list(range(n)),
        'block_grid_nodes': block_grid_nodes,
        'block_cell_data': block_cell_data,
        'block_boundaries': block_boundaries,
        'cells_per_block': cells_per_block
    }

    return result




def plot_blocks(data, fig=None, ax=None, linewidth=0.1, color='k', 
                x_range=None, y_range=None, title=None, figsize=(6.4, 4.8)):
    """
    Plot block boundaries using flattened data.
    Automatically reconstructs block structure and handles both cartesian and polar.
    """
    
    coord_system = data['coord_system']

    data = get_block_structure(data, variables=[])
    
    print('===============================')
    print(f"Started plotting block boundaries ({coord_system})")
    
    block_boundaries = data['block_data']['block_boundaries']
    N = block_boundaries.shape[0]
    coord_system = data['block_data']['coord_system']
    
    segments = []
    patches = []  # <- CHANGE: Store patches in a list first
    
    if coord_system.lower() == 'cartesian':
        # Original cartesian logic
        for i in range(N):
            block_corners = block_boundaries[i]  # Shape: (4, 2)
            
            x_coords = block_corners[:, 0]
            y_coords = block_corners[:, 1]
            
            if x_range is not None:
                if np.all(x_coords < x_range[0]) or np.all(x_coords > x_range[1]):
                    continue
            if y_range is not None:
                if np.all(y_coords < y_range[0]) or np.all(y_coords > y_range[1]):
                    continue
            
            segments.extend([
                [block_corners[0], block_corners[1]],  # Bottom edge
                [block_corners[1], block_corners[2]],  # Right edge
                [block_corners[2], block_corners[3]],  # Top edge
                [block_corners[3], block_corners[0]]   # Left edge
            ])
            
    elif coord_system.lower() == 'polar':
        # Polar wedge logic using Wedge patches
        wedges = []
        
        for i in range(N):
            rmin, rmax, theta_min, theta_max = block_boundaries[i]
            
            # Check spatial range in Cartesian
            if x_range is not None or y_range is not None:
                # Compute bounding box of the wedge
                corners_x = [
                    rmin * np.cos(theta_min), rmin * np.cos(theta_max),
                    rmax * np.cos(theta_min), rmax * np.cos(theta_max)
                ]
                corners_y = [
                    rmin * np.sin(theta_min), rmin * np.sin(theta_max),
                    rmax * np.sin(theta_min), rmax * np.sin(theta_max)
                ]
                xmin, xmax = min(corners_x), max(corners_x)
                ymin, ymax = min(corners_y), max(corners_y)
                
                if x_range is not None:
                    if xmax < x_range[0] or xmin > x_range[1]:
                        continue
                if y_range is not None:
                    if ymax < y_range[0] or ymin > y_range[1]:
                        continue
            
            # Convert theta from radians to degrees for Wedge
            theta1_deg = np.degrees(theta_min)
            theta2_deg = np.degrees(theta_max)
            
            # Create a Wedge patch centered at origin
            wedge = Wedge(center=(0, 0), r=rmax, theta1=theta1_deg, theta2=theta2_deg,
                         width=rmax - rmin, facecolor='none', edgecolor=color, linewidth=linewidth)
            wedges.append(wedge)
        
        patches = wedges  # <- CHANGE: Just store, don't add collection yet
    
    # FIX: Create figure and axis BEFORE adding collections
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    
    # FIX: Add collections AFTER axis exists
    if coord_system.lower() == 'cartesian' and segments:
        lc = LineCollection(segments, linewidths=linewidth, colors=color)
        ax.add_collection(lc)
    elif patches:  # <- CHANGE: Check if we have patches to add
        patch_collection = PatchCollection(patches, facecolors='none', 
                                          edgecolors=color, linewidths=linewidth,
                                          match_original=True)  # <- CHANGE: Added match_original
        ax.add_collection(patch_collection)
    
    ax.set_aspect('equal')
    
    if coord_system.lower() == 'cartesian':
        ax.set_xlabel('$x/L$')
        ax.set_ylabel('$y/L$')
        if title is None:
            t = data['time']
            title = f'$t/t_c = {t:.1f}$'
    elif coord_system.lower() == 'polar':
        ax.set_xlabel('$x/R_{\\rm NS}$')
        ax.set_ylabel('$z/R_{\\rm NS}$')
        if title is None:
            t = data['time']
            title = f'$t = {t:.1f}\\,$'+'$[R_{\\rm NS}/c]$'
    
    # Set plot limits
    if segments:
        all_points = np.array([pt for seg in segments for pt in seg])
        if x_range is None:
            ax.set_xlim(all_points[:, 0].min(), all_points[:, 0].max())
        else:
            ax.set_xlim(x_range)
        
        if y_range is None:
            ax.set_ylim(all_points[:, 1].min(), all_points[:, 1].max())
        else:
            ax.set_ylim(y_range)
    elif patches:  # <- CHANGE: Set limits for polar case
        # Set limits based on block boundaries
        if coord_system.lower() == 'polar':
            rmax_all = block_boundaries[:, 1].max()
            if x_range is None:
                ax.set_xlim(-rmax_all * 1.1, rmax_all * 1.1)
            else:
                ax.set_xlim(x_range)
            if y_range is None:
                ax.set_ylim(-rmax_all * 1.1, rmax_all * 1.1)
            else:
                ax.set_ylim(y_range)
    
    ax.set_title(title)
    
    print("Finished plotting block boundaries")
    print('===============================')
    
    return fig, ax
 
 
def plot_cells(data, variable, 
               fig=None, ax=None,
               figsize=(6.4, 4.8), vmin=None, vmax=None, 
               cmap=None, x_range=None, y_range=None, label='',
               title=None, save_path=None, dpi=300, colorbar=True,
               orientation='vertical', location='right', pad=0.1,
               edgecolors=None, linewidths=0.1, use_log_norm=False,
               aspect=20, shrink=0.95):
    """
    Plot exact simulation cells using preserved block structure from extracted data.
    Supports both cartesian and polar coordinate systems.
    
    Accepts:
      - variable: str (variable name),
      - variable: numpy array (1D flat of length total_cells),
      - variable: callable(data) -> flat array or (nblocks, nx, ny) array,
      - variable: callable(block_vars) -> (nx, ny) array for that block.
    """
    coord_system = data['coord_system']
    
    print('===============================')
    print(f"Started plotting ({coord_system}).")
    
    
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    
    # Ensure block structure available
    if 'block_data' not in data:
        if callable(variable):
            needed_vars = data['metadata']['variable_names']
        elif isinstance(variable, str):
            needed_vars = [variable]
        else:
            needed_vars = data['metadata']['variable_names']
        data = get_block_structure(data, needed_vars)

    block_info = data['block_data']
    nblocks = block_info['nblocks']
    active_size_x = block_info['active_size_x']
    active_size_y = block_info['active_size_y']
    cells_per_block = block_info['cells_per_block']
    total_cells = nblocks * cells_per_block

    # Normalized storage for per-block arrays (fast path)
    per_block_arrays = None
    use_global_precompute = False

    # ---------- handle numpy array / list input ----------
    if isinstance(variable, (np.ndarray, list, tuple)):
        arr = np.asarray(variable)
        if arr.ndim == 1:
            if arr.size == total_cells:
                # flatten for entire dataset -> reshape to (nblocks, nx, ny)
                per_block_arrays = arr.reshape(nblocks, active_size_x, active_size_y).copy()
            elif arr.size == cells_per_block:
                # single block provided -> broadcast to all blocks
                per_block_arrays = np.tile(arr.reshape(1, active_size_x, active_size_y), (nblocks, 1, 1)).copy()
            else:
                raise ValueError(f"Array input has length {arr.size}; expected {total_cells} (full) or {cells_per_block} (per-block).")
        elif arr.ndim == 3 and arr.shape[0] == nblocks and arr.shape[1] == active_size_x and arr.shape[2] == active_size_y:
            per_block_arrays = arr.copy()
        else:
            raise ValueError(f"Array input has unsupported shape {arr.shape}.")
        if np.isnan(per_block_arrays).any():
            print("Warning: input array contains NaNs.")
        if np.isinf(per_block_arrays).any():
            print("Warning: input array contains Infs.")
    
    # ---------- handle string variable ----------
    elif isinstance(variable, str):
        per_block_arrays = None

    # ---------- handle callable ----------
    elif callable(variable):
        try:
            out = variable(data)
            out = np.asarray(out)
            if out.ndim == 1 and out.size == total_cells:
                per_block_arrays = out.reshape(nblocks, active_size_x, active_size_y).copy()
                use_global_precompute = True
            elif out.ndim == 3 and out.shape == (nblocks, active_size_x, active_size_y):
                per_block_arrays = out.copy()
                use_global_precompute = True
            elif out.ndim == 2 and out.shape == (active_size_x, active_size_y):
                per_block_arrays = np.tile(out[None, ...], (nblocks, 1, 1)).copy()
                use_global_precompute = True
            else:
                per_block_arrays = None
                use_global_precompute = False
                print(f"Note: callable(data) returned shape {out.shape}; will try per-block calls instead.")
        except TypeError:
            per_block_arrays = None
            use_global_precompute = False
        except Exception as e:
            print("Warning: callable(data) raised an exception; will try per-block calls instead.")
            print("  Exception:", repr(e))
            per_block_arrays = None
            use_global_precompute = False
    else:
        raise ValueError("Unsupported 'variable' type. Must be str, ndarray, or callable.")

    # ---------- build cells + values ----------
    all_cells = []
    all_values = []

    for block_idx in range(nblocks):
        G = block_info['block_grid_nodes'][block_idx]   # (grid_nx, grid_ny, 2)
        block_vars = block_info['block_cell_data'][block_idx]

        # Extract var_data for this block
        if per_block_arrays is not None:
            var_data = per_block_arrays[block_idx]
        else:
            if isinstance(variable, str):
                if variable not in block_vars:
                    raise ValueError(f"Variable '{variable}' not found in block data; available: {list(block_vars.keys())}")
                var_data = block_vars[variable]
            else:
                try:
                    var_data = variable(block_vars)
                except Exception as e:
                    raise RuntimeError(f"Callable failed for block {block_idx}: {e}")
                var_data = np.asarray(var_data)
                if var_data.shape != (active_size_x, active_size_y):
                    raise ValueError(f"Callable(block_vars) must return shape ({active_size_x},{active_size_y}), got {var_data.shape}")

        if coord_system.lower() == 'cartesian':
            # Build rectangular polygons
            Xn = G[:, :, 0]
            Yn = G[:, :, 1]
            
            x0 = Xn[:-1, :-1].ravel()
            y0 = Yn[:-1, :-1].ravel()
            x1 = Xn[1:, :-1].ravel()
            y1 = Yn[1:, :-1].ravel()
            x2 = Xn[1:, 1:].ravel()
            y2 = Yn[1:, 1:].ravel()
            x3 = Xn[:-1, 1:].ravel()
            y3 = Yn[:-1, 1:].ravel()
            
            block_cells = np.stack([
                np.stack([x0, y0], axis=1),
                np.stack([x1, y1], axis=1),
                np.stack([x2, y2], axis=1),
                np.stack([x3, y3], axis=1)
            ], axis=1)
            
        elif coord_system.lower() == 'polar':
            # Build wedge parameters (rmin, rmax, theta_min, theta_max) for each cell
            Xn = G[:, :, 0]
            Yn = G[:, :, 1]
            
            # Extract four corners of each cell
            x0 = Xn[:-1, :-1]
            y0 = Yn[:-1, :-1]
            x1 = Xn[1:, :-1]
            y1 = Yn[1:, :-1]
            x2 = Xn[1:, 1:]
            y2 = Yn[1:, 1:]
            x3 = Xn[:-1, 1:]
            y3 = Yn[:-1, 1:]
            
            # Convert all corners to polar
            r0 = np.sqrt(x0**2 + y0**2)
            r1 = np.sqrt(x1**2 + y1**2)
            r2 = np.sqrt(x2**2 + y2**2)
            r3 = np.sqrt(x3**2 + y3**2)
            
            theta0 = np.arctan2(y0, x0)
            theta1 = np.arctan2(y1, x1)
            theta2 = np.arctan2(y2, x2)
            theta3 = np.arctan2(y3, x3)
            
            # Find radial bounds (straightforward)
            rmin = np.minimum(np.minimum(r0, r1), np.minimum(r2, r3))
            rmax = np.maximum(np.maximum(r0, r1), np.maximum(r2, r3))
            
            # FIX: Handle angle wrapping for theta bounds
            # Unwrap angles relative to theta0 to handle 2π wraparound
            theta1_unwrap = np.unwrap([theta0.ravel(), theta1.ravel()], axis=0)[1].reshape(theta0.shape)
            theta2_unwrap = np.unwrap([theta0.ravel(), theta2.ravel()], axis=0)[1].reshape(theta0.shape)
            theta3_unwrap = np.unwrap([theta0.ravel(), theta3.ravel()], axis=0)[1].reshape(theta0.shape)
            
            theta_min = np.minimum(np.minimum(theta0, theta1_unwrap), np.minimum(theta2_unwrap, theta3_unwrap))
            theta_max = np.maximum(np.maximum(theta0, theta1_unwrap), np.maximum(theta2_unwrap, theta3_unwrap))
            
            # Convert to degrees
            theta_min_deg = np.degrees(theta_min)
            theta_max_deg = np.degrees(theta_max)
            
            # Store as (rmin, rmax, theta_min_deg, theta_max_deg)
            block_cells = np.stack([rmin, rmax, theta_min_deg, theta_max_deg], axis=-1).reshape(-1, 4)
        
        else:
            raise ValueError(f"Unknown coordinate system: {coord_system}")
        
        block_values = var_data.ravel()
        
        # Mask by x_range / y_range
        mask = np.ones(len(block_cells), dtype=bool)
        
        if coord_system.lower() == 'cartesian':
            if x_range is not None:
                xmin, xmax = x_range
                mask &= (block_cells[:, :, 0].min(axis=1) <= xmax) & (block_cells[:, :, 0].max(axis=1) >= xmin)
            if y_range is not None:
                ymin, ymax = y_range
                mask &= (block_cells[:, :, 1].min(axis=1) <= ymax) & (block_cells[:, :, 1].max(axis=1) >= ymin)
        
        elif coord_system.lower() == 'polar':
            if x_range is not None or y_range is not None:
                rmin_cells = block_cells[:, 0]
                rmax_cells = block_cells[:, 1]
                theta_min_cells = np.radians(block_cells[:, 2])
                theta_max_cells = np.radians(block_cells[:, 3])
                
                corners_x = np.stack([
                    rmin_cells * np.cos(theta_min_cells),
                    rmin_cells * np.cos(theta_max_cells),
                    rmax_cells * np.cos(theta_min_cells),
                    rmax_cells * np.cos(theta_max_cells)
                ], axis=1)
                corners_y = np.stack([
                    rmin_cells * np.sin(theta_min_cells),
                    rmin_cells * np.sin(theta_max_cells),
                    rmax_cells * np.sin(theta_min_cells),
                    rmax_cells * np.sin(theta_max_cells)
                ], axis=1)
                
                if x_range is not None:
                    xmin, xmax = x_range
                    mask &= (corners_x.min(axis=1) <= xmax) & (corners_x.max(axis=1) >= xmin)
                if y_range is not None:
                    ymin, ymax = y_range
                    mask &= (corners_y.min(axis=1) <= ymax) & (corners_y.max(axis=1) >= ymin)
        
        all_cells.append(block_cells[mask])
        all_values.append(block_values[mask])

    if len(all_cells) == 0:
        raise ValueError("No cells found in specified range")
    
    all_cells = np.concatenate(all_cells, axis=0)
    all_values = np.concatenate(all_values, axis=0)
    
    if vmin is None:
        vmin = np.nanmin(all_values)
    if vmax is None:
        vmax = np.nanmax(all_values)
    if cmap is None:
        cmap = 'viridis'
    
    
    if use_log_norm:
        if vmin <= 0:
            print(f"Warning: LogNorm requested but vmin={vmin} <= 0. Using linear scale.")
            norm = None
        else:
            norm = LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = None
    
    # Create appropriate collection based on coordinate system
    if coord_system.lower() == 'cartesian':
        # Use PolyCollection for rectangles
        if norm is not None:
            poly_collection = PolyCollection(all_cells, norm=norm, array=all_values,
                                             cmap=cmap, edgecolors=edgecolors,
                                             linewidths=linewidths)
        else:
            poly_collection = PolyCollection(all_cells, array=all_values,
                                             cmap=cmap, edgecolors=edgecolors,
                                             clim=(vmin, vmax), linewidths=linewidths)
        ax.add_collection(poly_collection)
        collection = poly_collection
        
    
    elif coord_system.lower() == 'polar':
        if edgecolors == None:
            antialiased = False
        else:
            antialiased = True
        
        
        # Extract all parameters at once (vectorized)
        rmin = all_cells[:, 0]
        rmax = all_cells[:, 1]
        theta_min_deg = all_cells[:, 2]
        theta_max_deg = all_cells[:, 3]
        
        wedges = [Wedge(center=(0, 0), r=rmax[i], theta1=theta_min_deg[i], 
                    theta2=theta_max_deg[i], width=(rmax[i] - rmin[i]))
                for i in range(len(all_cells))]
        
        if norm is not None:
            collection = PatchCollection(wedges, norm=norm, array=all_values,
                                        cmap=cmap, edgecolors=edgecolors,
                                        linewidths=linewidths, match_original=True,
                                        antialiased=antialiased)
        else:
            collection = PatchCollection(wedges, array=all_values,
                                        cmap=cmap, edgecolors=edgecolors,
                                        linewidths=linewidths, clim=(vmin, vmax),
                                        match_original=True, antialiased=antialiased)
        ax.add_collection(collection)
    
    if colorbar:
        extend_type = 'neither'
        if np.any(all_values < vmin) and np.any(all_values > vmax):
            extend_type = 'both'
        elif np.any(all_values < vmin):
            extend_type = 'min'
        elif np.any(all_values > vmax):
            extend_type = 'max'
        
        plt.colorbar(collection, ax=ax, extend=extend_type,
                     label=label, orientation=orientation, location=location,
                     pad=pad, aspect=aspect, shrink=shrink)
    
    # Set axis limits and labels
    if coord_system.lower() == 'cartesian':
        all_x = all_cells[:, :, 0].ravel()
        all_y = all_cells[:, :, 1].ravel()
        ax.set_xlim(all_x.min(), all_x.max())
        ax.set_ylim(all_y.min(), all_y.max())
    elif coord_system.lower() == 'polar':
        if x_range is not None:
            ax.set_xlim(x_range)
        else:
            rmax_all = all_cells[:, 1].max()
            ax.set_xlim(-rmax_all * 1.1, rmax_all * 1.1)
        if y_range is not None:
            ax.set_ylim(y_range)
        else:
            rmax_all = all_cells[:, 1].max()
            ax.set_ylim(-rmax_all * 1.1, rmax_all * 1.1)
    
    ax.set_aspect('equal')
    
    if coord_system.lower() == 'cartesian':
        ax.set_xlabel('$x/L$')
        ax.set_ylabel('$y/L$')
        if title is None:
            t = data['time']
            title = f'$t/t_c = {t:.1f}$'
    elif coord_system.lower() == 'polar':
        ax.set_xlabel('$x/R_{\\rm NS}$')
        ax.set_ylabel('$z/R_{\\rm NS}$')
        if title is None:
            t = data['time']
            title = f'$t = {t:.1f}\\,$'+'$[R_{\\rm NS}/c]$'
    

    ax.set_title(title)
    
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
    
    print("Finished plotting.")
    print('===============================')
    
    return collection, fig, ax



def interpolate_var_to_grid(data, var, Ngrid_x=2048, Ngrid_y=2048, method='nearest', x_range=None, y_range=None):
    """
    Interpolates the specified variable from cell center data onto a uniform 2D grid.

    Parameters:
    - data (dict): A dictionary containing the data, including 'center_x', 'center_y', and the variable to interpolate.
    - var (str): The key in the `data` dictionary corresponding to the variable to be interpolated.
    - Ngrid_x (int): The number of grid points along the x-axis (default is 2048).
    - Ngrid_y (int): The number of grid points along the y-axis (default is 2048).
    - method (str): The interpolation method to use ('nearest', 'linear', or 'cubic'; default is 'nearest').
    - x_range (tuple, optional): A tuple (xmin, xmax) to limit the interpolation to the specified x bounds. If None, no limits are applied.
    - y_range (tuple, optional): A tuple (ymin, ymax) to limit the interpolation to the specified y bounds. If None, no limits are applied.

    Returns:
    - tuple: A tuple containing the grid points in the x-direction, grid points in the y-direction,
              and the interpolated variable on the uniform grid.
    """
    # print('===============================')
    # print(f"Started interpolating")
    # start_time = time.time()

    center_x, center_y = data["center_x"], data["center_y"]
    # Create initial mask for both x and y
    mask = np.ones(center_x.shape, dtype=bool)

    # Apply spatial filtering based on the provided x_range
    if x_range is not None:
        x_mask = (center_x >= x_range[0]) & (center_x <= x_range[1])
        mask &= x_mask  # Combine with the overall mask

    # Apply spatial filtering based on the provided y_range
    if y_range is not None:
        y_mask = (center_y >= y_range[0]) & (center_y <= y_range[1])
        mask &= y_mask  # Combine with the overall mask

    # Filter the center_x, center_y, and variable data based on the combined mask
    filtered_center_x = center_x[mask]
    filtered_center_y = center_y[mask]
    filtered_var_data = var[mask]  # Ensure var data is filtered according to the same mask

    # Create a uniform grid based on the range of filtered x and y
    grid_x, grid_y = np.linspace(filtered_center_x.min(), filtered_center_x.max(), Ngrid_x), \
                     np.linspace(filtered_center_y.min(), filtered_center_y.max(), Ngrid_y)
    grid_x, grid_y = np.meshgrid(grid_x, grid_y)

    # Interpolate point data onto the uniform grid
    if method == 'linear':
        # Using LinearNDInterpolator for faster linear interpolation
        interpolator = LinearNDInterpolator((filtered_center_x, filtered_center_y), filtered_var_data)
        interpolated_var = interpolator(grid_x, grid_y)

    else:
        interpolated_var = griddata((filtered_center_x, filtered_center_y), filtered_var_data, (grid_x, grid_y), method=method)

    # Fill NaNs if using linear interpolation
    if method == 'linear':
        interpolated_var = fill_nan(interpolated_var)

    # end_time = time.time()
    # elapsed_time = end_time - start_time
    # print(f"Finished interpolating")
    # print(f"Time taken to interpolate: {elapsed_time:.4f} seconds")
    # print('===============================')

    return grid_x, grid_y, interpolated_var




def determine_extend_from_plot(plot_obj):
    """
    Determines the 'extend' parameter for the colorbar based on the data of the plot object.

    Parameters:
    - plot_obj (ScalarMappable): The ScalarMappable object from which to determine the color range.

    Returns:
    - str: A string indicating which end(s) of the colorbar should be extended ('neither', 'both', 'min', 'max').
    """

    norm = plot_obj.norm
    vmin = norm.vmin
    vmax = norm.vmax
    data = plot_obj.get_array()
    
    if np.any(data < vmin) and np.any(data > vmax):
        return 'both'
    elif np.any(data < vmin):
        return 'min'
    elif np.any(data > vmax):
        return 'max'
    else:
        return 'neither'
    
    
    


def smooth_vect_pot(data, Ngrid_x = 2048, Ngrid_y = 2048, sigma=5, x_range=None, y_range=None):
    """
    Interpolates the magnetic fields Bx and By, integrates them to obtain the vector potential Az, 
    and applies Gaussian smoothing to the result.

    Parameters:
    - data (dict): A dictionary containing the magnetic field data with keys 'b1' and 'b2'.
    - Ngrid_x (int): The number of grid points along the x-axis (default is 2048).
    - Ngrid_y (int): The number of grid points along the y-axis (default is 2048).
    - sigma (float): The standard deviation for Gaussian smoothing (default is 5).
    - method (str): The interpolation method to use ('nearest', 'linear', or 'cubic'; default is 'nearest').

    Returns:
    - ndarray: The smoothed vector potential Az.
    """

    results, (grid_x, grid_y) = fast_interpolate_multiple_vars_to_grid(data, [data['b1'],data['b2']], Ngrid_x=Ngrid_x, Ngrid_y=Ngrid_y, x_range=x_range, y_range=y_range)
    Bx_interp, By_interp = results[0], results[1]

    F = cumtrapz(Bx_interp, grid_y, axis=1, initial=0)        
    G = cumtrapz(-By_interp, grid_x, axis=0, initial=0) - F        
    Az_computed = F + G

    sigma = sigma  
    Az_smooth = gaussian_filter(Az_computed, sigma=sigma)

    Az_computed = Az_smooth
    
    
    # print('all:', Az_computed)
    # print('min:', np.min(Az_computed))
    # print('max:', np.max(Az_computed))

    return grid_x, grid_y, Az_computed



def plot_cells_fast(data, variable, 
                    fig=None, ax=None,
                    figsize=(6.4, 4.8), vmin=None, vmax=None,
                    cmap=None, x_range=None, y_range=None, label='',
                    title=None, save_path=None, dpi=300, colorbar=True,
                    orientation='vertical', location='right', pad=0.1,
                    use_log_norm=False, scaling=1, k_neighbors=4,
                    leafsize=16, verbose=False, aspect=20, shrink=95):
    """
    Fast plotting using interpolation to uniform grid and imshow.
    Accepts same variable types as plot_cells but uses fast interpolation.
    
    Parameters:
    -----------
    data : dict
        Data dictionary from read_HDF5
    variable : str, ndarray, or callable
        - str: variable name from data
        - ndarray: 1D array of length total_cells or cells_per_block
        - callable(data) -> array or callable(block_vars) -> (nx,ny) array
    fig, ax : matplotlib objects
        Figure and axis (created if None)
    figsize : tuple
        Figure size if creating new figure
    vmin, vmax : float
        Color scale limits
    cmap : str or colormap
        Colormap for plotting
    x_range, y_range : tuple
        (min, max) bounds for spatial filtering
    label : str
        Colorbar label
    title : str
        Plot title (defaults to time)
    save_path : str
        Path to save figure
    dpi : int
        DPI for saved figure
    colorbar : bool
        Whether to add colorbar
    orientation, location : str
        Colorbar position
    pad : float
        Colorbar padding
    use_log_norm : bool
        Use logarithmic color scale
    scaling : int
        Resolution scaling factor (1=full resolution, 2=half, etc.)
    k_neighbors : int
        Number of neighbors for KD-tree interpolation
    leafsize : int
        Leaf size for KD-tree
    verbose : bool
        Print timing information
    
    Returns:
    --------
    im : image object
    fig : figure object
    ax : axis object
    """
    
    coord_system = data['coord_system']

    if verbose:
        print('===============================')
        print(f"Started fast plotting ({coord_system}).")
    start_time = time.time()
    
    # Ensure block structure if needed for callable variables
    if callable(variable) and 'block_data' not in data:
        # Try calling with full data first
        try:
            test_result = variable(data)
            # If successful, we don't need block structure
        except:
            # Need block structure for per-block callable
            data = get_block_structure(data, data['metadata']['variable_names'])
    
    # Prepare field data based on variable type
    if isinstance(variable, str):
        # Simple variable name
        if variable not in data:
            raise ValueError(f"Variable '{variable}' not found in data")
        field_data = data[variable]
        
    elif isinstance(variable, (np.ndarray, list, tuple)):
        # Direct array input
        field_data = np.asarray(variable).ravel()
        total_cells = data['metadata']['total_cells']
        cells_per_block = data['metadata']['active_cells_per_block']
        
        if field_data.size == cells_per_block:
            # Single block data - tile to all cells
            n_blocks = data['metadata']['nblocks_processed']
            field_data = np.tile(field_data, n_blocks)
        elif field_data.size != total_cells:
            raise ValueError(f"Array size {field_data.size} doesn't match total cells {total_cells}")
            
    elif callable(variable):
        # Function to compute derived variable
        try:
            # Try calling with full data dictionary
            result = variable(data)
            field_data = np.asarray(result).ravel()
        except:
            # Fall back to per-block computation
            if 'block_data' not in data:
                data = get_block_structure(data, data['metadata']['variable_names'])
            
            block_info = data['block_data']
            active_size = block_info['active_size']
            all_values = []
            
            for block_idx in range(block_info['nblocks']):
                block_vars = block_info['block_cell_data'][block_idx]
                block_result = variable(block_vars)
                
                if np.asarray(block_result).shape != (active_size, active_size):
                    raise ValueError(f"Block function must return ({active_size},{active_size}) array")
                
                all_values.append(np.asarray(block_result).ravel())
            
            field_data = np.concatenate(all_values)
    else:
        raise ValueError("Variable must be string, array, or callable")
    
    # Determine grid resolution based on data density and scaling
    center_x = data['center_x']
    center_y = data['center_y']
    
    
    if coord_system.lower() == 'cartesian':
        # Apply range filtering for resolution calculation
        mask = np.ones(len(center_x), dtype=bool)
        if x_range is not None:
            mask &= (center_x >= x_range[0]) & (center_x <= x_range[1])
        if y_range is not None:
            mask &= (center_y >= y_range[0]) & (center_y <= y_range[1])
            
    elif coord_system.lower() == 'polar':
        
        radii = np.sqrt(center_x**2 + center_y**2)
        rmin = np.min(radii)
        rmax = np.max(radii)
                
        # Create masks for x_range, y_range, and rmin
        if x_range is not None:
            x_min = x_range[0]
            x_max = x_range[1]
            x_mask = (center_x >= x_range[0]) & (center_x <= x_range[1])
        else:
            x_min = rmin
            x_max = rmax
            x_mask = np.ones_like(center_x, dtype=bool)

        if y_range is not None:
            y_min = y_range[0]
            y_max = y_range[1]
            y_mask = (center_y >= y_range[0]) & (center_y <= y_range[1])
        else:
            y_min = -rmax
            y_max = rmax
            y_mask = np.ones_like(center_x, dtype=bool)
            
        r_mask = radii >= rmin

        # Combine all masks
        mask = x_mask & y_mask & r_mask
        

    x_filtered = center_x[mask]
    y_filtered = center_y[mask]
    
    # Estimate grid resolution from data density
    if len(x_filtered) > 1:
        # Find typical cell spacing
        x_sorted = np.sort(np.unique(x_filtered))
        y_sorted = np.sort(np.unique(y_filtered))
        
        if len(x_sorted) > 1:
            dx_min = np.min(np.diff(x_sorted))
        else:
            dx_min = (x_filtered.max() - x_filtered.min()) / 100
            
        if len(y_sorted) > 1:
            dy_min = np.min(np.diff(y_sorted))
        else:
            dy_min = (y_filtered.max() - y_filtered.min()) / 100
        
        # Calculate grid dimensions
        x_span = x_filtered.max() - x_filtered.min()
        y_span = y_filtered.max() - y_filtered.min()
        
        nx = max(2, int(x_span / dx_min)) // scaling
        ny = max(2, int(y_span / dy_min)) // scaling
        
        # Cap at reasonable maximum
        nx = min(nx, 4096 // scaling)
        ny = min(ny, 4096 // scaling)
    else:
        # Fallback for very sparse data
        nx = ny = 512 // scaling
    
    if verbose:
        print(f"Grid resolution: {nx} x {ny}")
    
    # Fast interpolation to uniform grid
    interpolated_field, (grid_x, grid_y) = fast_interpolate_var_to_grid(
        data, field_data,
        Ngrid_x=nx, Ngrid_y=ny,
        x_range=x_range, y_range=y_range,
        k_neighbors=k_neighbors,
        leafsize=leafsize,
        verbose=verbose
    )
    
    
    if coord_system.lower() == 'polar':
        # Create a mask for the half-circle (if applicable)
        mask1 = grid_x**2 + grid_y**2 <= rmax**2
        mask2 = grid_x**2 + grid_y**2 > rmin**2
        
        mask = mask1 & mask2

        interpolated_field = np.ma.masked_where(~mask, interpolated_field)

    
    # Determine color limits
    if vmin is None:
        vmin = np.nanmin(interpolated_field)
    if vmax is None:
        vmax = np.nanmax(interpolated_field)
    
    # Create figure if needed
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    
    # Set up color normalization
    if use_log_norm:
        if vmin <= 0:
            print(f"Warning: LogNorm requested but vmin={vmin} <= 0. Using linear scale.")
            norm = None
        else:
            norm = LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = None
    
    # Get extent for imshow
    extent = [grid_x.min(), grid_x.max(), grid_y.min(), grid_y.max()]
    
    # Plot using imshow
    if cmap is None:
        cmap = 'viridis'
        
    if norm is not None:
        im = ax.imshow(interpolated_field.T,  # Transpose for correct orientation
                      origin='lower',
                      extent=extent,
                      cmap=cmap,
                      norm=norm,
                      aspect='equal',
                      interpolation='nearest')
    else:
        im = ax.imshow(interpolated_field.T,  # Transpose for correct orientation
                      origin='lower',
                      extent=extent,
                      vmin=vmin,
                      vmax=vmax,
                      cmap=cmap,
                      aspect='equal',
                      interpolation='nearest')
    
    # Add colorbar
    if colorbar:
        # Determine extend parameter
        extend_type = 'neither'
        if np.any(interpolated_field < vmin):
            extend_type = 'min'
        if np.any(interpolated_field > vmax):
            extend_type = 'max'
        if np.any(interpolated_field < vmin) and np.any(interpolated_field > vmax):
            extend_type = 'both'
        
        plt.colorbar(im, ax=ax, extend=extend_type,
                    label=label, orientation=orientation,
                    location=location, pad=pad,
                    aspect=aspect, shrink=shrink)
    
    if coord_system.lower() == 'cartesian':
        ax.set_xlabel('$x/L$')
        ax.set_ylabel('$y/L$')
        if title is None:
            t = data['time']
            title = f'$t/t_c = {t:.1f}$'
    elif coord_system.lower() == 'polar':
        ax.set_xlabel('$x/R_{\\rm NS}$')
        ax.set_ylabel('$z/R_{\\rm NS}$')
        if title is None:
            t = data['time']
            title = f'$t = {t:.1f}\\,$'+'$[R_{\\rm NS}/c]$'
    
    
    ax.set_title(title)
    
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_aspect('equal')
    
    # Save if requested
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
    
    if verbose:
        total_time = time.time() - start_time
        print(f"Total plotting time: {total_time:.4f} seconds")
        print('===============================')
    
    return im, fig, ax




def rescale_axis_labels(ax, xscale=1.0, yscale=1.0, xunit='', yunit=''):
    """
    Handy function for relabeling axes labels.
    """
    if xscale != 1.0:
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x * xscale:g}{xunit}"))
    if yscale != 1.0:
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y * yscale:g}{yunit}"))



def solve_vartheta_newton(theta, h, tol=1e-12, maxiter=100):
    """
    Vectorized solver for vartheta + (h/2) sin(2*vartheta) = theta.
    Uses Newton's method with full array operations - no Python loops.
    Falls back to Brent's method only for non-converged elements.
    
    Parameters:
    -----------
    theta : float or array-like
        Input angle(s)
    h : float
        Parameter in the equation
    tol : float
        Convergence tolerance
    maxiter : int
        Maximum iterations for Newton's method
        
    Returns:
    --------
    vartheta : float or ndarray
        Solution(s) matching input shape
    """
    scalar_input = np.isscalar(theta)
    thetas = np.atleast_1d(theta).astype(float)
    
    # Initial guess (vectorized)
    v = thetas - 0.5 * h * np.sin(2.0 * thetas)
    
    # Precompute constant
    h_half = 0.5 * h
    
    # Newton iterations (fully vectorized)
    for _ in range(maxiter):
        sin_2v = np.sin(2.0 * v)
        cos_2v = np.cos(2.0 * sin_2v)  # reuse for efficiency
        
        # f(v) = v + (h/2)*sin(2v) - theta
        f = v + h_half * sin_2v - thetas
        
        # f'(v) = 1 + h*cos(2v)
        df = 1.0 + 2.0 * h_half * np.cos(2.0 * v)
        
        # Newton update
        delta = f / df
        v = v - delta
        
        # Check convergence
        if np.all(np.abs(delta) < tol):
            break
    
    # Handle non-converged elements with Brent (rarely needed with good initial guess)
    mask = np.abs(v + h_half * np.sin(2.0 * v) - thetas) > tol
    if np.any(mask):
        for i in np.where(mask)[0]:
            th = thetas[i]
            a, b = th - np.pi, th + np.pi
            try:
                res = root_scalar(lambda x: x + h_half*np.sin(2.0*x) - th, 
                                bracket=[a, b], method='brentq', xtol=tol)
                if res.converged:
                    v[i] = res.root
            except:
                pass  # Keep Newton result
    
    return v.item() if scalar_input else v



def mks_to_spherical_vector(data, v_mks, h=0.0, R0=0.0, vartheta=None):
    """
    Convert vector components from Modified Kerr-Schild (MKS) coordinates 
    to normal spherical Kerr-Schild coordinates.
    
    Based on the transformation described in Porth et al. (2017):
    (arXiv:1611.09720)

    IMPORTANT TO NOTE THE CODE FOR SOME REASON USES THE ONE WITHOUT THE ANALYTICAL TRANSFORMATION EQN(50)
    CHECK CODE AND USE AT YOUR OWN RISK (mod_coord_mks.t)
    
    """
    x, z = data['center_x'], data['center_y']
    y = 0 * x
    r = np.sqrt(x**2 + z**2 + y**2)
    R = np.sqrt(x**2 + y**2)

    theta =  np.pi/2 - np.arctan2(z,R)
 
    v_s, v_theta_mks, v_phi = v_mks
    
    # eqn 54-55
    s = np.log(r - R0)
    dr_ds = np.exp(s) 
    
    # eqn 50

    if vartheta is None:
        vartheta = solve_vartheta_newton(theta, h=h)
    dtheta_dtheta_mks = 1 + h * np.cos(2 * vartheta)
    

    v_r = v_s * dr_ds
    v_theta = v_theta_mks * dtheta_dtheta_mks
    v_phi_sph = v_phi  
    
    return v_r, v_theta, v_phi_sph


def spherical_to_cart_vector(data, vsph):
    
    x = data['center_x']
    z = data['center_y']

    x, z = data['center_x'], data['center_y']
    y = 0 * x
    r = np.sqrt(x**2 + z**2 + y**2)
    R = np.sqrt(x**2 + y**2)
    costheta = z/r
    sintheta = R/r
    cosphi = x/R
    sinphi = y/R
    
    vr, vth, vphi = vsph
    
    vx = vr * sintheta * cosphi + vth * costheta * cosphi - vphi * sinphi
    vy = vr * sintheta * sinphi + vth * costheta * sinphi + vphi * cosphi
    vz = vr * costheta - vth * sintheta
    
    return vx, vy, vz
    
    


def smooth_flux_fcn_test(data, Ngrid_x = 2048, Ngrid_y = 2048, sigma=5,
                    x_range=None, y_range=None, k_neighbors=4,
                    leafsize=16,verbose=False, 
                    mks=False, h=0, R0=0, vartheta=None):
    """
    Interpolates the magnetic fields Bx and Bz, integrates them to obtain the flux function Psi, 
    and applies Gaussian smoothing to the result.

    Parameters:
    - data (dict): A dictionary containing the magnetic field data with keys 'b1' and 'b2'.
    - Ngrid_x (int): The number of grid points along the x-axis (default is 2048).
    - Ngrid_y (int): The number of grid points along the y-axis (default is 2048).
    - sigma (float): The standard deviation for Gaussian smoothing (default is 5).
    - method (str): The interpolation method to use ('nearest', 'linear', or 'cubic'; default is 'nearest').

    Returns:
    - ndarray: The smoothed flux function Psi.
    """
    
    x, z = data['center_x'], data['center_y']
    y = 0 * x
    r = np.sqrt(x**2 + z**2 + y**2)
    R = np.sqrt(x**2 + y**2)
    costheta = z/r
    sintheta = R/r
    cosphi = x/R
    sinphi = y/R
    
    rmin = np.min(r)
    rmax = np.max(r)
    
    theta =  np.pi/2 - np.arctan2(z,R)

        
    if mks:
        # convert from mks to bl and then contravarient to physical magnetic field components
        print('Using MKS transformation with ' + f'h={h}, and R0={R0}')
        if x_range == None or y_range == None:
            print('You will want to have some xmax and zmin and zmax if you go far in log r, or resolution will be bad')
        if vartheta is None:
            vartheta = solve_vartheta_newton(theta, h=h)
        br, bth, bphi = mks_to_spherical_vector(data,[data['b1'],data['b2'],data['b3']], h=h, R0=R0, vartheta=vartheta)
        br, bth, bphi = br, r*bth, r * sintheta * bphi
    else:
        print('Treating as spherical.')
        # convert from contravarient to physical magnetic field components
        br = data['b1']
        bth = r * data['b2']
        bphi = r * sintheta * data['b3']
    
    
    bx = br * sintheta * cosphi + bth * costheta * cosphi - bphi * sinphi
    by = br * sintheta * sinphi + bth * costheta * sinphi + bphi * cosphi
    bz = br * costheta - bth * sintheta
    
    interp_bx, (grid_x, grid_z) = fast_interpolate_var_to_grid(
        data, bx,
        Ngrid_x=Ngrid_x, Ngrid_y=Ngrid_y,
        x_range=x_range, y_range=y_range,
        k_neighbors=k_neighbors,
        leafsize=leafsize,
        verbose=verbose
    )
    interp_bz, (grid_x, grid_z) = fast_interpolate_var_to_grid(
        data, bz,
        Ngrid_x=Ngrid_x, Ngrid_y=Ngrid_y,
        x_range=x_range, y_range=y_range,
        k_neighbors=k_neighbors,
        leafsize=leafsize,
        verbose=verbose
    )
    interp_bx, interp_bz, grid_x, grid_z = interp_bx.T, interp_bz.T, grid_x.T, grid_z.T
    nz, nx = interp_bx.shape

    # Extract radial coordinate R (bottom row)
    R = grid_x[0, :]

    # Compute uniform grid spacings
    dR = R[1] - R[0]
    dz = grid_z[1, 0] - grid_z[0, 0]

    # 1) Baseline R*A_phi at z0: integrate R'*Bz dR
    RA0 = np.zeros(nx)
    RA0[1:] = np.cumsum(R[:-1] * interp_bz[0, :-1] * dR)

    # 2) Cumulative z-integration of R*Br
    #    integral from z0 to z: \int R * Br dz
    R_factor = np.broadcast_to(R, (nz, nx))
    integrand = R_factor * interp_bx
    cum_RBr = np.vstack((
        np.zeros((1, nx)),
        np.cumsum(integrand[:-1, :] * dz, axis=0)
    ))

    # 3) Full R*A_phi on grid
    RA = RA0[None, :] - cum_RBr

    # 4) Compute A_phi = RA / R (mask axis R=0)
    A_phi = np.zeros_like(RA)
    mask = grid_x > 0
    A_phi[mask] = RA[mask] / grid_x[mask]

    # 5) Flux Psi = R * A_phi (equals RA)
    Psi = RA
    Psi = gaussian_filter(Psi, sigma=sigma) # smooth out noise
    
    
    
    mask1 = grid_x**2 + grid_z**2 <= rmax**2
    mask2 = grid_x**2 + grid_z**2 > rmin**2
    
    mask = mask1 & mask2

    Psi = np.ma.masked_where(~mask, Psi)

    
    return grid_x, grid_z, Psi




