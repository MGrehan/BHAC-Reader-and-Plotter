

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


def read_HDF5(filename, variable_names=None, blocks=None):
    """
    Extract cell-centered simulation data from HDF5 file.
    Always returns flattened arrays for maximum memory efficiency.
    
    Parameters:
    -----------
    filename : str
        Path to HDF5 file
    variable_names : list, optional
        Names for variables. If None, uses default names.
    blocks : array-like, optional
        Specific block indices to read. If None, reads all blocks.
        
    Returns:
    --------
    dict : Data dictionary with flattened arrays and metadata
    """
    
    with h5py.File(filename, 'r') as f:
        # Bulk read datasets once
        cc_vars = f['Cell-centered Variables'][:]   # (nblocks, nvars, 12, 12)
        grid = f['Grid'][:]                         # (nblocks, 81, 2)
        levels = f['Levels'][:]                     # (nblocks,)
        
        # File attributes
        GZ = int(f.attrs.get('GhostZones', 2))
        time = f.attrs.get('Time', 0.0)
        iteration = f.attrs.get('Iteration', 0)
        
    # Work outside of file context (now all NumPy arrays in memory)
    nblocks, nvars, block_size_x, block_size_y = cc_vars.shape
    active_size = block_size_x - 2 * GZ  # Should be 8x8 active cells
    
    print(f"Found {nblocks} blocks, {nvars} variables")
    print(f"Ghost zones: {GZ}, Active cells per block: {active_size}x{active_size}")
    print(f"Time: {time}, Iteration: {iteration}")
    
    # Determine which blocks to process
    if blocks is None:
        block_indices = np.arange(nblocks)
    else:
        block_indices = np.array(blocks, dtype=int)
        
    # Set up variable names
    if variable_names is None:
        known_names = [
            'rho', 'u1', 'u2', 'u3', 'p',
            'b1', 'b2', 'b3', 'e1', 'e2', 'e3', 's',
            'lfac', 'xi'
        ]
        variable_names = [
            known_names[i] if i < len(known_names) else f'var{i}'
            for i in range(nvars)
        ]
    elif len(variable_names) != nvars:
        raise ValueError(f"Number of variable names ({len(variable_names)}) "
                         f"must match number of variables in data ({nvars})")
    
    # Base result dictionary
    result = {
        'time': time,
        'metadata': {
            'iteration': iteration,
            'nblocks_processed': len(block_indices),
            'total_cells': len(block_indices) * active_size * active_size,
            'active_cells_per_block': active_size * active_size,
            'ghost_zones': GZ,
            'levels': levels if blocks is None else levels[block_indices],
            'variable_names': variable_names,
            'filename': filename
        }
    }
    
    # Create flattened arrays
    total_cells = len(block_indices) * active_size * active_size
    center_x = np.zeros(total_cells)
    center_y = np.zeros(total_cells)
    data_arrays = {name: np.zeros(total_cells) for name in variable_names}
    
    cell_idx = 0
    
    for b in block_indices:
        # Extract grid for this block - reshape from (81,2) to (9,9,2)
        g = grid[b]
        s = int(np.sqrt(g.shape[0]))  # Should be 9
        G = g.reshape(s, s, 2)
        
        Xn = G[:, :, 0]
        Yn = G[:, :, 1]
        
        # Cell centers (8x8 active cells)
        Xc = 0.5 * (Xn[:-1, :-1] + Xn[1:, 1:])
        Yc = 0.5 * (Yn[:-1, :-1] + Yn[1:, 1:])
        
        cells_in_block = active_size * active_size
        for var_idx, var_name in enumerate(variable_names):
            var_data = cc_vars[b, var_idx, GZ:-GZ, GZ:-GZ]
            data_arrays[var_name][cell_idx:cell_idx + cells_in_block] = var_data.ravel()
        
        center_x[cell_idx:cell_idx + cells_in_block] = Xc.ravel()
        center_y[cell_idx:cell_idx + cells_in_block] = Yc.ravel()
        
        cell_idx += cells_in_block
    
    # Add flattened data to result
    result['center_x'] = center_x
    result['center_y'] = center_y
    result.update(data_arrays)
    
    return result


def get_block_structure(data, variables=None):
    """
    Reconstruct block structure from flattened data.
    This is the reverse operation of flattening - takes flattened arrays and 
    reshapes them back into block structure for exact cell plotting.
    
    Parameters:
    -----------
    data : dict
        Data dictionary with flattened arrays (preserve_blocks=False)
        Must contain grid information in metadata to reconstruct blocks
    variables : list, optional
        List of variable names to unflatten. If None, uses all variables.
        
    Returns:
    --------
    dict : Dictionary with block structure added
    """
    if 'block_data' in data:
        print("Warning: Data already has block structure")
        return data
    
    # We need to reconstruct the block structure from the HDF5 file
    # since we need the grid information
    filename = data['metadata']['filename']
    nblocks_processed = data['metadata']['nblocks_processed']
    ghost_zones = data['metadata']['ghost_zones']
    variable_names = data['metadata']['variable_names']
    
    if variables is None:
        variables = variable_names
    
    print(f"Reconstructing block structure from {filename}...")
    
    # Re-read only the grid information (lightweight)
    with h5py.File(filename, 'r') as f:
        grid = f['Grid'][:]  # (nblocks, 81, 2)
        levels = f['Levels'][:]  # (nblocks,)
    
    active_size = int(np.sqrt(data['metadata']['active_cells_per_block']))  # Should be 8
    
    # Reconstruct block data structure
    block_grid_nodes = []
    block_cell_data = []
    block_boundaries = []
    
    cells_per_block = active_size * active_size
    
    for block_idx in range(nblocks_processed):
        # Extract grid for this block
        g = grid[block_idx]
        s = int(np.sqrt(g.shape[0]))  # Should be 9
        G = g.reshape(s, s, 2)
        block_grid_nodes.append(G.copy())
        
        # Extract block boundaries (corners)
        corners = np.array([
            [G[0, 0, 0], G[0, 0, 1]],      # Bottom-left
            [G[-1, 0, 0], G[-1, 0, 1]],    # Bottom-right  
            [G[-1, -1, 0], G[-1, -1, 1]],  # Top-right
            [G[0, -1, 0], G[0, -1, 1]]     # Top-left
        ])
        block_boundaries.append(corners)
        
        # Unflatten variables for this block
        start_idx = block_idx * cells_per_block
        end_idx = start_idx + cells_per_block
        
        block_vars = {}
        for var in variables:
            if var in data:
                # Reshape flattened data back to (active_size, active_size)
                flattened_block = data[var][start_idx:end_idx]
                block_vars[var] = flattened_block.reshape(active_size, active_size)
        
        block_cell_data.append(block_vars)
    
    # Create a new data dictionary with block structure added
    result = data.copy()
    result['block_data'] = {
        'nblocks': nblocks_processed,
        'active_size': active_size,
        'ghost_zones': ghost_zones,
        'block_indices': list(range(nblocks_processed)),
        'block_grid_nodes': block_grid_nodes,
        'block_cell_data': block_cell_data,
        'block_boundaries': np.array(block_boundaries),
        'cells_per_block': cells_per_block
    }
    
    print(f"Reconstructed block structure for {nblocks_processed} blocks")
    return result


def plot_cells(data, variable, 
               fig=None, ax=None,
               figsize=(6.4, 4.8), vmin=None, vmax=None, 
               cmap=None, x_range=None, y_range=None, label='',
               title=None, save_path=None, dpi=300, colorbar=True,
               orientation='vertical', location='right', pad=0.1,
               edgecolors=None, linewidths=0.1, use_log_norm=False,
               derived_function=None):
    """
    Plot exact simulation cells using preserved block structure from extracted data.
    This avoids re-reading the HDF5 file and shows true mesh geometry.
    
    """
    
    # Automatically get block structure if not present
    if 'block_data' not in data:
        print("Reconstructing block structure for plotting...")
        # Determine which variables we need based on the variable parameter
        if callable(variable):
            # For function case, we need to reconstruct all variables since we don't know
            # which ones the function will use. User can optimize by passing specific variables.
            needed_vars = data['metadata']['variable_names']
        else:
            # For simple variable case, we only need that variable
            needed_vars = [variable] if isinstance(variable, str) else data['metadata']['variable_names']
        
        data = get_block_structure(data, needed_vars)

    # Handle deprecated parameter
    if derived_function is not None and callable(variable):
        raise ValueError("Cannot specify both 'variable' as function and 'derived_function'. Use only 'variable'.")
    if derived_function is not None:
        variable = derived_function
        print("Warning: 'derived_function' parameter is deprecated. Pass function as 'variable' instead.")
    
    print('===============================')
    print(f"Started plotting.")
    
    block_info = data['block_data']
    
    all_polygons = []
    all_values = []
    
    for block_idx in range(block_info['nblocks']):
        G = block_info['block_grid_nodes'][block_idx]   # (9,9,2)
        Xn = G[:, :, 0]
        Yn = G[:, :, 1]
        
        block_vars = block_info['block_cell_data'][block_idx]
        active_size = block_info['active_size']
        
        # Variable extraction
        if callable(variable):
            var_data = variable(block_vars)
            if var_data.shape != (active_size, active_size):
                raise ValueError(f"Function must return array of shape ({active_size}, {active_size}), got {var_data.shape}")
        else:
            if variable not in block_vars:
                available_vars = list(block_vars.keys())
                raise ValueError(f"Variable '{variable}' not found. Available variables: {available_vars}")
            var_data = block_vars[variable]
        
        # Vectorized polygon construction
        x0 = Xn[:-1, :-1].ravel()
        y0 = Yn[:-1, :-1].ravel()
        x1 = Xn[1:, :-1].ravel()
        y1 = Yn[1:, :-1].ravel()
        x2 = Xn[1:, 1:].ravel()
        y2 = Yn[1:, 1:].ravel()
        x3 = Xn[:-1, 1:].ravel()
        y3 = Yn[:-1, 1:].ravel()
        
        block_polygons = np.stack([
            np.stack([x0, y0], axis=1),
            np.stack([x1, y1], axis=1),
            np.stack([x2, y2], axis=1),
            np.stack([x3, y3], axis=1)
        ], axis=1)
        
        block_values = var_data.ravel()
        
        # Vectorized masking by x_range / y_range
        mask = np.ones(len(block_polygons), dtype=bool)
        if x_range is not None:
            xmin, xmax = x_range
            mask &= (block_polygons[:, :, 0].min(axis=1) <= xmax) & (block_polygons[:, :, 0].max(axis=1) >= xmin)
        if y_range is not None:
            ymin, ymax = y_range
            mask &= (block_polygons[:, :, 1].min(axis=1) <= ymax) & (block_polygons[:, :, 1].max(axis=1) >= ymin)
        
        # Keep only in-range cells
        all_polygons.append(block_polygons[mask])
        all_values.append(block_values[mask])
    
    # Concatenate once
    all_polygons = np.concatenate(all_polygons, axis=0)
    all_values = np.concatenate(all_values, axis=0)
    
    if len(all_polygons) == 0:
        raise ValueError("No cells found in specified range")
    
    # Color limits
    if vmin is None:
        vmin = np.nanmin(all_values)
    if vmax is None:
        vmax = np.nanmax(all_values)
    if cmap is None:
        cmap = 'viridis'
    
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    
    # Color normalization
    if use_log_norm:
        if vmin <= 0:
            print(f"Warning: LogNorm requested but vmin={vmin} <= 0. Using linear scale.")
            norm = None
        else:
            norm = LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = None
    
    if norm is not None:
        poly_collection = PolyCollection(all_polygons, norm=norm, array=all_values,
                                         cmap=cmap, edgecolors=edgecolors,
                                         linewidths=linewidths)
    else:
        poly_collection = PolyCollection(all_polygons, array=all_values,
                                         cmap=cmap, edgecolors=edgecolors,
                                         clim=(vmin, vmax), linewidths=linewidths)
    
    ax.add_collection(poly_collection)
    
    # Colorbar
    if colorbar:
        extend_type = 'neither'
        if np.any(all_values < vmin):
            extend_type = 'min'
        if np.any(all_values > vmax):
            extend_type = 'max'
        if np.any(all_values < vmin) and np.any(all_values > vmax):
            extend_type = 'both'
        
        plt.colorbar(poly_collection, ax=ax, extend=extend_type,
                     label=label, orientation=orientation, location=location,
                     pad=pad)
    
    # Axis limits
    all_x = all_polygons[:, :, 0].ravel()
    all_y = all_polygons[:, :, 1].ravel()
    ax.set_xlim(all_x.min(), all_x.max())
    ax.set_ylim(all_y.min(), all_y.max())
    ax.set_aspect('equal')
    ax.set_xlabel('$x/L$')
    ax.set_ylabel('$y/L$')
    
    if title is None:
        t = data['time']
        title = f'$t/t_c = {t:.1f}$'
    ax.set_title(title)
    
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
    
    print("Finished plotting.")
    print('===============================')
    
    return fig, ax, poly_collection


def plot_blocks(data, fig=None, ax=None, linewidth=0.1, color='k', 
                x_range=None, y_range=None, title=None, figsize=(6.4, 4.8)):
    """
    Plot block boundaries using flattened data.
    Automatically reconstructs block structure as needed.
    
    """
    
    # Automatically get block structure if not present
    if 'block_data' not in data:
        print("Reconstructing block structure for plotting...")
        # For block plotting, we don't need any variable data, just grid structure
        # Pass empty list to avoid unflattening variable data unnecessarily
        data = get_block_structure(data, variables=[])
    
    print('===============================')
    print("Started plotting block boundaries")
    
    block_boundaries = data['block_data']['block_boundaries']  # Shape: (nblocks, 4, 2)
    N = block_boundaries.shape[0]
    
    # Prepare segments for LineCollection
    segments = []
    for i in range(N):
        block_corners = block_boundaries[i]  # Shape: (4, 2)
        
        # Extract x and y coordinates
        x_coords = block_corners[:, 0]
        y_coords = block_corners[:, 1]
        
        # Check if the block falls within the provided spatial range
        if x_range is not None:
            if np.all(x_coords < x_range[0]) or np.all(x_coords > x_range[1]):
                continue  # Skip this block
        if y_range is not None:
            if np.all(y_coords < y_range[0]) or np.all(y_coords > y_range[1]):
                continue  # Skip this block
        
        # Create segments for each edge of the block
        # Order: bottom, right, top, left (to form closed rectangle)
        segments.extend([
            [block_corners[0], block_corners[1]],  # Bottom edge (BL -> BR)
            [block_corners[1], block_corners[2]],  # Right edge (BR -> TR)
            [block_corners[2], block_corners[3]],  # Top edge (TR -> TL)
            [block_corners[3], block_corners[0]]   # Left edge (TL -> BL)
        ])
    
    # Create figure and axis if not provided
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=figsize)
        
        
    # Create and add LineCollection
    lc = LineCollection(segments, linewidths=linewidth, colors=color)
    ax.add_collection(lc)
    
    
    
    ax.set_aspect('equal')
    ax.set_xlabel('$x/L$')
    ax.set_xlabel('$y/L$')
    
    
    # Set plot limits based on provided ranges or data bounds
    if x_range is None:
        ax.set_xlim(x_coords.min(), x_coords.max())
    else:
        ax.set_xlim(x_range)
    
    if y_range is None:
        ax.set_ylim(y_coords.min(), y_coords.max())
    else:
        ax.set_ylim(y_range)

    # Set title
    if title is None:
        t = data['time']
        title = f'$t/t_c = {t:.1f}$'
    ax.set_title(title)
    
    print("Finished plotting block boundaries")
    print('===============================')
    
    return fig, ax

