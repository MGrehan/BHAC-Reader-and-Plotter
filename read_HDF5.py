import numpy as np
import h5py
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import cmasher as cmr
import time


# filename = '/Users/michaelgrehan/Downloads/data0001.hdf5'
# variable_names = ['rho', 'u1', 'u2', 'u3', 'p', 'b1', 'b2', 'b3', 'e1', 'e2', 'e3', 's', 'lfac', 'xi',
#                   'jz', 'v1', 'v2', 'v3', 'divB', 'curlbz', 'curlvz', 'divE']
# data = read_HDF5(filename=filename, variable_names=variable_names, preserve_blocks=True)
# plot_exact_cells_from_data(data, lambda vars: vars['b1']**2+vars['b2']**2,
#                            x_range=(-0.5,0.5), y_range=(-0.5,0.5),
#                            vmin=0, vmax=1, cmap=cmr.ember, pad=0.0, label='$B_\\perp^2/B_0$')




def read_HDF5(filename, variable_names=None, blocks=None, preserve_blocks=False):
    """
    Extract cell-centered simulation data from HDF5 file and return as 1D arrays.
    Optimized to bulk-read datasets instead of per-block slicing from disk.
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
    
    # Total number of active cells
    total_cells = len(block_indices) * active_size * active_size
    
    # Initialize output arrays
    center_x = np.zeros(total_cells)
    center_y = np.zeros(total_cells)
    data_arrays = {name: np.zeros(total_cells) for name in variable_names}
    
    # If preserving blocks, also store block structure data
    if preserve_blocks:
        block_grid_nodes = []
        block_cell_data = []
        block_start_indices = []
    
    # Process selected blocks
    cell_idx = 0
    for b in block_indices:
        # Extract grid for this block - reshape from (81,2) to (9,9,2)
        g = grid[b]
        s = int(np.sqrt(g.shape[0]))  # Should be 9
        G = g.reshape(s, s, 2)
        Xn = G[:, :, 0]
        Yn = G[:, :, 1]
        
        if preserve_blocks:
            block_grid_nodes.append(G.copy())
            block_start_indices.append(cell_idx)
            
            block_vars = {
                var_name: cc_vars[b, var_idx, GZ:-GZ, GZ:-GZ].copy()
                for var_idx, var_name in enumerate(variable_names)
            }
            block_cell_data.append(block_vars)
        
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
    
    # Assemble result
    result = {
        'center_x': center_x,
        'center_y': center_y,
        'time': time,
        'metadata': {
            'iteration': iteration,
            'nblocks_processed': len(block_indices),
            'total_cells': total_cells,
            'active_cells_per_block': active_size * active_size,
            'ghost_zones': GZ,
            'levels': levels if blocks is None else levels[block_indices],
            'variable_names': variable_names,
            'filename': filename
        }
    }
    
    result.update(data_arrays)
    
    if preserve_blocks:
        result['block_data'] = {
            'nblocks': len(block_indices),
            'active_size': active_size,
            'ghost_zones': GZ,
            'block_indices': block_indices.tolist(),
            'block_grid_nodes': block_grid_nodes,
            'block_cell_data': block_cell_data,
            'block_start_indices': block_start_indices,
            'cells_per_block': active_size * active_size
        }
        print(f"Preserved block structure data for {len(block_indices)} blocks")
    
    return result








def plot_cells(data, variable, figsize=(6.4, 4.8), vmin=None, vmax=None, 
                              cmap=None, x_range=None, y_range=None, label='',
                              title=None, save_path=None, colorbar=True,
                              orientation='vertical', location='right', pad=0.1,
                              edgecolors=None, linewidths=0.1, use_log_norm=False,
                              derived_function=None):
    """
    Plot exact simulation cells using preserved block structure from extracted data.
    This avoids re-reading the HDF5 file and shows true mesh geometry.
    
    Parameters:
    -----------
    data : dict
        Data dictionary returned by read_HDF5() with preserve_blocks=True
    variable : str or callable
        Name of variable to plot, OR a function that takes block_vars dict and returns (8,8) array
    derived_function : callable, optional
        DEPRECATED: use variable parameter instead
        Function that takes a dictionary of variables and returns computed values
    figsize : tuple, optional
        Figure size (width, height)
    vmin, vmax : float, optional
        Color scale limits
    cmap : str, optional
        Colormap name
    x_range, y_range : tuple, optional
        Spatial ranges to plot (xmin, xmax), (ymin, ymax)
    label : str, optional
        Colorbar label
    title : str, optional
        Plot title
    save_path : str, optional
        Path to save figure
    colorbar : bool, optional
        Whether to show colorbar
    orientation, location, pad : optional
        Colorbar properties
    edgecolors, linewidths : optional
        Cell edge properties
    use_log_norm : bool, optional
        Whether to use logarithmic color scale
        
    Returns:
    --------
    fig, ax, poly_collection : matplotlib objects
    
    Examples:
    ---------
    # Plot a basic variable
    plot_exact_cells_from_data(data, 'b1')
    
    # Plot a derived quantity using lambda function
    plot_exact_cells_from_data(data, lambda vars: vars['b1']**2 + vars['b2']**2)
    
    # Plot a complex derived quantity
    def magnetic_pressure(vars):
        return (vars['b1']**2 + vars['b2']**2 + vars['b3']**2) / (2 * vars['rho'])
    
    plot_exact_cells_from_data(data, magnetic_pressure, label='$P_{mag}/P_0$')
    """

    if 'block_data' not in data:
        raise ValueError("Data must be extracted with preserve_blocks=True to use this function")
    
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
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    print("Finished plotting.")
    print('===============================')
    
    return fig, ax, poly_collection
