"""
BHAC VTU Reader

Author: Michael Patrick Grehan
Date: October 2024
Email: michael.grehan@mail.utoronto.ca

------------------------------------------------------------------------------------
This module provides functionality for reading and processing 2D BHAC VTU data files
that do not include GR effects. The module is designed to read 
VTU files, extract relevant data fields, and return them in a format suitable for 
numerical analysis and visualization.
------------------------------------------------------------------------------------


--------------
Main Features:
--------------
- Efficient reading of VTU files for 2D simulations
- Handles grid and solution data extraction
- Supports loading into NumPy arrays for further manipulation
- Interpolation of field data
- Plotting of vector potential contour lines and arrows
- Plotting of cell boundaries
- Plotting of block boundaries
- Plotting of raw simulation data
- 1d slicing of raw simulation data


--------------
Usage Example (interpolation):
--------------
filename = 'data0075.vtu'
data, names = fast_vtu_reader(filename, blocks=True)
Ngrid_x, Ngrid_y = 1024, 2048
grid_x, grid_y, interp_b1 = interpolate_var_to_grid(data, "b1", Ngrid_x=Ngrid_x, Ngrid_y=Ngrid_y)
_, _, interp_b2 = interpolate_var_to_grid(data, "b2", Ngrid_x=Ngrid_x, Ngrid_y=Ngrid_y)
_, _, interp_b3 = interpolate_var_to_grid(data, "b3", Ngrid_x=Ngrid_x, Ngrid_y=Ngrid_y)
_, _, interp_p = interpolate_var_to_grid(data, "p", Ngrid_x=Ngrid_x, Ngrid_y=Ngrid_y)
B2 = interp_b1**2 + interp_b2**2 + interp_b3**2
_, _, Az_computed = smooth_vect_pot(data, Ngrid_x = Ngrid_x, Ngrid_y = Ngrid_y)

fig, ax = plt.subplots()
p1 = ax.imshow(np.abs(2*interp_p/B2), cmap="hot", origin="lower",
               extent=[data['xpoint'].min(), data['xpoint'].max(), data['ypoint'].min(), data['ypoint'].max()],
               norm=LogNorm(vmax=1e2, vmin=1e-1))
xmin, xmax = -0.1, 0.1
ymin, ymax = -0.2, 0.2
n_levels = 300
contour = ax.contour(grid_x, grid_y, Az_computed, levels=n_levels, colors='w', linewidths=0.5)
for collection in contour.collections:
    for path in collection.get_paths():
        path_data = path.vertices
        x = path_data[:, 0]
        y = path_data[:, 1]
        line, = ax.plot(x, y, color=collection.get_edgecolor(), linewidth=0.5)
        line.set_visible(False)
        add_arrow(line)
ax.set_xlabel('$x/L$')
ax.set_ylabel('$y/L$')
# ax.set_aspect('equal')
cbar = fig.colorbar(p1, ax=ax, pad=0.05,  extend=determine_extend_from_plot(p1), orientation='horizontal',  location='top')
cbar.set_label('$\\beta$')
plot_cells(data, fig=fig, ax=ax, linewidth=0.25, color='w', x_range=(xmin,xmax), y_range=(ymin,ymax))
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin ,ymax)
plt.show()

--------------
Usage Example (raw data plotting):
--------------
filename = 'data0075.vtu'
data, names = fast_vtu_reader(filename, attr={'p', 'b1', 'b2', 'b3'}, blocks=False)
fig, ax = plt.subplots()
xmin, xmax = -0.01, -0.007
ymin, ymax = -0.001, 0.001
plot_raw_data_cells(data, 2*data['p']/(data['b1']**2 + data['b2']**2 + data['b3']**2), fig=fig, ax=ax, x_range=(xmin,xmax), y_range=(ymin,ymax), cmap='hot', label='$\\beta$', linewidths=0.1, edgecolors='k', orientation='horizontal',  location='top', use_log_norm=True)
ax.set_xlabel('$x/L$')
ax.set_ylabel('$y/L$')
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin ,ymax)
plt.show()

--------------
Usage Example (polar plotting):
--------------
filename = '/Users/michaelgrehan/Desktop/data0001.vtu'
data, names = fast_vtu_reader(filename, blocks=False)
xmin, xmax = 0, 20
ymin, ymax = -25, 25
Ngrid = 1024
fig, ax = plt.subplots(figsize=(5,7))
grid_x, grid_y, interp_b1 = interpolate_var_to_grid(data, "b1", Ngrid_x=Ngrid, 
                                                    Ngrid_y=Ngrid, x_range=(xmin,xmax), 
                                                    y_range=(ymin, ymax))
grid_x, grid_y, interp_b2 = interpolate_var_to_grid(data, "b2", Ngrid_x=Ngrid, 
                                                    Ngrid_y=Ngrid, x_range=(xmin,xmax), 
                                                    y_range=(ymin, ymax))
p1, _, _ = plot_polar_data_cells_continuous(data, 
                                ((data['center_x']**2 + data['center_y']**2)**(3/2)) * data['b1'], 
                                 fig=fig, ax=ax, label='$r^3 B^{x}/B_\\star$', 
                                 x_range=(xmin,xmax), y_range=(ymin, ymax), 
                                 resolution = Ngrid,
                                 colorbar=None, cmap=cmr.wildfire)
cbar = fig.colorbar(p1, ax=ax, pad=0.01,  
                    extend=determine_extend_from_plot(p1), label='$r^3 B^{x}/B_\\star$')
ax.streamplot(grid_x, grid_y, interp_b1, interp_b2, color='w', linewidth=0.75, 
              broken_streamlines=False, density=0.35)
plt.show()

--------------
Usage Example (plotting sliced data from 3D simulation):
--------------
filename = '/Users/michaelgrehan/Downloads/data_d3_x+0.00D+00_n0000.vtu'
data, names = fast_vtu_reader_ascii(filename, slicedir='z') 
fig, ax = plt.subplots()
xmin, xmax = -0.5, 0.5
ymin, ymax = -0.5, 0.5
plot_raw_data_cells(data, data['b1'], fig=fig, ax=ax, x_range=(xmin,xmax), y_range=(ymin,ymax), cmap='cmr.iceburn', label='$B_x$', linewidths=0.1, edgecolors='k', orientation='horizontal',  location='top', use_log_norm=False)
ax.set_xlabel('$x/L$')
ax.set_ylabel('$y/L$')
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin ,ymax)
plt.show() 


--------------
Usage Example (plotting sliced 3D polar data in 2D with slicing):
--------------

filename = '/Users/michaelgrehan/Downloads/data0003.vtu'
data, names = fast_vtu_reader(filename, blocks=False)
x, y, z = calculate_cell_centers_3d(data)
r = np.sqrt(x**2 + y**2 + z**2)
theta = np.arccos(z / r)
phi = np.arctan2(y, x)
t = data['time']
br = data['b1']*np.sin(theta)*np.cos(phi) + data['b2']*np.sin(theta)*np.sin(phi)+data['b3']*np.cos(theta)
fig, ax = plt.subplots(figsize=(7,6.5))
plot_interpolated_2d_slice(data,
                           br*(r**3),
                           slice_direction='xz', fig=fig, ax=ax, cmap='cmr.ember',
                           label='$\\left( r/R_{\\rm{NS}} \\right)^3 B^r/B_\\star$',
                                                      atol=5e-3,slice_position=0.0,
                           vmax=2.0, vmin=1.0)

ax.streamplot(grid_x, grid_y, interp_b1, interp_b3, color='w', linewidth=0.75, 
                broken_streamlines=False, density=0.35)  
ax.set_title(f'$t/t_c = {t:.2f}$') 
plt.show() 




--------------
Libraries Used:
--------------
- struct: For handling binary data.
- numpy: For efficient numerical operations.
- xml.etree.ElementTree: For parsing the VTU XML structure.
- time: For timing the file reading and processing steps.
- base64: For decoding base64-encoded data.
- scipy.integrate, scipy.ndimage, scipy.interpolate: For interpolation and smoothing.
- matplotlib.pyplot: For plotting the data.
"""


import struct
import numpy as np
import xml.etree.ElementTree as ET
import time
import base64
from scipy.integrate import cumtrapz
from scipy.ndimage import gaussian_filter
from scipy.interpolate import NearestNDInterpolator
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.collections import LineCollection
from matplotlib.collections import PolyCollection



def fast_vtu_reader(filename, attr='all', blocks=False):
    """
    Reads a VTU file produced by BHAC for 2D simulations.

    Parameters:
    - filename: str, path to the VTU file.
    - attr: list or 'all', attributes to extract (default is 'all').
    - blocks: bool, whether to compute block boundaries for visualization (default is False).

    Returns:
    - data: dict, containing extracted points, attributes, and calculated cell centers.
    - data_keys: list, names of the data arrays present in the file.
    """
    
    print('===============================')
    print(f"Starting to read file: {filename}")
    start_time = time.time()

    with open(filename, 'rb') as f:
        content = f.read()

    appended_data_start = content.find(b'<AppendedData encoding="raw">')
    if appended_data_start == -1:
        raise ValueError("AppendedData section not found")

    data_start = content.find(b'_', appended_data_start) + 1
    xml_content = content[:appended_data_start].decode('utf-8', errors='ignore')
    root = ET.fromstring(xml_content + '</VTKFile>')
    
    data = {}
    
    time_field = root.find(".//FieldData/DataArray[@Name='TIME']")
    if time_field is not None:
        data['time'] = float(time_field.text.strip())
        print(f"Extracted time: {data['time']}")
    else:
        print("No time field found in the file.")

    pieces = root.findall('.//Piece')
    num_pieces = len(pieces)
    print(f"Number of Pieces: {num_pieces}")

    cells_per_piece = int(pieces[0].get('NumberOfCells'))
    total_cells = cells_per_piece * num_pieces
    print(f"Cells per piece: {cells_per_piece}")
    print(f"Total number of cells: {total_cells}")

    # Get all unique DataArray names
    data_array_names = set()
    for piece in pieces:
        for data_array in piece.findall('.//DataArray'):
            data_array_names.add(data_array.get('Name'))

    # Read Points (x, y, z coordinates)
    points_data = []
    block_boundaries = []
    for piece in pieces:
        points_data_array = piece.find('.//Points/DataArray')
        if points_data_array is None:
            raise ValueError("Points data not found")

        dtype = points_data_array.get('type')
        format = points_data_array.get('format')

        if format == 'appended':
            offset = int(points_data_array.get('offset', '0'))
            size = struct.unpack('<I', content[data_start+offset:data_start+offset+4])[0]
            raw_data = content[data_start+offset+4:data_start+offset+4+size]
        elif format == 'ascii':
            raw_data = points_data_array.text.strip().split()
        else:  # Assume inline base64
            raw_data = base64.b64decode(points_data_array.text.strip())

        if dtype == 'Float32':
            parsed_data = np.frombuffer(raw_data, dtype=np.float32) if format != 'ascii' else np.array(raw_data, dtype=np.float32)
        elif dtype == 'Float64':
            parsed_data = np.frombuffer(raw_data, dtype=np.float64) if format != 'ascii' else np.array(raw_data, dtype=np.float64)
        else:
            raise ValueError(f"Unsupported data type for Points: {dtype}")
        points_data.append(parsed_data)
        
        if blocks:
            # Reshape the parsed data for this piece
            piece_points = parsed_data.reshape(-1, 3)

            # Vectorized min and max operations for this piece
            x_min, y_min = np.min(piece_points[:, :2], axis=0)
            x_max, y_max = np.max(piece_points[:, :2], axis=0)

            # Define block corners for this piece
            corners = np.array([
                [x_min, y_min],
                [x_max, y_min],
                [x_max, y_max],
                [x_min, y_max],
                [x_min, y_min]  # Close the loop
            ])

            # Create block boundaries for this piece
            piece_boundaries = np.array([corners[:-1], corners[1:]]).transpose(1, 0, 2)
            
            block_boundaries.append(piece_boundaries)

    if blocks:
        data['block_coord'] = np.array(block_boundaries)

        

    if points_data:
        points = np.concatenate(points_data).reshape(-1, 3)  # Assuming 3D points (x, y, z)
        data['xpoint'], data['ypoint'], data['zpoint'] = points[:, 0], points[:, 1], points[:, 2]
        print(f"Extracted {len(data['xpoint'])} points")



    # Handle attributes
    if attr == 'all':
        data_array_names.discard(None)
        data_array_names.discard('types')
    else:
        data_array_names = attr
        data_array_names.add('connectivity')
        data_array_names.add('offsets')



    for name in data_array_names:
        combined_data = []

        for piece in pieces:
            piece_data_array = piece.find(f".//DataArray[@Name='{name}']")
            if piece_data_array is None:
                continue

            dtype = piece_data_array.get('type')
            format = piece_data_array.get('format')

            if format == 'appended':
                offset = int(piece_data_array.get('offset', '0'))
                size = struct.unpack('<I', content[data_start+offset:data_start+offset+4])[0]
                raw_data = content[data_start+offset+4:data_start+offset+4+size]
            elif format == 'ascii':
                raw_data = piece_data_array.text.strip().split()
            else:
                raw_data = base64.b64decode(piece_data_array.text.strip())

            if dtype == 'Float32':
                parsed_data = np.frombuffer(raw_data, dtype=np.float32) if format != 'ascii' else np.array(raw_data, dtype=np.float32)
            elif dtype == 'Float64':
                parsed_data = np.frombuffer(raw_data, dtype=np.float64) if format != 'ascii' else np.array(raw_data, dtype=np.float64)
            elif dtype == 'Int32':
                parsed_data = np.frombuffer(raw_data, dtype=np.int32) if format != 'ascii' else np.array(raw_data, dtype=np.int32)
            elif dtype == 'Int64':
                parsed_data = np.frombuffer(raw_data, dtype=np.int64) if format != 'ascii' else np.array(raw_data, dtype=np.int64)
            else:
                raise ValueError(f"Unsupported data type: {dtype}")

            combined_data.append(parsed_data)

        if combined_data:
            data[name] = np.concatenate(combined_data)
            
            
 


    data["ncells"] = total_cells
    data["center_x"], data["center_y"] = calculate_cell_centers(data)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Finished reading file: {filename}")
    print(f"Time taken to read: {elapsed_time:.4f} seconds")
    print('===============================')


    return data, list(data.keys())


def calculate_cell_centers(data):
    """
    Calculate the cell centers by averaging the vertex coordinates.

    Parameters:
    - data: dict, containing vertex coordinates and connectivity information.

    Returns:
    - center_x: ndarray, x-coordinates of cell centers.
    - center_y: ndarray, y-coordinates of cell centers.
    """

    print('===============================')
    print(f"Started finding cell centers")
    start_time = time.time()
    
    x = data['xpoint']
    y = data['ypoint']
    ncells = data['ncells']
    
    offsets = data['offsets']
    connectivity = data['connectivity']
    # Create mod_conn array using broadcasting instead of a for loop
    base_conn = connectivity[:np.max(offsets)]  # Base mod_conn for the first set
    num_iterations = int(4 * ncells / np.max(offsets))  # Number of iterations

    # Use broadcasting to create mod_conn without a loop
    offsets_array = np.arange(num_iterations) * (np.max(base_conn) + 1)  # Calculate all offsets at once
    mod_conn = base_conn + offsets_array[:, None]  # Broadcast and add offsets
    
    # Flatten mod_conn to a 1D array
    mod_conn = mod_conn.ravel()[:ncells * 4]  # Only take enough entries for ncells

    # Reshape mod_conn to group cell vertices (ncells x 4)
    cell_vertices = mod_conn.reshape(ncells, 4)

    # Vectorized calculation of cell centers
    cell_centers_x = np.mean(x[cell_vertices], axis=1)
    cell_centers_y = np.mean(y[cell_vertices], axis=1)
    
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Finished finding cell centers")
    print(f"Time taken to get centers: {elapsed_time:.4f} seconds")
    print('===============================')

    return cell_centers_x, cell_centers_y


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
    print('===============================')
    print(f"Started interpolating")
    start_time = time.time()

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
    filtered_var_data = data[var][mask]  # Ensure var data is filtered according to the same mask

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

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Finished interpolating")
    print(f"Time taken to interpolate: {elapsed_time:.4f} seconds")
    print('===============================')

    return grid_x, grid_y, interpolated_var



def interpolate_vect_pot_to_grid(data, Az, Ngrid_x=2048, Ngrid_y=2048, method='nearest', x_range=None, y_range=None):
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
    print('===============================')
    print(f"Started interpolating")
    start_time = time.time()

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
    filtered_Az = Az[mask]  # Ensure var data is filtered according to the same mask

    # Create a uniform grid based on the range of filtered x and y
    grid_x, grid_y = np.linspace(filtered_center_x.min(), filtered_center_x.max(), Ngrid_x), \
                     np.linspace(filtered_center_y.min(), filtered_center_y.max(), Ngrid_y)
    grid_x, grid_y = np.meshgrid(grid_x, grid_y)

    # Interpolate point data onto the uniform grid
    if method == 'linear':
        # Using LinearNDInterpolator for faster linear interpolation
        interpolator = LinearNDInterpolator((filtered_center_x, filtered_center_y), filtered_Az)
        interpolated_var = interpolator(grid_x, grid_y)

    else:
        interpolated_var = griddata((filtered_center_x, filtered_center_y), filtered_Az, (grid_x, grid_y), method=method)

    # Fill NaNs if using linear interpolation
    if method == 'linear':
        interpolated_var = fill_nan(interpolated_var)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Finished interpolating")
    print(f"Time taken to interpolate: {elapsed_time:.4f} seconds")
    print('===============================')

    return grid_x, grid_y, interpolated_var

def fill_nan(grid):
    """
    Fills NaN values in a 2D array using nearest neighbor interpolation.

    Parameters:
    - grid (ndarray): A 2D array with NaN values that need to be filled.

    Returns:
    - ndarray: The input array with NaN values filled.
    """
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

def smooth_vect_pot(data, Ngrid_x = 2048, Ngrid_y = 2048, sigma=5, method='nearest',  x_range=None, y_range=None):
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

    grid_x, grid_y, Bx_interp = interpolate_var_to_grid(data, "b1", Ngrid_x=Ngrid_x, Ngrid_y=Ngrid_y, method=method, x_range=x_range, y_range=y_range)
    _, _, By_interp = interpolate_var_to_grid(data, "b2", Ngrid_x=Ngrid_x, Ngrid_y=Ngrid_y, method=method, x_range=x_range, y_range=y_range)

    F = cumtrapz(Bx_interp, grid_y, axis=0, initial=0)        
    G = cumtrapz(-By_interp, grid_x, axis=1, initial=0) - F        
    Az_computed = F + G
    # Enforce periodic boundary conditions for Az_computed
    # Az_computed[:, -1] = Az_computed[:, 0]
    # Az_computed[-1, :] = Az_computed[0, :]


    # Apply Gaussian smoothing with a standard deviation (sigma) of your choice
    sigma = sigma  # Adjust this value as needed
    Az_smooth = gaussian_filter(Az_computed, sigma=sigma)

    Az_computed = Az_smooth
    
    
    return grid_x, grid_y, Az_computed
                
def unsmooth_vect_pot(data, Ngrid_x = 2048, Ngrid_y = 2048, method='nearest',  x_range=None, y_range=None):
    """
    Interpolates the magnetic fields Bx and By, integrates them to obtain the vector potential Az 
    without applying any smoothing.

    Parameters:
    - data (dict): A dictionary containing the magnetic field data with keys 'b1' and 'b2'.
    - Ngrid_x (int): The number of grid points along the x-axis (default is 2048).
    - Ngrid_y (int): The number of grid points along the y-axis (default is 2048).
    - method (str): The interpolation method to use ('nearest', 'linear', or 'cubic'; default is 'nearest').

    Returns:
    - ndarray: The unsmoothed vector potential Az.
    """

    grid_x, grid_y, Bx_interp = interpolate_var_to_grid(data, "b1", Ngrid_x=Ngrid_x, Ngrid_y=Ngrid_y, method=method, x_range=x_range, y_range=y_range)
    _, _, By_interp = interpolate_var_to_grid(data, "b2", Ngrid_x=Ngrid_x, Ngrid_y=Ngrid_y, method=method, x_range=x_range, y_range=y_range)

    F = cumtrapz(Bx_interp, grid_y, axis=0, initial=0)        
    G = cumtrapz(-By_interp, grid_x, axis=1, initial=0) - F        
    Az_computed = F + G
    # Enforce periodic boundary conditions for Az_computed
    # Az_computed[:, -1] = Az_computed[:, 0]
    # Az_computed[-1, :] = Az_computed[0, :]
    
    return grid_x, grid_y, Az_computed

    
def add_arrow(line, position=None, direction='right', size=7, color=None):
    """
    Adds an arrow to a contour plot to indicate direction.

    Parameters:
    - line (Line2D): The line object to which the arrow will be added.
    - position (int, optional): The index of the position on the line to place the arrow (default is the middle).
    - direction (str): The direction of the arrow ('right' for forward, 'left' for backward; default is 'right').
    - size (int): The size of the arrow (default is 7).
    - color (str, optional): The color of the arrow (default is the color of the line).

    Returns:
    - None
    """

    if color is None:
        color = line.get_color()
    xdata = line.get_xdata()
    ydata = line.get_ydata()
    if position is None:
        position = xdata.size // 2
    if direction == 'right':
        dx = xdata[position + 1] - xdata[position]
        dy = ydata[position + 1] - ydata[position]
    else:
        dx = xdata[position - 1] - xdata[position]
        dy = ydata[position - 1] - ydata[position]
    line.axes.annotate('',
                       xytext=(xdata[position], ydata[position]),
                       xy=(xdata[position] + dx, ydata[position] + dy),
                       arrowprops=dict(arrowstyle="->", color=color),
                       size=size)
    
    
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
    
    

    

def format_sci_notation(value):
    """
    Formats a given numerical value into scientific notation for LaTeX.

    Parameters:
    - value (float): The numerical value to format.

    Returns:
    - str: A string representing the formatted value in scientific notation (e.g., "1.0 \times 10^{3}").
    """

    
    coeff, exp = f"{value:.1E}".split("E")
    exp = int(exp)
    if coeff == '1.0':
        return f"10^{{{exp}}}"
    else:
        return f"{coeff} \\times 10^{{{exp}}}"
    
    
def plot_blocks(data, fig=None, ax=None, linewidth=0.1, color='k', 
                          x_range=None, y_range=None):
    """
    Plot block boundaries efficiently using LineCollection with optional spatial range filtering.
    
    Parameters:
    - data: dict containing 'block_coord' numpy array of shape (N, 4, 2, 2)
    - fig: matplotlib figure object (optional)
    - ax: matplotlib axis object (optional)
    - linewidth: width of the boundary lines
    - color: color of the boundary lines
    - x_range: tuple (xmin, xmax) to limit plotted blocks within x bounds (optional)
    - y_range: tuple (ymin, ymax) to limit plotted blocks within y bounds (optional)
    """
    print('===============================')
    print("Started plotting block boundaries")
    start_time = time.time()

    block_coords = data['block_coord']
    N = block_coords.shape[0]  # Number of blocks

    # Prepare segments for LineCollection
    segments = []

    for i in range(N):
        block = block_coords[i]
        
        # Extract the x and y coordinates of the block
        x_block = block[:, 0, 0]  # Shape (4,)
        y_block = block[:, 0, 1]  # Shape (4,)
        
        # Check if the block falls within the provided spatial range (x_range, y_range)
        if x_range is not None:
            if np.all(x_block < x_range[0]) or np.all(x_block > x_range[1]):
                continue  # Skip this block if it lies entirely outside the x_range

        if y_range is not None:
            if np.all(y_block < y_range[0]) or np.all(y_block > y_range[1]):
                continue  # Skip this block if it lies entirely outside the y_range

        # Create segments for each side of the block
        segments.extend([
            [block[0, 0], block[1, 0]],  # Bottom edge
            [block[1, 0], block[2, 0]],  # Right edge
            [block[2, 0], block[3, 0]],  # Top edge
            [block[3, 0], block[0, 0]]   # Left edge
        ])

    # Create figure and axis if not provided
    if fig is None or ax is None:
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        ax.set_xlabel('$x/L$')
        ax.set_ylabel('$y/L$')

    # Create and add LineCollection
    lc = LineCollection(segments, linewidths=linewidth, colors=color)
    ax.add_collection(lc)

    # Set plot limits based on provided ranges or data bounds
    if x_range is None:
        ax.set_xlim(block_coords[:, :, 0, 0].min(), block_coords[:, :, 0, 0].max())
    else:
        ax.set_xlim(x_range)
    
    if y_range is None:
        ax.set_ylim(block_coords[:, :, 0, 1].min(), block_coords[:, :, 0, 1].max())
    else:
        ax.set_ylim(y_range)


    


    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Finished plotting block boundaries")
    print(f"Time taken: {elapsed_time:.4f} seconds")
    print('===============================')

    return fig, ax


def plot_cells(data, fig=None, ax=None, linewidth=0.1, color='k', 
                             x_range=None, y_range=None):
    """
    Optimized plotting of grid cells using LineCollection with optional spatial range filtering.
    
    Parameters:
    - data: dictionary containing 'xpoint', 'ypoint', 'ncells', 'offsets', 'connectivity'
    - fig: matplotlib figure object (optional)
    - ax: matplotlib axis object (optional)
    - linewidth: width of the boundary lines
    - color: color of the boundary lines
    - x_range: tuple (xmin, xmax) to limit plotted cells within x bounds (optional)
    - y_range: tuple (ymin, ymax) to limit plotted cells within y bounds (optional)
    """
    print('===============================')
    print("Started plotting grid cells")
    start_time = time.time()

    x = data['xpoint']
    y = data['ypoint']
    ncells = data['ncells']
    offsets = data['offsets']
    connectivity = data['connectivity']

    # Create mod_conn array using broadcasting
    base_conn = connectivity[:np.max(offsets)]
    num_iterations = int(4 * ncells / np.max(offsets))
    offsets_array = np.arange(num_iterations) * (np.max(base_conn) + 1)
    mod_conn = (base_conn + offsets_array[:, None]).ravel()[:ncells * 4]
    cell_vertices = mod_conn.reshape(ncells, 4)

    # Extract x and y coordinates for all cells at once
    x_vals = x[cell_vertices]
    y_vals = y[cell_vertices]

    # Apply spatial filtering based on the provided x_range and y_range
    if x_range is not None:
        x_mask = (x_vals.min(axis=1) >= x_range[0]) & (x_vals.max(axis=1) <= x_range[1])
    else:
        x_mask = np.ones(ncells, dtype=bool)

    if y_range is not None:
        y_mask = (y_vals.min(axis=1) >= y_range[0]) & (y_vals.max(axis=1) <= y_range[1])
    else:
        y_mask = np.ones(ncells, dtype=bool)

    # Combine masks to filter cells
    valid_cells = x_mask & y_mask
    x_vals = x_vals[valid_cells]
    y_vals = y_vals[valid_cells]
    filtered_ncells = len(x_vals)

    # Create segments for LineCollection
    segments = []
    for i in range(filtered_ncells):
        # Define the four corners of each cell
        corners = [
            (x_vals[i, 0], y_vals[i, 0]),
            (x_vals[i, 1], y_vals[i, 0]),
            (x_vals[i, 1], y_vals[i, 2]),
            (x_vals[i, 0], y_vals[i, 2]),
            (x_vals[i, 0], y_vals[i, 0])  # Close the loop
        ]
        # Add segments for each side of the cell
        segments.extend([
            [corners[0], corners[1]],
            [corners[1], corners[2]],
            [corners[2], corners[3]],
            [corners[3], corners[0]]
        ])

    # Create figure and axis if not provided
    if fig is None or ax is None:
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        ax.set_xlabel('$x/L$')
        ax.set_ylabel('$y/L$')

    # Create and add LineCollection
    lc = LineCollection(segments, linewidths=linewidth, colors=color)
    ax.add_collection(lc)

    # Set plot limits
    ax.set_xlim(x_range if x_range is not None else (x.min(), x.max()))
    ax.set_ylim(y_range if y_range is not None else (y.min(), y.max()))

    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Finished plotting grid cells")
    print(f"Time taken: {elapsed_time:.4f} seconds")
    print('===============================')

    return fig, ax



def plot_raw_data_cells(data, field_data, fig=None, ax=None, x_range=None, y_range=None, 
                        vmin=None, vmax=None, cmap='viridis', label=None, 
                        edgecolors=None, linewidths=0.1, orientation='vertical',  
                        location='right', use_log_norm=False, pad=0.1, colorbar=True):
    """
    Plots raw simulation data by filling each cell based on the field value associated with that cell.

    Parameters:
    - data (dict): A dictionary containing the following keys:
        - 'xpoint' (ndarray): 1D array of x-coordinates for cell vertices.
        - 'ypoint' (ndarray): 1D array of y-coordinates for cell vertices.
        - 'ncells' (int): The total number of cells to be plotted.
        - 'offsets' (ndarray): 1D array specifying the starting index of each cell's vertices in the connectivity array.
        - 'connectivity' (ndarray): 1D array that defines the connections between cell vertices.
    - field_data (ndarray): 1D array of field values corresponding to each cell. These values determine the color of each cell.
    - fig (matplotlib.figure.Figure, optional): A Matplotlib figure object. If not provided, a new figure is created.
    - ax (matplotlib.axes.Axes, optional): A Matplotlib axis object. If not provided, a new axis is created.
    - x_range (tuple, optional): A tuple (xmin, xmax) to limit plotted cells within specified x bounds. If None, no limits are applied.
    - y_range (tuple, optional): A tuple (ymin, ymax) to limit plotted cells within specified y bounds. If None, no limits are applied.
    - vmin (float, optional): Minimum value for color mapping. If None, the minimum value from filtered_field_data is used.
    - vmax (float, optional): Maximum value for color mapping. If None, the maximum value from filtered_field_data is used.
    - cmap (str, optional): Colormap used to color the cells. Default is 'viridis'.
    - label (str, optional): Label for the colorbar.
    - edgecolors (str, optional): Color of the edges of the cells. Default is None, which uses the default edge color.
    - linewidths (float, optional): Width of the cell boundaries. Default is 0.1.
    - orientation (str, optional): Orientation of the colorbar. Default is 'vertical'.
    - location (str, optional): Location of the colorbar. Default is 'right'.
    - use_log_norm (bool, optional): If True, applies logarithmic normalization to the field data for color mapping. Default is False.
    - pad: The colorbar padding. Default is 0.1.
    Returns:
    - fig (matplotlib.figure.Figure): The figure object containing the plot.
    - ax (matplotlib.axes.Axes): The axis object with the plotted data.
    """
    print('===============================')
    print("Started plotting raw data cells")
    start_time = time.time()

    x = data['xpoint']
    y = data['ypoint']
    ncells = data['ncells']
    offsets = data['offsets']
    connectivity = data['connectivity']

    # Create mod_conn array using broadcasting
    base_conn = connectivity[:np.max(offsets)]
    num_iterations = int(4 * ncells / np.max(offsets))
    offsets_array = np.arange(num_iterations) * (np.max(base_conn) + 1)
    mod_conn = (base_conn + offsets_array[:, None]).ravel()[:ncells * 4]
    cell_vertices = mod_conn.reshape(ncells, 4)

    # Extract x and y coordinates for all cells at once
    x_vals = x[cell_vertices]
    y_vals = y[cell_vertices]

    # Apply spatial filtering based on the provided x_range and y_range
    if x_range is not None:
        x_mask = (x_vals.min(axis=1) >= x_range[0]) & (x_vals.max(axis=1) <= x_range[1])
    else:
        x_mask = np.ones(ncells, dtype=bool)

    if y_range is not None:
        y_mask = (y_vals.min(axis=1) >= y_range[0]) & (y_vals.max(axis=1) <= y_range[1])
    else:
        y_mask = np.ones(ncells, dtype=bool)


    # Combine masks to filter cells
    valid_cells = x_mask & y_mask
    x_vals = x_vals[valid_cells]
    y_vals = y_vals[valid_cells]
    filtered_field_data = field_data[valid_cells]
    filtered_ncells = len(x_vals)

    # Create polygons for PolyCollection (one polygon per cell)
    polygons = []
    for i in range(filtered_ncells):
        polygon = [
            (x_vals[i, 0], y_vals[i, 0]),
            (x_vals[i, 1], y_vals[i, 0]),
            (x_vals[i, 1], y_vals[i, 2]),
            (x_vals[i, 0], y_vals[i, 2])
        ]
        polygons.append(polygon)

    # Create figure and axis if not provided
    if fig is None or ax is None:
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        ax.set_xlabel('$x/L$')
        ax.set_ylabel('$y/L$')

    if vmin == None:
        vmin = np.min(filtered_field_data)
        
    if vmax == None:
        vmax = np.max(filtered_field_data)
        
    if use_log_norm:
        # Create a log normalization instance
        norm = LogNorm(vmin=vmin, 
                   vmax=vmax)
        poly_collection = PolyCollection(polygons, norm=norm, array=filtered_field_data, 
                                         cmap=cmap, edgecolors=edgecolors, 
                                         clim=(vmin,vmax), linewidths=linewidths)

    else:
        # Create a PolyCollection with the polygons and color them by field data
        poly_collection = PolyCollection(polygons, array=filtered_field_data, 
                                         cmap=cmap, edgecolors=edgecolors, 
                                         clim=(vmin,vmax), linewidths=linewidths)
    ax.add_collection(poly_collection)
    


    # Add colorbar with 'extend' parameter determined from the data
    if colorbar:
        # Determine extend based on comparisons
        extend_type = 'neither'  # Default
        if np.any(filtered_field_data < vmin):
            extend_type = 'min'
        if np.any(filtered_field_data > vmax):
            extend_type = 'max'
        if np.any(filtered_field_data < vmin) and np.any(filtered_field_data > vmax):
            extend_type = 'both'
        cbar = plt.colorbar(poly_collection, ax=ax, extend=extend_type, 
                            label=label, orientation=orientation, location=location, 
                            pad=pad)

    # Set plot limits
    # ax.set_xlim(x_range if x_range is not None else (x.min(), x.max()))
    # ax.set_ylim(y_range if y_range is not None else (y.min(), y.max()))
    ax.set_xlim(x_vals.min(), x_vals.max())
    ax.set_ylim(y_vals.min(), y_vals.max())


    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Finished plotting raw data cells")
    print(f"Time taken: {elapsed_time:.4f} seconds")
    print('===============================')

    return poly_collection, fig, ax



def extract_1d_slice(data, axis, slice_value, vars, tol=1e-5):
    """
    Extract a 1D slice of the data from a 2D grid, along a specified axis.

    Args:
    - data (dict): A dictionary containing the data, including 'center_x', 'center_y', 
    and the variable to interpolate.
    - axis (int): the axis to slice along (0 for x, 1 for y).
    - slice_value (float): the value of the axis where the slice should be taken.
    - vars (dict): any number of cell-centered data arrays names to slice eg {"b1", "b2"}.
    - tol: tolerance of finding cells near slice_value along slice. 
        
    Returns:
    - tuple: tuple containing the coordinates along the slice and the sliced variables 
    as a list of lists of more than one variable is sliced
    """
    
    print('===============================')
    print("Started 1d slice")
    start_time = time.time()

    
    center_x, center_y = data["center_x"], data["center_y"]

    # Choose the axis for slicing (0 = x-axis, 1 = y-axis)
    if axis == 0:  # Slice along x (i.e., fixed x = slice_value)
        slice_indices = np.isclose(center_x, slice_value, atol=tol)
        slice_coord = center_y[slice_indices]
    elif axis == 1:  # Slice along y (i.e., fixed y = slice_value)
        slice_indices = np.isclose(center_y, slice_value, atol=tol)
        slice_coord = center_x[slice_indices]
    else:
        raise ValueError("Invalid axis. Use 0 for x and 1 for y.")
    
    if not np.any(slice_indices):
        raise ValueError(f"No data found for slice at {slice_value} along axis {axis}.")
    
    # Extract the slice from each variable
    sliced_vars = [data[var][slice_indices] for var in vars]
    
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Finished 1d slice")
    print(f"Time taken: {elapsed_time:.4f} seconds")
    print('===============================')

    return slice_coord, sliced_vars



def plot_polar_data_cells_continuous(data, field_data, fig=None, ax=None, x_range=None,
                                     y_range=None, vmin=None, vmax=None, cmap='viridis', label=None, 
                                     orientation='vertical', location='right', use_log_norm=False, 
                                     pad=0.1, colorbar=True, resolution=512):
    """
    Plots polar simulation data on a Cartesian grid as a continuous field, resulting in a half-circle shape.
    Interpolates only over the specified x and y ranges.

    Parameters:
    - data (dict): A dictionary containing the following keys:
        - 'xpoint' (ndarray): 1D array of x-coordinates for cell vertices.
        - 'ypoint' (ndarray): 1D array of y-coordinates for cell vertices.
        - 'ncells' (int): The total number of cells to be plotted.
        - 'offsets' (ndarray): 1D array specifying the starting index of each cell's vertices in the connectivity array.
        - 'connectivity' (ndarray): 1D array that defines the connections between cell vertices.
    - field_data (ndarray): 1D array of field values corresponding to each cell. These values determine the color of each cell.
    - fig (matplotlib.figure.Figure, optional): A Matplotlib figure object. If not provided, a new figure is created.
    - ax (matplotlib.axes.Axes, optional): A Matplotlib axis object. If not provided, a new axis is created.
    - x_range (tuple, optional): A tuple (xmin, xmax) to limit plotted cells within specified x bounds. If None, no limits are applied.
    - y_range (tuple, optional): A tuple (ymin, ymax) to limit plotted cells within specified y bounds. If None, no limits are applied.
    - vmin (float, optional): Minimum value for color mapping. If None, the minimum value from filtered_field_data is used.
    - vmax (float, optional): Maximum value for color mapping. If None, the maximum value from filtered_field_data is used.
    - cmap (str, optional): Colormap used to color the cells. Default is 'viridis'.
    - label (str, optional): Label for the colorbar.
    - orientation (str, optional): Orientation of the colorbar. Default is 'vertical'.
    - location (str, optional): Location of the colorbar. Default is 'right'.
    - use_log_norm (bool, optional): If True, applies logarithmic normalization to the field data for color mapping. Default is False.
    - pad: The colorbar padding. Default is 0.1.

    Returns:
    - im: The image object from pcolormesh.
    - fig, ax: Matplotlib figure and axis objects.
    """
    print('===============================')
    print("Started plotting continuous polar data")
    start_time = time.time()

    # Extract data
    center_x = data['center_x']
    center_y = data['center_y']
        
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
        # y_mask = ((center_y >= y_range[0][0]) & (center_y <= y_range[0][1])) | \
        #          ((center_y >= y_range[1][0]) & (center_y <= y_range[1][1]))
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

    # Filter data points
    filtered_x = center_x[mask]
    filtered_y = center_y[mask]
    filtered_field_data = field_data[mask]
    
    xi = np.linspace(x_min, x_max, resolution)
    yi = np.linspace(y_min, y_max, resolution)
    XI, YI = np.meshgrid(xi, yi)


    # Interpolate the data onto the grid
    points = np.column_stack((filtered_x, filtered_y))
    grid_z = griddata(points, filtered_field_data, (XI, YI), method='nearest')

    # Create a mask for the half-circle (if applicable)
    mask1 = XI**2 + YI**2 <= rmax**2
    mask2 = XI**2 + YI**2 > rmin**2
    
    mask = mask1 & mask2

    grid_z = np.ma.masked_where(~mask, grid_z)

    # Create figure and axis if not provided
    if fig is None or ax is None:
        fig, ax = plt.subplots()
    ax.set_xlabel('$x/R_{\\rm{NS}}$')
    ax.set_ylabel('$y/R_{\\rm{NS}}$')
    ax.set_aspect('equal')

    # Set vmin and vmax if not provided
    if vmin is None:
        vmin = np.nanmin(grid_z)
    if vmax is None:
        vmax = np.nanmax(grid_z)

    # Apply log normalization if requested
    if use_log_norm:
        norm = LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = None

    # Plot the interpolated data
    im = ax.pcolormesh(XI, YI, grid_z, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)

    # Add colorbar
    if colorbar:
        extend_type = 'neither'
        if np.any(grid_z < vmin):
            extend_type = 'min'
        if np.any(grid_z > vmax):
            extend_type = 'max'
        if np.any(grid_z < vmin) and np.any(grid_z > vmax):
            extend_type = 'both'
        
        cbar = plt.colorbar(im, ax=ax, extend=extend_type, label=label, 
                            orientation=orientation, location=location, pad=pad)

    # Set plot limits
    ax.set_xlim(x_range)
    ax.set_ylim(y_range)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Finished plotting continuous polar data")
    print(f"Time taken: {elapsed_time:.4f} seconds")
    print('===============================')

    return im, fig, ax







def fast_vtu_reader_ascii(filename, attr='all', blocks=False, slicedir='z'):
    """
    Reads a VTU file produced by BHAC for 2D simulations, assuming the data is in ASCII format.
    Need to tell it direction of slice, assumes along z.
    It will always output cell centers x and y but that is just so it will work with the 
    plotting function, in reality center_x and center_y are the coordinates not sliced along,
    eg if slice along x then center_x is y and center_y is z.
    This is only for if the vtu is not cell centered.
    
    Parameters:
    - filename: str, path to the VTU file.
    - attr: list or 'all', attributes to extract (default is 'all').
    - blocks: bool, whether to compute block boundaries for visualization (default is False).
    - blocks: str, tell the direction of slice (default is 'z').

    Returns:
    - data: dict, containing extracted points, attributes, and calculated cell centers.
    - data_keys: list, names of the data arrays present in the file.
    """
    
    print('===============================')
    print(f"Starting to read file: {filename}")
    start_time = time.time()

    with open(filename, 'r') as f:
        content = f.read()

    # Parse the XML part of the file
    xml_content_start = content.find('<VTKFile')
    xml_content_end = content.find('</VTKFile>') + len('</VTKFile>')
    xml_content = content[xml_content_start:xml_content_end]
    
    root = ET.fromstring(xml_content)
    
    data = {}
    
    time_field = root.find(".//FieldData/DataArray[@Name='TIME']")
    if time_field is not None:
        data['time'] = float(time_field.text.strip())
        print(f"Extracted time: {data['time']}")
    else:
        print("No time field found in the file.")

    pieces = root.findall('.//Piece')
    num_pieces = len(pieces)
    print(f"Number of Pieces: {num_pieces}")

    cells_per_piece = int(pieces[0].get('NumberOfCells'))
    total_cells = cells_per_piece * num_pieces
    print(f"Cells per piece: {cells_per_piece}")
    print(f"Total number of cells: {total_cells}")

    # Get all unique DataArray names
    data_array_names = set()
    for piece in pieces:
        for data_array in piece.findall('.//DataArray'):
            data_array_names.add(data_array.get('Name'))

    # Read Points (x, y, z coordinates)
    points_data = []
    block_boundaries = []
    for piece in pieces:
        points_data_array = piece.find('.//Points/DataArray')
        if points_data_array is None:
            raise ValueError("Points data not found")

        dtype = points_data_array.get('type')
        format = points_data_array.get('format')

        # Read ASCII data
        raw_data = points_data_array.text.strip().split()
        
        if dtype == 'Float32':
            parsed_data = np.array(raw_data, dtype=np.float32)
        elif dtype == 'Float64':
            parsed_data = np.array(raw_data, dtype=np.float64)
        else:
            raise ValueError(f"Unsupported data type for Points: {dtype}")
        
        points_data.append(parsed_data)
        
        if blocks:
            # Reshape the parsed data for this piece
            piece_points = parsed_data.reshape(-1, 3)

            # Vectorized min and max operations for this piece
            x_min, y_min = np.min(piece_points[:, :2], axis=0)
            x_max, y_max = np.max(piece_points[:, :2], axis=0)

            # Define block corners for this piece
            corners = np.array([
                [x_min, y_min],
                [x_max, y_min],
                [x_max, y_max],
                [x_min, y_max],
                [x_min, y_min]  # Close the loop
            ])

            # Create block boundaries for this piece
            piece_boundaries = np.array([corners[:-1], corners[1:]]).transpose(1, 0, 2)
            
            block_boundaries.append(piece_boundaries)

    if blocks:
        data['block_coord'] = np.array(block_boundaries)

    if points_data:
        points = np.concatenate(points_data).reshape(-1, 3)  # Assuming 3D points (x, y, z)
        data['xpoint'], data['ypoint'], data['zpoint'] = points[:, 0], points[:, 1], points[:, 2]
        print(f"Extracted {len(data['xpoint'])} points")

    # Handle attributes
    if attr == 'all':
        data_array_names.discard(None)
        data_array_names.discard('types')
    else:
        data_array_names = attr
        data_array_names.add('connectivity')
        data_array_names.add('offsets')

    for name in data_array_names:
        combined_data = []

        for piece in pieces:
            piece_data_array = piece.find(f".//DataArray[@Name='{name}']")
            if piece_data_array is None:
                continue

            dtype = piece_data_array.get('type')
            format = piece_data_array.get('format')

            # Read ASCII data
            raw_data = piece_data_array.text.strip().split()

            if dtype == 'Float32':
                parsed_data = np.array(raw_data, dtype=np.float32)
            elif dtype == 'Float64':
                parsed_data = np.array(raw_data, dtype=np.float64)
            elif dtype == 'Int32':
                parsed_data = np.array(raw_data, dtype=np.int32)
            elif dtype == 'Int64':
                parsed_data = np.array(raw_data, dtype=np.int64)
            else:
                raise ValueError(f"Unsupported data type: {dtype}")

            combined_data.append(parsed_data)

        if combined_data:
            data[name] = np.concatenate(combined_data)

    data["ncells"] = total_cells
    
    if slicedir == 'z':
        x,y,z = data['xpoint'], data['ypoint'], data['zpoint']
        data["center_x"], data["center_y"] = calculate_cell_centers(data)
    
    if slicedir == 'y':
        x,y,z = data['xpoint'], data['ypoint'], data['zpoint']
        data['ypoint'] = z
        data["center_x"], data["center_y"] = calculate_cell_centers(data)
    
    if slicedir == 'x':
        x,y,z = data['xpoint'], data['ypoint'], data['zpoint']
        data['xpoint'] = y
        data['ypoint'] = z
        data["center_x"], data["center_y"] = calculate_cell_centers(data)
    
    
    
    # List of keys to exclude
    exclude_keys = ['time', 'xpoint', 'ypoint', 'zpoint', 'offsets', 'connectivity',
                     'ncells', 'center_x', 'center_y', 'center_z', 'block_coord']

    # Create a new dictionary with the excluded keys removed
    data_dict = {key: value for key, value in data.items() if key not in exclude_keys}

    # For some reason the field data is on the cell corners, need to move it the centers
    ncells = data['ncells']
    
    offsets = data['offsets']
    connectivity = data['connectivity']
    # Create mod_conn array using broadcasting instead of a for loop
    base_conn = connectivity[:np.max(offsets)]  # Base mod_conn for the first set
    num_iterations = int(4 * ncells / np.max(offsets))  # Number of iterations

    # Use broadcasting to create mod_conn without a loop
    offsets_array = np.arange(num_iterations) * (np.max(base_conn) + 1)  # Calculate all offsets at once
    mod_conn = base_conn + offsets_array[:, None]  # Broadcast and add offsets
    
    # Flatten mod_conn to a 1D array
    mod_conn = mod_conn.ravel()[:ncells * 4]  # Only take enough entries for ncells

    # Reshape mod_conn to group cell vertices (ncells x 4)
    cell_vertices = mod_conn.reshape(ncells, 4)

    for key in data_dict:
        data[key] = np.mean(data[key][cell_vertices], axis=1)
        
    # # Vectorized calculation of cell centers
    # cell_centers_b1 = np.mean(b1[cell_vertices], axis=1)

    

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Finished reading file: {filename}")
    print(f"Time taken to read: {elapsed_time:.4f} seconds")
    print('===============================')

    return data, list(data.keys())





def fast_vtu_reader_ascii_cc(filename, attr='all', blocks=False, slicedir='z'):
    """
    Reads a VTU file produced by BHAC for 2D simulations, assuming the data is in ASCII format.
    Need to tell it direction of slice, assumes along z.
    It will always output cell centers x and y but that is just so it will work with the 
    plotting function, in reality center_x and center_y are the coordinates not sliced along,
    eg if slice along x then center_x is y and center_y is z.
    This is only for if the vtu is cell centered.
    
    Parameters:
    - filename: str, path to the VTU file.
    - attr: list or 'all', attributes to extract (default is 'all').
    - blocks: bool, whether to compute block boundaries for visualization (default is False).
    - blocks: str, tell the direction of slice (default is 'z').

    Returns:
    - data: dict, containing extracted points, attributes, and calculated cell centers.
    - data_keys: list, names of the data arrays present in the file.
    """
    
    print('===============================')
    print(f"Starting to read file: {filename}")
    start_time = time.time()

    with open(filename, 'r') as f:
        content = f.read()

    # Parse the XML part of the file
    xml_content_start = content.find('<VTKFile')
    xml_content_end = content.find('</VTKFile>') + len('</VTKFile>')
    xml_content = content[xml_content_start:xml_content_end]
    
    root = ET.fromstring(xml_content)
    
    data = {}
    
    time_field = root.find(".//FieldData/DataArray[@Name='TIME']")
    if time_field is not None:
        data['time'] = float(time_field.text.strip())
        print(f"Extracted time: {data['time']}")
    else:
        print("No time field found in the file.")

    pieces = root.findall('.//Piece')
    num_pieces = len(pieces)
    print(f"Number of Pieces: {num_pieces}")

    cells_per_piece = int(pieces[0].get('NumberOfCells'))
    total_cells = cells_per_piece * num_pieces
    print(f"Cells per piece: {cells_per_piece}")
    print(f"Total number of cells: {total_cells}")

    # Get all unique DataArray names
    data_array_names = set()
    for piece in pieces:
        for data_array in piece.findall('.//DataArray'):
            data_array_names.add(data_array.get('Name'))

    # Read Points (x, y, z coordinates)
    points_data = []
    block_boundaries = []
    for piece in pieces:
        points_data_array = piece.find('.//Points/DataArray')
        if points_data_array is None:
            raise ValueError("Points data not found")

        dtype = points_data_array.get('type')
        format = points_data_array.get('format')

        # Read ASCII data
        raw_data = points_data_array.text.strip().split()
        
        if dtype == 'Float32':
            parsed_data = np.array(raw_data, dtype=np.float32)
        elif dtype == 'Float64':
            parsed_data = np.array(raw_data, dtype=np.float64)
        else:
            raise ValueError(f"Unsupported data type for Points: {dtype}")
        
        points_data.append(parsed_data)
        
        if blocks:
            # Reshape the parsed data for this piece
            piece_points = parsed_data.reshape(-1, 3)

            # Vectorized min and max operations for this piece
            x_min, y_min = np.min(piece_points[:, :2], axis=0)
            x_max, y_max = np.max(piece_points[:, :2], axis=0)

            # Define block corners for this piece
            corners = np.array([
                [x_min, y_min],
                [x_max, y_min],
                [x_max, y_max],
                [x_min, y_max],
                [x_min, y_min]  # Close the loop
            ])

            # Create block boundaries for this piece
            piece_boundaries = np.array([corners[:-1], corners[1:]]).transpose(1, 0, 2)
            
            block_boundaries.append(piece_boundaries)

    if blocks:
        data['block_coord'] = np.array(block_boundaries)

    if points_data:
        points = np.concatenate(points_data).reshape(-1, 3)  # Assuming 3D points (x, y, z)
        data['xpoint'], data['ypoint'], data['zpoint'] = points[:, 0], points[:, 1], points[:, 2]
        print(f"Extracted {len(data['xpoint'])} points")

    # Handle attributes
    if attr == 'all':
        data_array_names.discard(None)
        data_array_names.discard('types')
    else:
        data_array_names = attr
        data_array_names.add('connectivity')
        data_array_names.add('offsets')

    for name in data_array_names:
        combined_data = []

        for piece in pieces:
            piece_data_array = piece.find(f".//DataArray[@Name='{name}']")
            if piece_data_array is None:
                continue

            dtype = piece_data_array.get('type')
            format = piece_data_array.get('format')

            # Read ASCII data
            raw_data = piece_data_array.text.strip().split()

            if dtype == 'Float32':
                parsed_data = np.array(raw_data, dtype=np.float32)
            elif dtype == 'Float64':
                parsed_data = np.array(raw_data, dtype=np.float64)
            elif dtype == 'Int32':
                parsed_data = np.array(raw_data, dtype=np.int32)
            elif dtype == 'Int64':
                parsed_data = np.array(raw_data, dtype=np.int64)
            else:
                raise ValueError(f"Unsupported data type: {dtype}")

            combined_data.append(parsed_data)

        if combined_data:
            data[name] = np.concatenate(combined_data)

    data["ncells"] = total_cells
    
    if slicedir == 'z':
        x,y,z = data['xpoint'], data['ypoint'], data['zpoint']
        data["center_x"], data["center_y"] = calculate_cell_centers(data)
    
    if slicedir == 'y':
        x,y,z = data['xpoint'], data['ypoint'], data['zpoint']
        data['ypoint'] = z
        data["center_x"], data["center_y"] = calculate_cell_centers(data)
    
    if slicedir == 'x':
        x,y,z = data['xpoint'], data['ypoint'], data['zpoint']
        data['xpoint'] = y
        data['ypoint'] = z
        data["center_x"], data["center_y"] = calculate_cell_centers(data)
    

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Finished reading file: {filename}")
    print(f"Time taken to read: {elapsed_time:.4f} seconds")
    print('===============================')

    return data, list(data.keys())




def calculate_cell_centers_3d(data):
    """
    Calculate the cell centers by averaging the vertex coordinates for 3D data.

    Parameters:
    - data: dict, containing vertex coordinates and connectivity information.

    Returns:
    - center_x: ndarray, x-coordinates of cell centers.
    - center_y: ndarray, y-coordinates of cell centers.
    - center_z: ndarray, z-coordinates of cell centers.
    """

    print('===============================')
    print(f"Started finding cell centers (3D)")
    start_time = time.time()
    
    x = data['xpoint']
    y = data['ypoint']
    z = data['zpoint']
    ncells = data['ncells']
    
    offsets = data['offsets']
    connectivity = data['connectivity']
    
    # Create mod_conn array using broadcasting instead of a for loop
    base_conn = connectivity[:np.max(offsets)]  # Base mod_conn for the first set
    num_iterations = int(8 * ncells / np.max(offsets))  # Number of iterations for 3D (8 vertices per cell)

    # Use broadcasting to create mod_conn without a loop
    offsets_array = np.arange(num_iterations) * (np.max(base_conn) + 1)  # Calculate all offsets at once
    mod_conn = base_conn + offsets_array[:, None]  # Broadcast and add offsets
    
    # Flatten mod_conn to a 1D array
    mod_conn = mod_conn.ravel()[:ncells * 8]  # Only take enough entries for ncells (8 vertices per cell)

    # Reshape mod_conn to group cell vertices (ncells x 8)
    cell_vertices = mod_conn.reshape(ncells, 8)

    # Vectorized calculation of cell centers
    cell_centers_x = np.mean(x[cell_vertices], axis=1)
    cell_centers_y = np.mean(y[cell_vertices], axis=1)
    cell_centers_z = np.mean(z[cell_vertices], axis=1)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Finished finding cell centers (3D)")
    print(f"Time taken to get centers: {elapsed_time:.4f} seconds")
    print('===============================')

    return cell_centers_x, cell_centers_y, cell_centers_z





def interp_2d_slice(data, var, slice_direction='xy', slice_position=0, atol=5e-1, 
                    Ngrid_x=256, Ngrid_y=256, method='nearest', x_range=None, y_range=None):
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
    print('===============================')
    print(f"Started interpolating")
    start_time = time.time()

    x, y, z = calculate_cell_centers_3d(data)
    coords = {'x': x, 'y': y, 'z': z}

    # Define slice based on slice_direction
    if slice_direction == 'xy':
        mask = np.isclose(coords['z'], slice_position, atol=atol)
        x, y = coords['x'][mask], coords['y'][mask]
    elif slice_direction == 'yz':
        mask = np.isclose(coords['x'], slice_position, atol=atol)
        x, y = coords['y'][mask], coords['z'][mask]
    elif slice_direction == 'xz':
        mask = np.isclose(coords['y'], slice_position, atol=atol)
        x, y = coords['x'][mask], coords['z'][mask]
    else:
        raise ValueError("Invalid slice_direction. Use 'xy', 'yz', or 'xz'.")
    
    center_x, center_y = x, y
    
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
    filtered_var_data = data[var][mask]  # Ensure var data is filtered according to the same mask

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


    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Finished interpolating")
    print(f"Time taken to interpolate: {elapsed_time:.4f} seconds")
    print('===============================')

    return grid_x, grid_y, interpolated_var



def plot_interpolated_2d_slice(data, field_data, slice_direction='xy', slice_position=0, 
                               resolution=256, fig=None, ax=None, vmin=None, vmax=None, 
                               cmap='viridis', label=None, colorbar=True, use_log_norm=False, 
                               atol=1e-2, method='nearest',
                               x_range=None, y_range=None,
                               orientation='vertical', location='right',
                               pad=0.1):
    """
    Plots the interpolated data from a 2D slice of a 3D simulation.

    Parameters:
    - data (dict): Contains 3D coordinate arrays ('x', 'y', 'z') of the simulation data.
    - field_data (ndarray): 1D array of the field values for interpolation.
    - slice_direction (str): Direction to slice ('xy', 'yz', 'xz').
    - slice_position (float): Position along the orthogonal axis for slicing.
    - resolution (int): Resolution of the interpolated grid.
    - fig, ax (optional): Matplotlib figure and axis objects.
    - vmin, vmax (optional): Min/Max values for color mapping.
    - cmap (str, optional): Colormap used to color the cells.
    - label (str, optional): Label for the colorbar.
    - colorbar (bool, optional): Whether to display the colorbar. Default is True.
    - use_log_norm (bool, optional): If True, apply log normalization.

    Returns:
    - im: Image object from pcolormesh.
    - fig, ax: Matplotlib figure and axis objects.
    """
    print('===============================')
    print(f"Plotting a 2D slice in the {slice_direction}-plane at {slice_position}")
    
    x, y, z = calculate_cell_centers_3d(data)
    coords = {'x': x, 'y': y, 'z': z}

    # Define slice based on slice_direction
    if slice_direction == 'xy':
        mask = np.isclose(coords['z'], slice_position, atol=atol)
        x, y = coords['x'][mask], coords['y'][mask]
    elif slice_direction == 'yz':
        mask = np.isclose(coords['x'], slice_position, atol=atol)
        x, y = coords['y'][mask], coords['z'][mask]
    elif slice_direction == 'xz':
        mask = np.isclose(coords['y'], slice_position, atol=atol)
        x, y = coords['x'][mask], coords['z'][mask]
    else:
        raise ValueError("Invalid slice_direction. Use 'xy', 'yz', or 'xz'.")
    
        
    radii = np.sqrt(x**2 + y**2)
    rmin = np.min(radii)
    rmax = np.max(radii)
    

    # Calculate the wedge angle from the data
    angles = np.arctan2(y, x)
    wedge_min_angle = np.min(angles)
    wedge_max_angle = np.max(angles)
    # print(f"Calculated wedge angle: [{wedge_min_angle}, {wedge_max_angle}] radians")

    # Create masks for x_range, y_range, and rmin
    if x_range is not None:
        x_min = x_range[0]
        x_max = x_range[1]
        x_mask = (x >= x_range[0]) & (x <= x_range[1])
    else:
        if np.min(x) < 0:
            x_min = np.min(x) #-rmax
        else:
            x_min = 0
        if np.max(x) >= 0:
            x_max = np.max(x) #rmax
        else:
            x_max = 0
        x_mask = np.ones_like(x, dtype=bool)

    if y_range is not None:
        y_min = y_range[0]
        y_max = y_range[1]
        y_mask = (y >= y_range[0]) & (y <= y_range[1])
    else:
        if np.min(y) < 0:
            y_min = np.min(y) #-rmax
        else:
            y_min = 0
        if np.max(y) >= 0:
            y_max = np.max(y) #rmax
        else:
            y_max = 0
        y_mask = np.ones_like(y, dtype=bool)

    field_slice = field_data[mask]

    r_mask = radii >= rmin

    # Combine all masks
    mask = x_mask & y_mask & r_mask

    # Filter data points
    filtered_x = x[mask]
    filtered_y = y[mask]
    filtered_field_data = field_slice[mask]

    xi = np.linspace(x_min, x_max, resolution)
    yi = np.linspace(y_min, y_max, resolution)
    XI, YI = np.meshgrid(xi, yi)

    points = np.column_stack((filtered_x, filtered_y))
    grid_z = griddata(points, filtered_field_data, (XI, YI), method=method)
    
    # Create a mask for the wedge-shaped region
    angles_grid = np.arctan2(YI, XI)
    mask1 = XI**2 + YI**2 <= rmax**2
    mask2 = XI**2 + YI**2 > rmin**2
    mask3 = (angles_grid >= wedge_min_angle) & (angles_grid <= wedge_max_angle)

    mask = mask1 & mask2 & mask3

    grid_z = np.ma.masked_where(~mask, grid_z)

    # Create figure and axis if not provided
    if fig is None or ax is None:
        fig, ax = plt.subplots()

    # Apply log normalization if requested
    if use_log_norm:
        norm = LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = None

    # Set vmin and vmax if not provided
    if vmin is None:
        vmin = np.nanmin(grid_z)
    if vmax is None:
        vmax = np.nanmax(grid_z)

    im = ax.pcolormesh(XI, YI, grid_z, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)

    # Add colorbar
    if colorbar:
        extend_type = 'neither'
        if np.any(grid_z < vmin):
            extend_type = 'min'
        if np.any(grid_z > vmax):
            extend_type = 'max'
        if np.any(grid_z < vmin) and np.any(grid_z > vmax):
            extend_type = 'both'

        cbar = plt.colorbar(im, ax=ax, extend=extend_type, label=label, 
                            orientation=orientation, location=location, pad=pad)

    # Set plot limits
    ax.set_xlim(x_range)
    ax.set_ylim(y_range)

    # Set labels
    ax.set_xlabel(f'${slice_direction[0]}$''$/R_{\\rm{NS}}$')
    ax.set_ylabel(f'${slice_direction[1]}$''$/R_{\\rm{NS}}$')
    ax.set_aspect('equal')

    print(f"Finished plotting slice in {slice_direction}-plane.")
    print('===============================')

    return im, fig, ax
