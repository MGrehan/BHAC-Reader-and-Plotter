BHAC VTU Reader for 2D Data (No General Relativity)

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
Az_computed = smooth_vect_pot(data, Ngrid_x = Ngrid_x, Ngrid_y = Ngrid_y)
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
data, names = fast_vtu_reader(filename, attr={'p', 'b1', 'b2', 'b3'}, blocks=Fale)
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
Libraries Used:
--------------
- struct: For handling binary data.
- numpy: For efficient numerical operations.
- xml.etree.ElementTree: For parsing the VTU XML structure.
- time: For timing the file reading and processing steps.
- base64: For decoding base64-encoded data.
- scipy.integrate, scipy.ndimage, scipy.interpolate: For interpolation and smoothing.
- matplotlib.pyplot: For plotting the data.
