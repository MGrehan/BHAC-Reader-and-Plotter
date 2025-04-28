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
- Plotting of raw simulation data
- 1d slicing of raw simulation data


--------------
Usage Example (interpolation):
--------------
```
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
cbar = fig.colorbar(p1, ax=ax, pad=0.05,  extend=determine_extend_from_plot(p1), orientation='horizontal',  location='top')
cbar.set_label('$\\beta$')
plot_cells(data, fig=fig, ax=ax, linewidth=0.25, color='w', x_range=(xmin,xmax), y_range=(ymin,ymax))
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin ,ymax)
plt.show()
```
--------------
Usage Example (raw data plotting):
--------------
```
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
```


--------------
Usage Example (polar plotting):
--------------
```
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
```


--------------
Usage Example (plotting sliced data from 3D simulation):
--------------
```
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
```

--------------
Usage Example (plotting sliced 3D polar data in 2D with slicing):
--------------
```
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
```

--------------
Usage Example (plotting using fast plotting):
--------------
```
fig, ax = plt.subplots()
xmin, xmax = -0.5, 0.5
ymin, ymax = -0.5, 0.5

plot_cell_centers_fast(data, data['rho'], fig=fig, ax=ax, 
                       x_range=(xmin,xmax), y_range=(ymin,ymax), cmap='cmr.iceburn', 
                       label='$\\rho$', orientation='vertical',  location='right', 
                       use_log_norm=False, pad=0.0, workers=-1, scaling=1)
# plot_blocks(data, fig=fig, ax=ax, color='w')
ax.set_xlabel('$x/L$')
ax.set_ylabel('$y/L$')
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin ,ymax)
ax.set_title(f'$t/t_c = {t:.1f}$')
plt.show()
```


--------------
Usage Example (polar plotting with flux function):
--------------
```
vmin, vmax = -0.5, 0.5
xmin, xmax = 0.1, 5.0 
ymin, ymax = -2.5, 2.5 
figsize=(11,8)
fig, ax = plt.subplots(figsize=figsize)
grid_x, grid_z, Psi = smooth_flux_fcn(data, Ngrid_x = 4096, Ngrid_y = 4096, sigma=5, method='nearest',  x_range=(xmin,xmax), y_range=(ymin,ymax))
contour = ax.contour(grid_x, grid_z, Psi, levels=35, colors='w', linewidths=1.0, linestyles='-')
        
p1, _, _ = plot_polar_data_cells_continuous(data, 
                            (r**(3/2)) * bphi, 
                            fig=fig, ax=ax, label='$(r/R_{\\rm{NS}})^3 B^{\\theta}/B_\\star$', 
                            x_range=(xmin,xmax), y_range=(ymin, ymax), 
                            resolution = 4096,
                            colorbar=None, cmap=cmr.iceburn,
                            vmin = vmin, vmax = vmax)


cbar = fig.colorbar(p1, ax=ax, pad=0.01,  
                    extend=determine_extend_from_plot(p1), label='$(r/R_{\\rm{NS}})^{3/2} B^{\\phi}/B_\\star$')
ax.tick_params(axis='both', which='both', length=0)

rescale_axis_labels(ax, xscale=1/RLC, yscale=1/RLC)
ax.set_xlabel('$x/R_{\\rm{LC}}$')
ax.set_ylabel('$z/R_{\\rm{LC}}$')
ax.set_title(f'$t/P = {t*Omega/(2*np.pi):.3f}$')  
```

--------------
Usage Example (polar plotting with flux function and parallelized functions):
--------------

```
vmin, vmax = -0.25, 0.25
xmin, xmax = 0.1, 5.0 
ymin, ymax = -2.5, 2.5 
figsize=(11,8)

fig, ax = plt.subplots(figsize=figsize)

grid_x, grid_z, Psi = smooth_flux_fcn_fast(data, Ngrid_x = 2048, Ngrid_y = 2048, sigma=5, method='nearest',  x_range=(xmin,xmax), y_range=(ymin,ymax), workers=-1)
contour = ax.contour(grid_x, grid_z, Psi, levels=35, colors='w', linewidths=1.0, linestyles='-')

p1, _, _ = plot_polar_cell_centers_fast(data, (r**(3/2)) * bphi, fig=fig, ax=ax,
                              x_range=(xmin,xmax), y_range=(ymin,ymax), vmin=vmin, vmax=vmax,
                              cmap=cmr.iceburn, label=None, orientation='vertical',
                              location='right', use_log_norm=False, pad=0.1, workers=-1, 
                              resolution = 2048, colorbar=False)

cbar = fig.colorbar(p1, ax=ax, pad=0.01,  
                    extend=determine_extend_from_plot(p1), label='$(r/R_{\\rm{NS}})^{3/2} B^{\\phi}/B_\\star$')
ax.tick_params(axis='both', which='both', length=0)

rescale_axis_labels(ax, xscale=1/RLC, yscale=1/RLC)
ax.set_xlabel('$x/R_{\\rm{LC}}$')
ax.set_ylabel('$z/R_{\\rm{LC}}$')
ax.set_title(f'$t/P = {t*Omega/(2*np.pi):.3f}$')  
```

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
