"""
BHAC VTU aux functions

Author: Michael Patrick Grehan
Date: January 2025
Email: michael.grehan@mail.utoronto.ca

------------------------------------------------------------------------------------
This module provides functionality for processing BHAC VTU data files
that do not include GR effects. The module is designed to read 
VTU files, extract relevant data fields, and return them in a format suitable for 
numerical analysis and visualization.
------------------------------------------------------------------------------------

------------------------------------------------------------------------------------
A largue number of the power spectra analysis functions were taken from James
Beattie's tools (https://github.com/AstroJames/PLASMAtools/tree/main).
------------------------------------------------------------------------------------


--------------
Usage Example (2D Spectra):
--------------
filename = '/Users/michaelgrehan/Downloads/data0500.vtu'
data, names = fast_vtu_reader(filename) 
t = data['time']
power_2D_b1 = compute_power_spectrum_2D(data, data['b1'], Ngrid_x= Ngrid,Ngrid_y=Ngrid)
power_2D_b2 = compute_power_spectrum_2D(data, data['b2'], Ngrid_x= Ngrid,Ngrid_y=Ngrid)
power_2D = power_2D_b1 + power_2D_b2
k_modes, power_1D = spherical_integrate_2D(power_2D)

plt.loglog(k_modes, power_1D, '-')
plt.title(f'$t/t_c = {t:.1f}$')
plt.xlim(left=1)
plt.xlabel("$k_\\perp L/2 \\pi$")
plt.ylabel("$B_\\perp^2(k_\\perp)$")
plt.show()

--------------
Usage Example (3D Spectra):
--------------
filename = '/Users/michaelgrehan/Downloads/data0450.vtu'
data, names = fast_vtu_reader(filename, blocks=False) 
t = data['time']
Ngrid = 128
power_3D_b1 = compute_power_spectrum_3D(data, data['b1'], Ngrid_x= Ngrid,Ngrid_y=Ngrid,Ngrid_z=Ngrid)
power_3D_b2 = compute_power_spectrum_3D(data, data['b2'], Ngrid_x= Ngrid,Ngrid_y=Ngrid,Ngrid_z=Ngrid)
power_3D = power_3D_b1 + power_3D_b2
k_modes, power = spherical_integrate_3D(power_3D)
kperp, kpar, power_cyl = cylindrical_integrate(power_3D)

plt.loglog(k_modes, power, '-')
plt.title(f'$t/t_c = {t:.1f}$')
plt.xlim(left=1)
plt.xlabel("$k L/2 \\pi$")
plt.ylabel("$B_\\perp^2(k)$")
plt.show()

plt.loglog(kperp, power_cyl[:,0], '-')
plt.title(f'$t/t_c = {t:.1f}$')
plt.xlim(left=1)
plt.xlabel("$k_\\perp L/2 \\pi$")
plt.ylabel("$B_\\perp^2(k_\\perp)$")
plt.show()

plt.loglog(kpar, power_cyl[0,:], '-')
plt.title(f'$t/t_c = {t:.1f}$')
plt.xlim(left=1)
plt.xlabel("$k_\\parallel L/2 \\pi$")
plt.ylabel("$B_\\perp^2(k_\\parallel)$")
plt.show()

"""




import numpy as np
import finufft

def compute_power_spectrum_2D(data: dict, field: np.ndarray, 
                              Ngrid_x: int = 256, Ngrid_y: int = 256, norm: str = "forward",
                              x_range: tuple = None, y_range: tuple = None,
                              cores:int=1) -> np.ndarray:
    """
    Computes power spectrum of a 2D scalar field sampled at non-uniform (x, y) points using NUFFT.
    
    Args:
        data (dict): bhac vtu reader output
        field (np.ndarray): 1D array of field values at (x, y) points.
        Ngrid_x (int): The number of grid points along the x-axis (default is 256).
        Ngrid_y (int): The number of grid points along the y-axis (default is 256).        
        norm (str): Normalization mode ("forward" or "backward").
        x_range (tuple, optional): A tuple (xmin, xmax) to limit the interpolation to the specified x bounds. If None, no limits are applied.
        y_range (tuple, optional): A tuple (ymin, ymax) to limit the interpolation to the specified y bounds. If None, no limits are applied.
        cores (int): Number of cores to use in parallel when interpolating. Use 0 or -1 to use all cores.

    Returns:
        power_spectrum_2D (np.ndarray): 2D power spectrum (shifted to center).
    """
    grid_size = (Ngrid_x, Ngrid_y)
    x = data['center_x']  # 1D array
    y = data['center_y']  # 1D array
    
    # Create initial mask for both x and y
    mask = np.ones(x.shape, dtype=bool)

    # Apply spatial filtering based on the provided x_range
    if x_range is not None:
        x_mask = (x >= x_range[0]) & (x <= x_range[1])
        mask &= x_mask  # Combine with the overall mask

    # Apply spatial filtering based on the provided y_range
    if y_range is not None:
        y_mask = (y >= y_range[0]) & (y <= y_range[1])
        mask &= y_mask  # Combine with the overall mask

    # Filter the center_x, center_y, and variable data based on the combined mask
    x = x[mask]
    y = y[mask]
    field = field[mask] 
    
    diffx = np.abs(np.diff(x))
    diffx = np.min(diffx[diffx != 0])
    diffy = np.abs(np.diff(y))
    diffy = np.min(diffy[diffy != 0])
    ncellsx = np.abs(x.max() - x.min() )/diffx
    ncellsy = np.abs(y.max() - y.min() )/diffy

    
    # Non-uniform case (using finufft)
    # Normalize coordinates to [-π, π]
    x_norm = 2 * np.pi * (x - x.min()) / (x.ptp()) - np.pi
    y_norm = 2 * np.pi * (y - y.min()) / (y.ptp()) - np.pi
        
    # Compute NUFFT
    fft_result = finufft.nufft2d1(x_norm, y_norm, field.astype(np.complex64), grid_size, modeord=1, isign=-1,nthreads=cores)
    
    # Apply normalization
    if norm == "forward":
        fft_result /= ncellsx * ncellsy  # Match FFT's "forward" scaling
                                         # Perserves Parseval's theorem
    
    # Shift and compute power
    power_spectrum_shifted = np.fft.fftshift(np.abs(fft_result)**2)
    return power_spectrum_shifted


def compute_power_spectrum_3D(data: dict, field: np.ndarray, 
                              Ngrid_x: int = 256, Ngrid_y: int = 256, Ngrid_z: int = 256, norm: str = "forward",
                              x_range: tuple = None, y_range: tuple = None, z_range: tuple = None,
                              cores:int=1) -> np.ndarray:
    """
    Computes power spectrum of a 2D scalar field sampled at non-uniform (x, y) points using NUFFT.
    
    Args:
        data (dict): bhac vtu reader output
        field (np.ndarray): 1D array of field values at (x, y) points.
        Ngrid_x (int): The number of grid points along the x-axis (default is 256).
        Ngrid_y (int): The number of grid points along the y-axis (default is 256).
        Ngrid_z (int): The number of grid points along the z-axis (default is 256).                
        norm (str): Normalization mode ("forward" or "backward").
        x_range (tuple, optional): A tuple (xmin, xmax) to limit the interpolation to the specified x bounds. If None, no limits are applied.
        y_range (tuple, optional): A tuple (ymin, ymax) to limit the interpolation to the specified y bounds. If None, no limits are applied.
        z_range (tuple, optional): A tuple (zmin, zmax) to limit the interpolation to the specified y bounds. If None, no limits are applied.
        cores (int): Number of cores to use in parallel when interpolating. Use 0 or -1 to use all cores.

    Returns:
        power_spectrum_3D (np.ndarray): 3D power spectrum (shifted to center).
    """
    grid_size = (Ngrid_x, Ngrid_y, Ngrid_z)
    x, y, z = calculate_cell_centers_3d(data)
    
    # Create initial mask for both x and y
    mask = np.ones(x.shape, dtype=bool)

    # Apply spatial filtering based on the provided x_range
    if x_range is not None:
        x_mask = (x >= x_range[0]) & (x <= x_range[1])
        mask &= x_mask  # Combine with the overall mask

    # Apply spatial filtering based on the provided y_range
    if y_range is not None:
        y_mask = (y >= y_range[0]) & (y <= y_range[1])
        mask &= y_mask  # Combine with the overall mask
        
    # Apply spatial filtering based on the provided z_range
    if z_range is not None:
        z_mask = (z >= z_range[0]) & (z <= z_range[1])
        mask &= z_mask  # Combine with the overall mask

    # Filter the center_x, center_y, and variable data based on the combined mask
    x = x[mask]
    y = y[mask]
    z = z[mask]
    field = field[mask] 
    
    diffx = np.abs(np.diff(x))
    diffx = np.min(diffx[diffx != 0])
    diffy = np.abs(np.diff(y))
    diffy = np.min(diffy[diffy != 0])
    diffz = np.abs(np.diff(z))
    diffz = np.min(diffz[diffz != 0])
    ncellsx = np.abs(x.max() - x.min() )/diffx
    ncellsy = np.abs(y.max() - y.min() )/diffy
    ncellsz = np.abs(z.max() - z.min() )/diffz

    
    # Non-uniform case (using finufft)
    # Normalize coordinates to [-π, π]
    x_norm = 2 * np.pi * (x - x.min()) / (x.ptp()) - np.pi
    y_norm = 2 * np.pi * (y - y.min()) / (y.ptp()) - np.pi
    z_norm = 2 * np.pi * (z - z.min()) / (z.ptp()) - np.pi
        
    # Compute NUFFT
    fft_result = finufft.nufft3d1(x_norm, y_norm, z_norm, field.astype(np.complex64), grid_size, modeord=1, isign=-1, nthreads=cores)
    
    # Apply normalization
    if norm == "forward":
        fft_result /= ncellsx * ncellsy * ncellsz  # Match FFT's "forward" scaling
                                         # Perserves Parseval's theorem
    
    # Shift and compute power
    power_spectrum_shifted = np.fft.fftshift(np.abs(fft_result)**2)
    return power_spectrum_shifted


def spherical_integrate_2D(Data: np.ndarray, 
                           bins: int = None) -> tuple:
    """
    Integrates the 2D power spectrum over spherical shells of constant k.
    
    Args:
        data: The 2D power spectrum.
        bins: Number of bins for radial integration. Uses Nyquist limit if None.
    
    Returns:
        k_modes: Bin centers of the k modes.
        radial_sum: Integrated power in each bin.
        
    Author: James Beattie
    """
    y, x = np.indices(Data.shape)
    center = np.array([(i - 1) / 2.0 for i in Data.shape])
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)  

    N = Data.shape[0]
    if bins is None:
        bins = N // 2

    bin_edges = np.linspace(0.5, bins, bins + 1)
    bin_indices = np.digitize(r, bin_edges)

    radial_sum = np.zeros(bins)
    for i in range(1, bins + 1):
        mask = bin_indices == i
        radial_sum[i-1] = np.sum(Data[mask])

    k_modes = (bin_edges[:-1] + bin_edges[1:]) / 2

    return k_modes, radial_sum




def spherical_integrate_3D(Data: np.ndarray, 
                        bins: int = None) -> tuple:
    """
    The spherical integrate function takes the 3D power spectrum and integrates
    over spherical shells of constant k. The result is a 1D power spectrum.
    
    It has been tested to reproduce the 1D power spectrum of an input 3D Gaussian
    random field.
    
    It has been tested to maintain the correct normalisation of the power spectrum
    i.e., the integral over the spectrum. For small grids (number of bins) the normalisation
    will be off by a small amount (roughly factor 2 for 128^3 with postive power-law indexes). 
    This is because the frequencies beyond the Nyquist limit are not included in the radial 
    integration. This is not a problem for grids of 256^3 or larger, or for k^-a style spectra,
    which are far more commonly encountered, and the normalisation is closer to 1/10,000 numerical
    error
    
    Args:
        data: The 3D power spectrum
        bins: The number of bins to use for the radial integration. 
              If not specified, the Nyquist limit is used (as should always be the case, anyway).

    Returns:
        k_modes: The k modes corresponding to the radial integration
        radial_sum: The radial integration of the 3D power spectrum (including k^2 correction)
        
    Author: James Beattie
    """
    z, y, x = np.indices(Data.shape)
    center = np.array([(i - 1) / 2.0 for i in Data.shape])
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2 + (z - center[2])**2)

    N = Data.shape[0]
    if not bins:
        bins = N // 2

    bin_edges = np.linspace(0.5, bins, bins+1)

    # Use np.digitize to assign each element to a bin
    bin_indices = np.digitize(r, bin_edges)

    # Compute the radial profile
    radial_sum = np.zeros(bins)
    for i in range(1, bins+1):
        mask = bin_indices == i
        radial_sum[i-1] = np.sum(Data[mask])

    # Generate the spatial frequencies with dk=1
    # Now k_modes represent the bin centers
    k_modes = np.ceil((bin_edges[:-1] + bin_edges[1:])/2)

    return k_modes, radial_sum



def cylindrical_integrate_3D(data: np.ndarray, 
                          bins_perp: int = 0,
                          bins_para:  int = 0) -> tuple:
    """
    The cylindrical integrate function takes the 3D power spectrum and integrates
    over cylindrical shells of constant k in a plane, and then a 1D spectrum along
    the the remaining dimension. The result is a 2D power spectrum of k_perp and k_par
    modes.
    
    
   Args:
        data     : The 3D power spectrum
        bins_perp: Number of bins for k_perpendicular
        bins_para: Number of bins for k_parallel

    Returns:
        k_perp_modes    : k_perpendicular modes
        k_para_modes    : k_parallel modes
        cylindrical_sum : cylindrical integration of the 3D power spectrum
                          note that k_perp is on axis 0, k_para is on axis 1
                          
    Author: James Beattie
    """
    
    z, y, x = np.indices(data.shape)
    center  = np.array([(i - 1) / 2.0 for i in data.shape])
    k_perp  = np.sqrt((x - center[0])**2 + (y - center[1])**2)  # Cylindrical radius
    k_para  = np.abs(z - center[2])                             # Distance from the plane

    N = data.shape[0]
    if bins_perp == 0:
        bins_perp = N // 2
    if bins_para == 0:
        bins_para = N // 2

    # initailize cylindrical sum
    cylindrical_sum = np.zeros((bins_perp, bins_para))

    # define bin edges (note starting at 0.5 to avoid binning the zero mode)
    bin_edges_perp = np.linspace(0, bins_perp, bins_perp+1)
    bin_edges_para = np.linspace(0, bins_para, bins_para+1)

    # Vectorized bin assignment
    bin_indices_perp = np.digitize(k_perp, bin_edges_perp) - 1
    bin_indices_para = np.digitize(k_para, bin_edges_para) - 1

    # Create 2D linear indices
    linear_indices = bin_indices_perp + bin_indices_para * bins_perp

    # Use np.bincount for efficient summation
    cylindrical_sum = np.bincount(linear_indices.ravel(), 
                                  weights=data.ravel(), 
                                  minlength=bins_perp * bins_para)
    
    # Ensure that the length matches the expected size
    cylindrical_sum = cylindrical_sum[:bins_perp * bins_para]
    cylindrical_sum = cylindrical_sum.reshape((bins_perp, bins_para))
    k_perp_modes    = (bin_edges_perp[:-1] + bin_edges_perp[1:]) / 2
    k_para_modes    = (bin_edges_para[:-1] + bin_edges_para[1:]) / 2

    return k_perp_modes, k_para_modes, cylindrical_sum




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
    

    return cell_centers_x, cell_centers_y, cell_centers_z
