from read_MPG import *
import numpy as np
import matplotlib.pyplot as plt


filename = 'data0075.vtu'
data, names = fast_vtu_reader(filename, attr={'p', 'b1', 'b2', 'b3', 'rho'}, blocks=False)


fig, ax = plt.subplots()
xmin, xmax = -0.1, 0.1
ymin, ymax = -0.05, 0.05

# plot_raw_data_cells(data, 2*data['p']/(data['b1']**2 + data['b2']**2 + data['b3']**2), fig=fig, ax=ax, x_range=(xmin,xmax), 
#                      y_range=(ymin,ymax), cmap='hot', label='$\\beta$', linewidths=0.2, edgecolors=None, orientation='horizontal',  
#                      location='top', use_log_norm=True, pad=0.0)

plot_raw_data_cells(data, data['rho'], fig=fig, ax=ax, x_range=(xmin,xmax), 
                     y_range=(ymin,ymax), cmap='hot', label='$\\rho$', linewidths=0.2, edgecolors=None, orientation='horizontal',  
                     location='top', use_log_norm=True, pad=0.0)
ax.set_xlabel('$x/L$')
ax.set_ylabel('$y/L$')
# ax.set_aspect('equal')
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin ,ymax)

plt.savefig("test.png")
