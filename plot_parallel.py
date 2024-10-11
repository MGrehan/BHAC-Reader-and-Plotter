from read_MPG import *
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
import sys
import os

PATHTO = './'


def process_file(i):
    
    filenumber = '{}'.format(i).zfill(4)
    filename = 'data'+filenumber+'.vtu'
    data, names = fast_vtu_reader(filename, attr={'p', 'b1', 'b2', 'b3', 'rho'}, blocks=False)
    fig, ax = plt.subplots()
    
    xmin, xmax = -0.1, 0.1
    ymin, ymax = -0.05, 0.05

    # plot_raw_data_cells(data, 2*data['p']/(data['b1']**2 + data['b2']**2 + data['b3']**2), fig=fig, ax=ax, x_range=(xmin,xmax), 
    #                      y_range=(ymin,ymax), cmap='hot', label='$\\beta$', linewidths=0.2, edgecolors=None, orientation='horizontal',  
    #                      location='top', use_log_norm=True, pad=0.0)

    plot_raw_data_cells(data, data['rho'], fig=fig, ax=ax, x_range=(xmin,xmax), 
                        y_range=(ymin,ymax), cmap='hot', label='$\\rho$', linewidths=0.2, edgecolors=None, orientation='vertical',  
                        location='right', use_log_norm=True, pad=0.0, colorbar=True, vmin=1e-2, vmax=1e0)
    ax.set_xlabel('$x/L$')
    ax.set_ylabel('$y/L$')
    # ax.set_aspect('equal')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin ,ymax)
    dt = 0.1
    t = dt * i
    ax.set_title(f'$t/t_c = {t:.1f}$')    
    save_path = os.path.join(PATHTO+f"rho_"+filenumber+".png")
    while True:
        try:
            plt.savefig(save_path, format="png", bbox_inches="tight", dpi=300)
            plt.close()
            break
        except Exception as e:
            print_error(f"Error saving the plot for i={i}. Retrying...")
            time.sleep(1)



def plotxy_vtu():
    i0 = int(sys.argv[1])
    i1 = int(sys.argv[2])
    
    # Create a pool of worker processes
    num_workers = 16 #os.cpu_count()
    pool = Pool(processes=num_workers)
    
    for i in range(i0, i1 + 1):
        pool.apply_async(process_file, (i,))

    # Close the pool and wait for all processes to finish
    pool.close()
    pool.join()

if __name__ == "__main__":
    plotxy_vtu()
