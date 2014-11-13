import numpy as np
import h5py


def plot_figure3B_top(ax, ylim, cmap, vlim):
    """plots the emd bifurcation.

    returns image data from imshow
    """
    try:
        hdf = h5py.File('data_figure3Btop.hdf5', 'r')
    except IOError:
        print "requires 'data_figure3Btop.hdf5'. please run calculate_figure_data.py first."
        exit(1)

    data = hdf['average_fg_velocity'].value
    fg_pos = hdf['fg_positions'].value
    bg_vel = hdf['bg_velocities'].value

    fg_pos = np.degrees(fg_pos)
    bg_vel = np.degrees(bg_vel)

    extent = [bg_vel[0], bg_vel[-1], fg_pos[0], fg_pos[-1]]

    im = ax.imshow(data, vmin=-vlim, vmax=vlim, aspect='auto', origin='lower',
                   extent=extent, interpolation='nearest', cmap=cmap)

    ax.set_xlim(0, -5)
    ax.set_ylim(*ylim)

    return im


if __name__ == '__main__':

    import matplotlib.pyplot as plt

    ax = plt.subplot2grid((1, 1), (0, 0))
    im = plot_figure3B_top(ax, (0, 180), plt.get_cmap('PiYG'), vlim=0.4)
    plt.colorbar(im)

    plt.show()
