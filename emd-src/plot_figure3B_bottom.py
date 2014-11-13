#!/usr/bin/env python
import numpy as np
import h5py
import emd.wrapping


def plot_figure3B_bottom(ax, ylim, cmap, vlim=None):
    """plots figure 3B bottom
    """
    try:
        hdf = h5py.File('data_figure3Bbottom.hdf5', 'r')
    except IOError:
        print "requires 'data_figure3Bbottom.hdf5'. please run calculate_figure_data.py first."
        exit(1)
    # get data
    bg_vel = hdf['bg_velocity'].value
    fg_pos_flies = hdf['fg_position'].value

    bg_vel = np.degrees(bg_vel)
    fg_pos_flies = np.degrees(emd.wrapping.wrap_to(fg_pos_flies, -np.pi, np.pi))

    BINY = 180
    BINX = 35

    H_sum = np.zeros((BINY, BINX))
    for fg_pos in fg_pos_flies:
        # histogram2d needs a monotonic increase in range
        H, _, _ = np.histogram2d(fg_pos, bg_vel,
                                 range=[[-180, 180], [-5, 0]],
                                 bins=(BINY, BINX), normed=True)
        H_sum += H
    H_sum = H_sum / fg_pos_flies.shape[0]
    H_sum = H_sum[::-1, ::-1]  # flip back (histogram2d requirements)

    extent = [0, -5, 180, -180]

    im = ax.imshow(H_sum, extent=extent, interpolation='nearest',
                   vmin=0, vmax=vlim, origin='lower', aspect='auto', cmap=cmap)
    ax.set_ylim(*ylim)
    ax.set_xlim(0, -5)

    return im


if __name__ == '__main__':

    import matplotlib.pyplot as plt

    ax = plt.subplot2grid((1, 1), (0, 0))
    plot_figure3B_bottom(ax, (-20, 60), plt.get_cmap('Greys'), vlim=0.015)

    plt.show()
