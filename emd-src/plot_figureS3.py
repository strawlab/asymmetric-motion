import numpy as np
import h5py

from emd.receptivefield import ReceptiveField

# Receptive Field parameters
# ==========================
EMD_SIGMA_RF = {'sigma_theta': np.radians(33),
                'sigma_phi_a': np.radians(45),
                'sigma_phi_b': np.radians(102)}


def plot_figure_S3A(ax0, ax1, cmap):
    """plots the two recpetive fields from the emd model

    """
    receptive_field_left = ReceptiveField(phi_c=np.radians(15), **EMD_SIGMA_RF)
    receptive_field_right = ReceptiveField(phi_c=np.radians(-15), **EMD_SIGMA_RF)

    # generate grid data
    dphi = np.radians(2)
    phi = np.arange(-np.pi, np.pi + dphi, dphi)
    theta = np.linspace(-25 * dphi, 25 * dphi, 51)
    PHI, THETA = np.meshgrid(phi, theta)

    extent = map(np.degrees, [phi[0], phi[-1], theta[0], theta[-1]])

    ax0.imshow(receptive_field_left(PHI, THETA), extent=extent, interpolation='nearest', cmap=cmap)
    ax0.set_xlim([-180,180])
    ax0.set_xticklabels([])
    ax0.set_xticks([-180, -90, 0, 90, 180])
    ax0.set_ylim([-50,50])
    ax0.set_yticks([-50, 0, 50])
    ax0.set_ylabel('Elevation (deg)')

    im = ax1.imshow(receptive_field_right(PHI, THETA), extent=extent, interpolation='nearest', cmap=cmap)
    ax1.set_xlim([-180,180])
    ax1.set_xticks([-180, -90, 0, 90, 180])
    ax1.set_xlabel('Azimuth(deg)')
    ax1.set_ylim([-50,50])
    ax1.set_yticks([-50, 0, 50])
    ax1.set_ylabel('Elevation (deg)')

    return im


def plot_figure_S3B(ax):
    """plots the motion response vs temporal frequency

    """
    try:
        hdf = h5py.File('data_figureS3B.hdf5', 'r')
    except IOError:
        print "requires 'data_figureS3B.hdf5'. please run calculate_figure_data.py first."
        exit(1)
    # assgin data
    tf = hdf['temporal_frequency'].value
    hsl = hdf['average_hs_cell_response'].value
    hsl /= hsl.max()  # normalize
    hsr = hsl[::-1]  # the right cell is mirror symmetric to the left

    ax.plot(tf, hsl, lw=2, color='r', marker='o')
    ax.plot(tf, hsr, lw=2, color='b', marker='o')
    ax.set_xlim(min(tf), max(tf))
    ax.set_ylim(-0.5, 1.)
    ax.set_xticks([-50, -25, 0, 25, 50])
    ax.set_yticks([-0.5, 0, 0.5, 1])
    ax.set_ylabel('Normalized Response (a.u.)')
    ax.set_xlabel('Temporal Frequency (Hz)')

    return


if __name__ == '__main__':

    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax0 = plt.subplot2grid((2, 1), (0, 0))
    ax1 = plt.subplot2grid((2, 1), (1, 0))
    fig.subplots_adjust(bottom=0.4)
    cbar_ax = fig.add_axes([0.9, 0.4, 0.03, 0.5])

    im = plot_figure_S3A(ax0, ax1, plt.get_cmap('gray'))

    cbar = plt.colorbar(im, cax=cbar_ax, ticks=[0, 1])
    cbar.ax.set_yticklabels([0, 1])
    cbar.ax.set_ylabel('Response (a.u.)')

    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    plot_figure_S3B(ax)

    plt.show()
