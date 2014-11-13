#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')

from plot_figure3B_bottom import plot_figure3B_bottom as p3Bb
from plot_figure3B_top import plot_figure3B_top as p3Bt
from plot_figureS3 import plot_figure_S3A as pS3A
from plot_figureS3 import plot_figure_S3B as pS3B


if __name__ == '__main__':

    import matplotlib.pyplot as plt

    plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    p3Bb(ax, (-20, 60), plt.get_cmap('Greys'), vlim=0.015)
    plt.savefig('figure_3B_bottom.pdf')

    plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    im = p3Bt(ax, (0, 180), plt.get_cmap('PiYG'), vlim=0.4)
    plt.colorbar(im)
    plt.savefig('figure_3B_top.pdf')

    fig = plt.figure()
    ax0 = plt.subplot2grid((2, 1), (0, 0))
    ax1 = plt.subplot2grid((2, 1), (1, 0))
    fig.subplots_adjust(bottom=0.4)
    cbar_ax = fig.add_axes([0.9, 0.4, 0.03, 0.5])
    im = pS3A(ax0, ax1, plt.get_cmap('gray'))
    cbar = plt.colorbar(im, cax=cbar_ax, ticks=[0, 1])
    cbar.ax.set_yticklabels([0, 1])
    cbar.ax.set_ylabel('Response (a.u.)')
    plt.savefig('figure_S3A.pdf')

    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    pS3B(ax)
    plt.savefig('figure_S3B.pdf')

