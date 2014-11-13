#!/usr/bin/env python
"""calculate_figure_data.py

Calculates the data required for generating the plots for the EMD model.
Specifically: Fig 3B, Fig S3B

All data files are saved in the hdf5 file format. To generate them just
run:

    python calculate_figure_data.py  # generates all data files

or, dependent on which dataset you want to generate:

    python calculate_figure_data.py --fig3Btop
                                    --fig3Bbottom
                                    --figS3B
"""
import itertools
import multiprocessing
import os
import signal

import h5py
import numpy as np
import progressbar

from emd.stimulus import CheckerboardFGCorrect, SinusoidalFigureGroundStimulus
from emd.model import VisualSystemModel

# Global variables
# ================
N_CPUS = max(1, multiprocessing.cpu_count() - 1)


# Helper functions
# ================
def euler_maruyama_integrate(model, ic_fg, ic_bg, T_max, noise_sigma,
                             coupling, sample_rate_hz, omega_fg_lim, cells):
    """euler_maruyama_integrate(...)

    Euler-Maruyama integration of the model output torque.

    Parameters
    ----------
    model : emd.model.VisualSystemModel
      an instance of the model
    ic_fg : (float, float, float)
      the initial conditions for the figure (x, v, a)
    ic_bg : (float, float, float)
      the initial conditions for the background (x, v, a)
    T_max : float
      the simulation time limit (seconds)
    noise_sigma : float
      the standard deviation of the (gaussian) noise
    coupling : float
      the coupling constant used to calculate angular velocity from torque
      v[i+1] = v[i] + coupling * a[i]
    sample_rate : float
      the sampling rate of the model (Hz)
    omega_fg_lim : float
      a cutoff limit for the figure velocity (use np.inf for no limit)
    cells : (bool, bool)
      specifies which cells to use for calculating the torque

    Returns
    -------
    (t, x, v, a) :
      time, position, velocity and torque

    """
    dt = 1. / sample_rate_hz
    t = np.arange(0., T_max, dt)

    # generate the noise array
    noise = np.random.normal(0, scale=noise_sigma, size=t.shape)

    # preallocate empty arrays for position, velocity and acceleration
    # for figure and background
    phi_fg = np.zeros_like(t)
    omega_fg = np.zeros_like(t)
    alpha_fg = np.zeros_like(t)

    # set starting position for figure
    phi_fg[0] = ic_fg[0]
    omega_fg[0] = ic_fg[1]
    alpha_fg[0] = ic_fg[2]
    phi_bg = ic_bg[0]
    omega_bg = ic_bg[1]
    alpha_bg = ic_bg[2]

    # integrate the model torque
    for i, (current_time, dt) in enumerate(zip(t[1:], np.diff(t))):
        hsl, hsr = model.step(phi_fg[i], phi_bg, current_time)
        alpha_fg[i+1] = cells[0] * hsl - cells[1] * hsr
        omega_fg[i+1] = omega_fg[i] + coupling * (alpha_fg[i] * dt + noise[i] * np.sqrt(dt))
        # limit velocity to omega_fg_lim
        if not (-omega_fg_lim < omega_fg[i+1] < omega_fg_lim):
            omega_fg[i+1] = max(min(omega_fg[i+1], omega_fg_lim), -omega_fg_lim)
        phi_fg[i+1] = phi_fg[i] + omega_fg[i] * dt
        omega_bg += alpha_bg * dt
        phi_bg += omega_bg * dt

    return t, phi_fg, omega_fg, alpha_fg


def _parallel_wrapper(args):
    """parallel_wrapper(args)

    Helper function for multiprocessing.Pool.imap
    Needed because only a single argument can be passed via imap.

    """
    # unpack argument.
    stimulus_params, model_params, ic_fg, ic_bg, simulation_params = args

    # create model instance
    StimulusClass = stimulus_params.pop('type')
    stimulus = StimulusClass(**stimulus_params)
    model = VisualSystemModel(stimulus, simulation_params['sample_rate_hz'], **model_params)

    # integrate
    t, xf, vf, af = euler_maruyama_integrate(model, ic_fg, ic_bg, **simulation_params)
    # complete arrays are returned here.
    # Note: postprocessing should be done already here, to reduce overhead.
    return (t, xf, vf, af)


def calculate_parallel(args, postprocess=lambda x: x, label=''):
    """calculate_parallel(args, postprocessing)

    Calculates the simulation in parallel processes

    Parameters
    ----------
    args : [tuple, ...]
      List of arguments used in _parallel_wrapper
    postprocess : function
      Postprocessing function for time, position, velocity, torque tuple

    Returns
    -------
    data : list
      List of accumulated results from _parallel_wrapper (Same order as args)

    """
    # pretty display: create a progressbar
    widgets = ['Calculating ', label, ': ', progressbar.Percentage(),
               ' ', progressbar.Bar(), ' ', progressbar.ETA()]
    pbar = progressbar.ProgressBar(maxval=len(args), widgets=widgets).start()

    # required init function for child processes
    def init_worker():
        """make childprocess ignore SIGINT"""
        signal.signal(signal.SIGINT, signal.SIG_IGN)

    data = []
    pool = multiprocessing.Pool(processes=N_CPUS, initializer=init_worker)
    # start the processing and stop on CTRL-C
    # Note: catching CTRL-C by the parent might take a while if many processes
    #       are running
    try:
        for trajectory in pbar(pool.imap(_parallel_wrapper, args)):
            data.append(postprocess(*trajectory))
    except KeyboardInterrupt:
        print "\nCaught KeyboardInterrupt (CTRL-C). Terminating worker processes."
        pool.terminate()
        exit(1)

    return data


# Data functions
# ==============
def generate_fig3Btop_dataframe(simulation_params, model_params, stimulus_params):
    # simulation parameters for 3Btop
    x_fg = np.linspace(-np.pi, np.pi, 50)
    v_bg = np.linspace(0., -0.1, 25)

    # generate the list of argument tuples for parallel processing
    pargs = []
    for phi_fg, omega_bg in itertools.product(x_fg, v_bg):
        pargs.append((stimulus_params, model_params, (phi_fg, 0, 0), (0, omega_bg, 0), simulation_params))

    print "Fig3B top: This will take approximatley %0.1f minutes" % (1150. / N_CPUS)
    data = calculate_parallel(pargs, label="Fig3B top", postprocess=lambda t, x, v, a: np.mean(v))

    if os.path.exists('data_figure3Btop.hdf5'):
        os.unlink('data_figure3Btop.hdf5')
    hdf = h5py.File('data_figure3Btop.hdf5', 'w-')
    hdf.create_dataset('fg_positions', data=x_fg).attrs['unit'] = 'rad'
    hdf.create_dataset('bg_velocities', data=v_bg).attrs['unit'] = 'rad/s'
    hdf.create_dataset('average_fg_velocity', data=np.array(data), shape=(len(x_fg), len(v_bg))).attrs['unit'] = 'rad/s'
    hdf.close()
    print 'data_figure3Btop.hdf5 saved.'


def generate_fig3Bbottom_dataframe(simulation_params, model_params, stimulus_params):
    # simulation parameters for 3Bbottom
    x_fg = 2*np.pi * np.random.random(20) - np.pi
    a_bg = [-0.0004]

    # generate the list of argument tuples for parallel processing
    pargs = []
    for phi_fg, alpha_bg in itertools.product(x_fg, a_bg):
        pargs.append((stimulus_params, model_params, (phi_fg, 0, 0), (0, 0, alpha_bg), simulation_params))

    print "Fig3B bottom: This will take approximatley %0.1f minutes" % (2500. / N_CPUS)
    data = calculate_parallel(pargs, label="Fig3B bottom", postprocess=lambda t, x, v, a: (t, x))

    if os.path.exists('data_figure3Bbottom.hdf5'):
        os.unlink('data_figure3Bbottom.hdf5')
    hdf = h5py.File('data_figure3Bbottom.hdf5', 'w-')
    hdf.create_dataset('time', data=data[0][0]).attrs['unit'] = 's'
    hdf.create_dataset('bg_velocity', data=data[0][0]*a_bg[0]).attrs['unit'] = 'rad/s'
    hdf.create_dataset('fg_position', data=np.array([x for t, x in data])).attrs['unit'] = 'rad'
    hdf.close()
    print 'data_figure3Bbottom.hdf5 saved.'


def generate_figS3B_dataframe(simulation_params, model_params, stimulus_params):
    # simulation parameters for S3B
    tf_bg = np.array([-50, -40, -30, -25, -20, -15, -10, -5, -2.5, -1, -0.5,
                      -0.1, 0.1, 0.5, 1, 2.5, 5, 10, 15, 20, 25, 30, 40, 50])
    v_bg = tf_bg * stimulus_params['ground_spatial_wavelength']

    # generate the list of argument tuples for parallel processing
    pargs = []
    for tf, omega_bg in zip(tf_bg, v_bg):
        simp = simulation_params.copy()
        simp['T_max'] = 10. / abs(tf)
        pargs.append((stimulus_params, model_params, (0, 0, 0), (0, omega_bg, 0), simp))

    print "FigS3B: This will take approximatley %0.1f minutes" % (500. / N_CPUS)
    data = calculate_parallel(pargs, label="FigS3B", postprocess=lambda t, x, v, a: np.mean(a[len(a) / 2:]))

    if os.path.exists('data_figureS3B.hdf5'):
        os.unlink('data_figureS3B.hdf5')
    hdf = h5py.File('data_figureS3B.hdf5', 'w-')
    hdf.create_dataset('temporal_frequency', data=tf_bg).attrs['unit'] = 'Hz'
    hdf.create_dataset('average_hs_cell_response', data=data).attrs['unit'] = 'a.u.'
    hdf.close()
    print 'data_figureS3B.hdf5 saved.'


if __name__ == '__main__':

    import sys

    # Default parameters
    # ==================
    MODEL_PARAMS = {
            'variant': 7,
            'Ei': -0.5,
            }
    STIMULUS_PARAMS = {
            'type': CheckerboardFGCorrect,
            'stripe_width_radians': np.radians(20),
            'stripe_spatial_wavelength': np.radians(20),
            'ground_spatial_wavelength': np.radians(20),
            }
    SIMULATION_PARAMS = {
            'sample_rate_hz': 1000.,
            'T_max': 20.,
            'noise_sigma': 0.001,
            'coupling': -3000,
            'omega_fg_lim': np.inf,
            'cells': (True, True)
            }

    # No cmdline parameters means generate all data
    args = sys.argv[1:]
    RUNALL = not args

    # Ensures that blurmatrix is created in this process
    _ = VisualSystemModel(SinusoidalFigureGroundStimulus(**{k:v for k,v in STIMULUS_PARAMS.items() if k != 'type'}),
                          SIMULATION_PARAMS['sample_rate_hz'], **MODEL_PARAMS)

    if '--fig3Btop' in args or RUNALL:
        np.random.seed(0)
        SIMULATION_PARAMS['T_max'] = 2.
        generate_fig3Btop_dataframe(SIMULATION_PARAMS, MODEL_PARAMS, STIMULUS_PARAMS)

    if '--fig3Bbottom' in args or RUNALL:
        np.random.seed(0)
        SIMULATION_PARAMS['T_max'] = 250.
        generate_fig3Bbottom_dataframe(SIMULATION_PARAMS, MODEL_PARAMS, STIMULUS_PARAMS)

    if '--figS3B' in args or RUNALL:
        np.random.seed(0)
        STIMULUS_PARAMS['type'] = SinusoidalFigureGroundStimulus
        STIMULUS_PARAMS['ground_spatial_wavelength'] = np.radians(10.)
        STIMULUS_PARAMS['stripe_width_radians'] = 0.
        SIMULATION_PARAMS['cells'] = (True, False)
        generate_figS3B_dataframe(SIMULATION_PARAMS, MODEL_PARAMS, STIMULUS_PARAMS)
