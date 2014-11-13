
import numpy as np
import scipy.signal

from filters import FilterMaker
from blurmatrix import FastBlurMatrix as BlurMatrix
from receptivefield import ReceptiveField

import functools


class VisualSystemModel(object):

    def __init__(self, stimulus, hz, **kwargs):
        """

        Parameters
        ----------
        stimulus :
        hz :


        Keyword Arguments
        -----------------

        """

        def setkwarg(x, v):
            if not x in kwargs: kwargs[x] = v

        VARIANT = kwargs.pop('variant', 7)
        if VARIANT == 1:
            setkwarg('filter_type_earlyvis', 'LP')
            setkwarg('tau_lp_earlyvis', 0.008)
            setkwarg('filter_type_arm1', 'LP')
            setkwarg('tau_lp_arm1', 0.035)
            setkwarg('filter_type_arm2', 'NONE')
            setkwarg('tau_hp_arm2', 0.035)
            setkwarg('output_type_spatialweighting', 'LINEAR')
            setkwarg('output_type_summation', 'LIN')
            setkwarg('output_type_lowpass', 'NONE')
            setkwarg('tau_mem', 0.008)
        elif VARIANT == 2:
            setkwarg('filter_type_earlyvis', 'LP')
            setkwarg('tau_lp_earlyvis', 0.008)
            setkwarg('filter_type_arm1', 'LP')
            setkwarg('tau_lp_arm1', 0.035)
            setkwarg('filter_type_arm2', 'HP')
            setkwarg('tau_hp_arm2', 0.035)
            setkwarg('output_type_spatialweighting', 'LINEAR')
            setkwarg('output_type_summation', 'LIN')
            setkwarg('output_type_lowpass', 'NONE')
            setkwarg('tau_mem', 0.008)
        elif VARIANT == 3:
            setkwarg('filter_type_earlyvis', 'LP')
            setkwarg('tau_lp_earlyvis', 0.008)
            setkwarg('filter_type_arm1', 'LP')
            setkwarg('tau_lp_arm1', 0.006)
            setkwarg('filter_type_arm2', 'NONE')
            setkwarg('tau_hp_arm2', 0.035)
            setkwarg('output_type_spatialweighting', 'LINEAR')
            setkwarg('output_type_summation', 'LIN')
            setkwarg('output_type_lowpass', 'NONE')
            setkwarg('tau_mem', 0.008)
        elif VARIANT == 4:
            setkwarg('filter_type_earlyvis', 'LMC')
            setkwarg('tau_lp_earlyvis', 0.008)
            setkwarg('filter_type_arm1', 'LP')
            setkwarg('tau_lp_arm1', 0.015)
            setkwarg('filter_type_arm2', 'NONE')
            setkwarg('tau_hp_arm2', 0.035)
            setkwarg('output_type_spatialweighting', 'LINEAR')
            setkwarg('output_type_summation', 'LIN')
            setkwarg('output_type_lowpass', 'NONE')
            setkwarg('tau_mem', 0.008)
        elif VARIANT == 5:
            setkwarg('filter_type_earlyvis', 'LMC')
            setkwarg('tau_lp_earlyvis', 0.008)
            setkwarg('filter_type_arm1', 'LP')
            setkwarg('tau_lp_arm1', 0.010)
            setkwarg('filter_type_arm2', 'HP')
            setkwarg('tau_hp_arm2', 0.060)
            setkwarg('output_type_spatialweighting', 'LINEAR')
            setkwarg('output_type_summation', 'LIN')
            setkwarg('output_type_lowpass', 'NONE')
            setkwarg('tau_mem', 0.008)
        elif VARIANT == 6:
            setkwarg('filter_type_earlyvis', 'LMC')
            setkwarg('tau_lp_earlyvis', 0.008)
            setkwarg('filter_type_arm1', 'LP')
            setkwarg('tau_lp_arm1', 0.010)
            setkwarg('filter_type_arm2', 'HP')
            setkwarg('tau_hp_arm2', 0.060)
            setkwarg('output_type_spatialweighting', 'NONLINEAR')
            setkwarg('output_type_summation', 'MEM')
            setkwarg('output_type_lowpass', 'NONE')
            setkwarg('tau_mem', 0.008)
        elif VARIANT == 7:
            setkwarg('filter_type_earlyvis', 'LMC')
            setkwarg('tau_lp_earlyvis', 0.008)
            setkwarg('filter_type_arm1', 'LP')
            setkwarg('tau_lp_arm1', 0.010)
            setkwarg('filter_type_arm2', 'HP')
            setkwarg('tau_hp_arm2', 0.060)
            setkwarg('output_type_spatialweighting', 'NONLINEAR')
            setkwarg('output_type_summation', 'MEM')
            setkwarg('output_type_lowpass', 'C')
            setkwarg('tau_mem', 0.008)

        FILTER_TYPE_EARLYVIS    = kwargs.pop('filter_type_earlyvis')
        TAU_LP_EARLYVIS         = kwargs.pop('tau_lp_earlyvis')
        FILTER_TYPE_ARM1        = kwargs.pop('filter_type_arm1')
        TAU_LP_ARM1             = kwargs.pop('tau_lp_arm1')
        FILTER_TYPE_ARM2        = kwargs.pop('filter_type_arm2')
        TAU_HP_ARM2             = kwargs.pop('tau_hp_arm2')
        OUTPUT_TYPE_SPATIALWEIGHTING = kwargs.pop('output_type_spatialweighting')
        OUTPUT_TYPE_SUMMATION   = kwargs.pop('output_type_summation')
        OUTPUT_TYPE_LOWPASS     = kwargs.pop('output_type_lowpass')
        TAU_MEM                 = kwargs.pop('tau_mem')

        g0 = kwargs.pop('g0', 1295.0)
        E0 = kwargs.pop('E0', 0.0)
        Ee = kwargs.pop('Ee', 1.0)
        Ei = kwargs.pop('Ei', -0.95)

        FILTERS_NORMALIZE = kwargs.pop('filters_normalize', True)

        # stimulus -> .calculate_luminance( fg_angle, bg_angle )
        self.stimulus = stimulus

        # photoreceptor positions
        dphi = np.radians(2)
        self.phi = np.arange(-np.pi, np.pi, dphi)
        self.theta = np.linspace(-25*dphi, 25*dphi, 51)
        self.PHI, self.THETA = np.meshgrid( self.phi, self.theta )
        self.receptor_grid_shape = (len(self.theta), len(self.phi))

        # photoreceptor stimulus downsampling
        self.bm = BlurMatrix(self.stimulus.phi, self.stimulus.theta, self.phi, self.theta, dphi,
                             cachefile='.emd_blurmatrix_cached.npz')

        # receptive fields
        RFL = ReceptiveField(phi_c=np.radians( 15), sigma_theta=np.radians(33),
                             sigma_phi_a=np.radians(102), sigma_phi_b=np.radians( 45))
        RFR = ReceptiveField(phi_c=np.radians(-15), sigma_theta=np.radians(33),
                             sigma_phi_a=np.radians( 45), sigma_phi_b=np.radians(102))
        self.receptive_field_left = RFL(self.PHI, self.THETA)
        self.receptive_field_right = RFR(self.PHI, self.THETA)

        # filter functions
        fm = FilterMaker(hz, normalize=FILTERS_NORMALIZE, test_filter=True,
                         incremental_increase_till_valid=True)
        n_receptors = len(self.theta) * len(self.phi)

        """Setup the Earlyvision Filter"""

        if FILTER_TYPE_EARLYVIS == "LMC":
            b, a = fm.james_lmc()
        elif FILTER_TYPE_EARLYVIS == "LP":
            b, a = fm.iir_lowpass1(TAU_LP_EARLYVIS)
            b *= -1
        elif FILTER_TYPE_EARLYVIS == "NONE":
            b, a = (1,), (1,) # -> id()
        else:
            raise ValueError("filter_type_earlyvis not in ['LMC', 'LP', 'NONE']")

        self.filter_zi_earlyvis = np.vstack((scipy.signal.lfiltic(b, a, [0.]).astype(np.float64),)*n_receptors)
        self.filter_func_earlyvis = functools.partial(scipy.signal.lfilter, b, a, axis=1)

        """Setup both filter arms"""

        if FILTER_TYPE_ARM1 == "LP":
            b, a = fm.iir_lowpass1(TAU_LP_ARM1)
        elif FILTER_TYPE_ARM1 == "NONE":
            b, a = (1,), (1,) # -> id()
        else:
            raise ValueError("filter_type_arm1 not in ['LP', 'NONE']")

        self.filter_zi_filterA1 = np.vstack((scipy.signal.lfiltic(b, a, [0.]).astype(np.float64),)*n_receptors)
        self.filter_func_filterA1 = functools.partial(scipy.signal.lfilter, b, a, axis=1)

        if FILTER_TYPE_ARM2 == "HP":
            b, a = fm.iir_highpass1(TAU_HP_ARM2)
        elif FILTER_TYPE_ARM2 == "NONE":
            b, a = (1,), (1,) # -> id()
        else:
            raise ValueError("filter_type_arm2 not in ['HP', 'NONE']")

        self.filter_zi_filterA2 = np.vstack((scipy.signal.lfiltic(b, a, [0.]).astype(np.float64),)*n_receptors)
        self.filter_func_filterA2 = functools.partial(scipy.signal.lfilter, b, a, axis=1)

        """Setup output stage"""

        if OUTPUT_TYPE_SPATIALWEIGHTING == "NONLINEAR":
            self.spatial_weighting = lambda weights, hcout: np.sum( weights * np.clip(hcout, 0, np.inf) )
        elif OUTPUT_TYPE_SPATIALWEIGHTING == "LINEAR":
            self.spatial_weighting = lambda weights, hcout: np.sum( weights * hcout )
        else:
            raise ValueError("output_type_spatialweighting not in ['NONLINEAR', 'LINEAR']")

        if OUTPUT_TYPE_SUMMATION == "LIN":
            self.output_summation = lambda sge, sgi: sge - sgi
        elif OUTPUT_TYPE_SUMMATION == "MEM":
            self.output_summation = lambda sge, sgi: (E0 * g0 + Ee * sge + Ei * sgi) / (g0 + sge + sgi)
        else:
            raise ValueError("output_type_summation not in ['LIN', 'MEM']")

        if OUTPUT_TYPE_LOWPASS == "C":
            b, a = fm.iir_lowpass1(TAU_MEM)
        elif OUTPUT_TYPE_LOWPASS == "NONE":
            b, a = (1,), (1,) # -> id()
        else:
            raise ValueError("output_type_lowpass not in ['C', 'NONE']")

        self.filter_zi_output_left = scipy.signal.lfiltic(b, a, [0.]).astype(np.float64)
        self.filter_zi_output_right = scipy.signal.lfiltic(b, a, [0.]).astype(np.float64)
        self.filter_func_output = functools.partial(scipy.signal.lfilter, b, a)


    def step(self, stripe_angle, bg_angle, t, return_all=False):

        # Get Stimulus luminance
        luminance_input_2d = self.stimulus.calculate_luminance( stripe_angle, bg_angle, t=t)
        retinal_image = self.bm.blur(luminance_input_2d)
        retinal_image_flat = retinal_image.reshape((-1, 1))

        # earlyvision filter
        earlyvisA, self.filter_zi_earlyvis = self.filter_func_earlyvis(retinal_image_flat, zi=self.filter_zi_earlyvis)
        # filter arm 1
        filteredA1, self.filter_zi_filterA1 = self.filter_func_filterA1(earlyvisA, zi=self.filter_zi_filterA1)
        # filter arm 2
        filteredA2, self.filter_zi_filterA2 = self.filter_func_filterA2(earlyvisA, zi=self.filter_zi_filterA2)

        # generate halfcorrelators
        filteredA1_grid = filteredA1.reshape(self.receptor_grid_shape) # delayed
        filteredA2_grid = filteredA2.reshape(self.receptor_grid_shape) # undelayed
        half_correlators_rightward = filteredA1_grid * np.roll(filteredA2_grid, 1, axis=1)
        half_correlators_leftward  = filteredA2_grid * np.roll(filteredA1_grid, 1, axis=1)

        # LEFT HS CELL
        sge_left = self.spatial_weighting( self.receptive_field_left, half_correlators_leftward )
        sgi_left = self.spatial_weighting( self.receptive_field_left, half_correlators_rightward )

        # RIGHT HS CELL
        sge_right = self.spatial_weighting( self.receptive_field_right, half_correlators_rightward )
        sgi_right = self.spatial_weighting( self.receptive_field_right, half_correlators_leftward )

        # sum excitatory and inhibitory outputs        
        left_hs = np.atleast_1d(self.output_summation(sge_left, sgi_left))
        right_hs = np.atleast_1d(self.output_summation(sge_right, sgi_right))

        # filter output
        left_output, self.filter_zi_output_left = self.filter_func_output(left_hs, zi=self.filter_zi_output_left)
        right_output, self.filter_zi_output_right = self.filter_func_output(right_hs, zi=self.filter_zi_output_right)

        if not return_all:
            return left_output, right_output
        else:
            return {'stimulus_photoreceptor' : retinal_image,
                    'stimulus_earlyvis'      : earlyvisA.reshape(retinal_image.shape),
                    'stimulus_lp'            : filteredA1.reshape(retinal_image.shape),
                    'stimulus_hp'            : filteredA2.reshape(retinal_image.shape),
                    'half_correlators_rightward': half_correlators_rightward,
                    'half_correlators_leftward': half_correlators_leftward,
                    'right_halfcorrelator_e' : sge_right,
                    'right_halfcorrelator_i' : sgi_right,
                    'left_halfcorrelator_e'  : sge_left,
                    'left_halfcorrelator_i'  : sgi_left,
                    'right_hs'               : right_output,
                    'left_hs'                : left_output,
                    }



if __name__ == '__main__':

    from stimulus import SinusoidalFigureGroundStimulus
    import matplotlib.pyplot as plt
    import progressbar

    np.seterr(divide='ignore')

    # global variables for simulation
    SAMPLE_RATE = 3000 # Hz

    dt = 1.0/SAMPLE_RATE
    t = np.arange(0.0, 2.0, dt)

    stimulus = SinusoidalFigureGroundStimulus(t_start=0.3, t_stop=1.8, figure_width_radians=np.radians(22.5), bg_spatial_wavelength=np.radians(45))
    model = VisualSystemModel(stimulus, SAMPLE_RATE)

    fg_angles = t * 8 * np.pi / max(t)
    bg_angles = np.zeros_like(fg_angles)

    HSL, HSR = [], []

    pbar = progressbar.ProgressBar(maxval=t.size)
    pbar.start()

    for i, (fg, bg, current_time) in enumerate(zip(fg_angles, bg_angles, t)):
        pbar.update(i)
        hsl, hsr = model.step( fg, bg, current_time )
        HSL.append(hsl)
        HSR.append(hsr)

    pbar.finish()

    plt.plot(t, HSL)
    plt.plot(t, HSR)

    plt.show()

