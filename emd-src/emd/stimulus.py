
import numpy as np
import scipy.signal

from .wrapping import wrap_to


class CheckerboardFGCorrect(object):

    def __init__(self, stripe_width_radians, 
                       ground_spatial_wavelength,
                       stripe_spatial_wavelength=np.radians(30),
                       t_start=-np.inf,
                       t_stop=np.inf,
                       elevation_extent=np.radians(60),
                       azimuth_samples=512):

        self.ground_spatial_wavelength = float(ground_spatial_wavelength)
        self.ground_circular_frequency = 2.*np.pi / self.ground_spatial_wavelength
        assert self.ground_circular_frequency % 1 == 0, "BG wavelength!"

        self.stripe_spatial_wavelength = float(stripe_spatial_wavelength)
        self.stripe_circular_frequency = 2.*np.pi / self.stripe_spatial_wavelength
        assert self.stripe_circular_frequency % 1 == 0, "FG wavelength!"

        self.stripe_width = float(stripe_width_radians)

        # stimulus start and end
        self.t_start = float(t_start)
        self.t_stop = float(t_stop)

        # phi and theta grid
        self.phi, self.dphi = np.linspace(-np.pi, np.pi, int(azimuth_samples), endpoint=False, retstep=True)
        _theta = np.arange(0., float(elevation_extent), self.dphi)
        self.theta = np.hstack((-1*_theta[:1:-1], _theta))


    def squarewave_edge_aware(self, x, dx, wavelength, zero_position):
        """returns a squarewave function starting at 'zero_position' with
        'wavelength' and amplitude between -1 and 1

        Array elements that contain an edge are correctly adjusted!
        """
        # assert that a checkerboard field is bigger than one array field
        array_field_width = dx
        ckbrd_field_width = wavelength / 2.
        assert ckbrd_field_width > array_field_width

        # calculate the checkerboard field number and the ratiometric position within a field
        ckbrd_field_ratio, ckbrd_field_number = np.modf((x - zero_position) / ckbrd_field_width)

        # calculate the alternating checkerboard "color" value 1 or -1
        ckbrd_field_sign = (((ckbrd_field_number.astype(np.int) & 1) << 1) - 1)

        # calculate the fill ratio of each element
        array_field_ratio = (array_field_width / 2.) / ckbrd_field_width
        X = np.clip(ckbrd_field_ratio, -array_field_ratio, array_field_ratio) / array_field_ratio

        return ckbrd_field_sign * X


    def checkerboard_edge_aware(self, x, dx, y, dy, wavelength, zero_position_x, zero_position_y):
        X = self.squarewave_edge_aware(x, dx, wavelength, zero_position_x)
        Y = self.squarewave_edge_aware(y, dy, wavelength, zero_position_y)
        return np.outer(Y, X)


    def mask_vertical_stripe(self, stripe_angle):
        """this creates a horizontal mask for a vertical stripe, where
        the edges are weighted correctly

        But note:
        There is a weird corner case where this approach is wrong, it occurs when
        the figure edge and one or more checkerboard edges are in one array field.

        TODO:
        consider background wavelength, stripe wavelength and both positions!
        """
        X = np.abs(wrap_to(self.phi - stripe_angle, -np.pi, np.pi))
        # = 0 @ stripewidth/2 + self.dphi/2.
        X0 = ((self.stripe_width + self.dphi)/2.)
        # make =1 @ stripewidth/2 - self.dphi/2.
        return np.clip((X - X0) / (-self.dphi), 0, 1)


    def calculate_luminance(self, stripe_angle, ground_angle, t=None):
        if t is not None and not self.t_start <= t < self.t_stop:
            return np.zeros((self.theta.shape[0], self.phi.shape[0]), dtype=np.float64)

        ground = self.checkerboard_edge_aware(self.phi, self.dphi, self.theta, self.dphi,
                                              self.ground_spatial_wavelength, ground_angle, 0)
        stripe = self.checkerboard_edge_aware(self.phi, self.dphi, self.theta, self.dphi,
                                              self.stripe_spatial_wavelength, stripe_angle, 0)

        A = self.mask_vertical_stripe(stripe_angle)

        return A * stripe + (1 - A) * ground


class CheckerboardFigureGroundStimulus(object):

    def __init__(self, stripe_width_radians,
                       ground_spatial_wavelength,
                       stripe_spatial_wavelength=np.radians(40),
                       t_start=-np.inf,
                       t_stop=np.inf,
                       elevation_extent=np.radians(60),
                       azimuth_samples=512):

        self.bg_spatial_wavelength = float(ground_spatial_wavelength)
        self.bg_circular_frequency = 2.*np.pi / self.bg_spatial_wavelength
        assert self.bg_circular_frequency % 1 == 0, "BG wavelength!"
        self.fg_spatial_wavelength = float(stripe_spatial_wavelength)
        self.fg_circular_frequency = 2.*np.pi / self.fg_spatial_wavelength
        self.fg_width = float(stripe_width_radians)

        # stimulus start and end
        self.t_start = float(t_start)
        self.t_stop = float(t_stop)

        # phi and theta grid
        self.phi, self.dphi = np.linspace(-np.pi, np.pi, int(azimuth_samples), endpoint=False, retstep=True)
        _theta = np.arange(0., float(elevation_extent), self.dphi)
        self.theta = np.hstack((-1*_theta[:1:-1], _theta))

        self._bg = self._prepare_bg()
        self._fg, self._fg_az_mask = self._prepare_fg()

    def _prepare_bg(self):
        bg_az = scipy.signal.square(self.bg_circular_frequency * self.phi)
        bg_el = scipy.signal.square(self.bg_circular_frequency * self.theta)
        bg = np.outer( bg_el, bg_az )
        return np.clip(bg, 0, 1)

    def _prepare_fg(self):
        fg_az = scipy.signal.square(self.fg_circular_frequency * self.phi)
        fg_az_mask = np.abs(self.phi) >= self.fg_width/2.
        fg_az[fg_az_mask] = 0.
        fg_el = scipy.signal.square(self.fg_circular_frequency * self.theta)
        fg = np.outer( fg_el, fg_az )
        return np.clip(fg, 0, 1), fg_az_mask

    def calculate_luminance(self, fg_angle, bg_angle, t=None):
        if t is not None and not self.t_start <= t < self.t_stop:
            return np.zeros((self.theta.shape[0], self.phi.shape[0]), dtype=np.float64)

        bg_rolls = int(np.round(wrap_to(bg_angle, -np.pi, np.pi) / self.dphi))
        fg_rolls = int(np.round(wrap_to(fg_angle, -np.pi, np.pi) / self.dphi))

        result = ( np.roll(self._bg, bg_rolls, axis=1) * np.roll(self._fg_az_mask, fg_rolls)
                 + np.roll(self._fg, fg_rolls, axis=1) * ~np.roll(self._fg_az_mask, fg_rolls) )

        return result



class SinusoidalFigureGroundStimulus(CheckerboardFigureGroundStimulus):

    def __init__(self, *args, **kwargs):
        self.A = kwargs.pop('amplitude', 1.)
        self.C = kwargs.pop('offset', 0.)
        super(SinusoidalFigureGroundStimulus, self).__init__(*args, **kwargs)


    def _prepare_bg(self):
        return None

    def _prepare_fg(self):
        return None, None

    def calculate_luminance(self, fg_angle, bg_angle, t=None):
        if t is not None and not self.t_start <= t < self.t_stop:
            return np.zeros((self.theta.shape[0], self.phi.shape[0]), dtype=np.float64)

        A = self.A*np.sin(self.bg_circular_frequency*(self.phi - bg_angle)) + self.C
        B = np.ones_like(self.theta)

        return np.outer( B, A )

class VerticalSlit(SinusoidalFigureGroundStimulus):


    def _prepare_bg(self):
        return None

    def _prepare_fg(self):
        return None, None

    def calculate_luminance(self, fg_angle, bg_angle, t=None):
        if t is not None and not self.t_start <= t < self.t_stop:
            return np.zeros((self.theta.shape[0], self.phi.shape[0]), dtype=np.float64)

        A = self.A*np.sin(self.fg_circular_frequency*(self.phi - fg_angle)) + self.C
        maximum = np.radians(15)

        A[self.phi < (maximum - self.fg_width/2.)] = self.C
        A[self.phi > (maximum + self.fg_width/2.)] = self.C
        B = np.ones_like(self.theta)
        return np.outer( B, A )

class HorizontalSlit(SinusoidalFigureGroundStimulus):


    def _prepare_bg(self):
        return None

    def _prepare_fg(self):
        return None, None

    def calculate_luminance(self, fg_angle, bg_angle, t=None):
        if t is not None and not self.t_start <= t < self.t_stop:
            return np.zeros((self.theta.shape[0], self.phi.shape[0]), dtype=np.float64)

        A = self.A*np.sin(self.fg_circular_frequency*(self.phi - fg_angle)) + self.C
        B = np.ones_like(self.theta)

        arr = np.outer( B, A )

        arr[self.theta < -self.fg_width/2. , : ] = self.C
        arr[self.theta > self.fg_width/2. , : ] = self.C

        return arr



if __name__ == '__main__':

    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    stimulus = CheckerboardFGCorrect( np.radians(40), np.radians(60), np.radians(30), azimuth_samples=512 )
    extent = map(np.degrees, [stimulus.phi[0], stimulus.phi[-1], stimulus.theta[0], stimulus.theta[-1]])

    phi_bg = np.linspace(0, -2*np.pi, 128)
    phi_fg = np.linspace(0, 2*np.pi, 128)

    fig = plt.figure()
    plt.xlabel('azimuth (deg)')
    plt.ylabel('elevation (deg)')

    img = []
    for bg, fg in zip(phi_bg, phi_fg):
        stim = stimulus.calculate_luminance(fg, bg)
        img.append([plt.imshow(stim, cmap=plt.get_cmap('gray'), extent=extent)])

    ani = animation.ArtistAnimation(fig, img, interval=100, blit=True, repeat_delay=0)

    plt.show()


