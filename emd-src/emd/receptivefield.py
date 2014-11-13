
import numpy as np


class ReceptiveField(object):
    """ReceptiveField for the EMD model"""
    def __init__( self, phi_c, sigma_theta, sigma_phi_a, sigma_phi_b):
        self.phi_c = phi_c
        self.sigma_theta = sigma_theta
        self.sigma_phi_a = sigma_phi_a
        self.sigma_phi_b = sigma_phi_b

    def __call__(self, PHI, THETA):
        a = np.exp(-((1./self.sigma_theta * THETA)**2)) * np.exp(-(1./self.sigma_phi_a*(PHI - self.phi_c)**2))
        b = np.exp(-((1./self.sigma_theta * THETA)**2)) * np.exp(-(1./self.sigma_phi_b*(PHI - self.phi_c)**2))
        cond = PHI >= self.phi_c
        result = np.where(cond, a, b)
        return result


if __name__ == "__main__":

    import matplotlib.pyplot as plt

    receptive_field_left = ReceptiveField( phi_c=np.radians( 15), sigma_theta=np.radians(33), sigma_phi_a=np.radians(102), sigma_phi_b=np.radians( 45))
    receptive_field_right= ReceptiveField( phi_c=np.radians(-15), sigma_theta=np.radians(33), sigma_phi_a=np.radians( 45), sigma_phi_b=np.radians(102))

    dphi = np.radians(2)
    phi = np.arange(-np.pi, np.pi, dphi)
    theta = np.linspace(-25*dphi, 25*dphi, 51)
    PHI, THETA = np.meshgrid(phi, theta)

    extent = map(np.degrees, [phi[0], phi[-1], theta[0], theta[-1]])
    ax0 = plt.subplot2grid((2,1), (0,0))
    ax0.imshow(receptive_field_left(PHI, THETA), extent=extent, interpolation='nearest')
    ax0.set_xlabel('phi (deg)')
    ax0.set_ylabel('theta (deg)')

    ax1 = plt.subplot2grid((2,1), (1,0))
    ax1.imshow(receptive_field_right(PHI, THETA), extent=extent, interpolation='nearest')
    ax1.set_xlabel('phi (deg)')
    ax1.set_ylabel('theta (deg)')

    plt.show()

