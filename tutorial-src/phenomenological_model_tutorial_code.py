#!/usr/bin/env python
"""
# Supplemental Material Python Code S1 #

This code can run as a standalone Python script, but is originally intended to
be viewed and executed inside an IPython Notebook.

It is recommended to get the full and up-to-date program-code from:

http://strawlab.org/asymmetric-motion


## Additional Notes regarding the format ##

This file will be seperated into snippets enclosed by special comments:

    # <TAG>
    ...
    # </TAG>

These can be used to seperate the code into snippents which can be embedded in
an IPython notebook.


For saving all plots as svgs, run:

    python phenomenological_model_tutorial_code.py

For interactively displaying all plots, run:

    python phenomenological_model_tutorial_code.py --show

"""
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--show", action="store_true")

import matplotlib
# When running in a normal python interpreter without the --show cmd line
# option use the Agg backend.
if not parser.parse_args().show:
    matplotlib.use("Agg")
# No text as paths in svgs. Assume font installed.
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.size'] = 14

# This is some interface magic, so that this code can run in python, and the
# IPython notebook after it has been reduced to the snippets.
try:
    from IPython import get_ipython
    if get_ipython() is None:
        raise ImportError
except ImportError:
    def get_ipython():
        class FakeMagic(object):
            def magic(self, *x, **y):
                pass
        return FakeMagic()

#~~~~~~~~~~~~~~~~~~~~
# START ALL SNIPPETS
#~~~~~~~~~~~~~~~~~~~~

# <imports>
import itertools

get_ipython().magic("matplotlib inline")  # required for ipynb compability
import matplotlib.pyplot as plt
import matplotlib.colors
import mpl_toolkits.axes_grid1
import numpy as np
import scipy.integrate
import scipy.interpolate
import sympy
import sympy.abc
# </imports>

print "Generating Figure S2A"
# <figureS2A>
fig = plt.figure()
fig.suptitle("Preferred Direction Motion Response")

# Motion response data from Schnell et al. 2010 Figure 2B
SCHNELL_MR2BA = np.array([[0.113636, 0.541573],
                          [0.520833, 0.676404],
                          [0.757575, 0.721348],
                          [1.022727, 0.885393],
                          [1.212121, 0.883146],
                          [1.505681, 0.795505],
                          [1.998106, 0.730337],
                          [3.011363, 0.528089],
                          [5.009469, 0.337078]]).T
A_wavelength = 22
MR2BA_tf, MR2BA = SCHNELL_MR2BA
MR2BA_av, MR2BA_n = MR2BA_tf * A_wavelength, MR2BA / MR2BA.max()

SCHNELL_MR2BB = np.array([[0.104166, 0.447191],
                          [0.511363, 0.768539],
                          [0.748106, 0.764044],
                          [1.013257, 0.883146],
                          [1.193181, 0.802247],
                          [1.496212, 0.782022],
                          [2.007575, 0.678651],
                          [3.020833, 0.476404],
                          [5.018939, 0.377528]]).T
B_wavelength = 45
MR2BB_tf, MR2BB = SCHNELL_MR2BB
MR2BB_av, MR2BB_n = MR2BB_tf * B_wavelength, MR2BB / MR2BB.max()

# Plot the experimental data against angular velocity
ax_omega = plt.subplot2grid((1, 2), (0, 0))
ax_omega.set_xlabel("Angular Velocity (deg/s)")
ax_omega.set_ylabel("Response (a.u.)")
ax_omega.set_xlim(0, 250.)
ax_omega.set_ylim(0, 1.1)
ax_omega.plot(MR2BA_av, MR2BA_n, color="black", marker="s", lw=2, markersize=3.)
ax_omega.plot(MR2BB_av, MR2BB_n, color="grey", marker="o", lw=2, markersize=3.)


# Preferred direction motion response approximation
def MotionResponsePD(temporal_frequency, optimal_frequency):
    """motion response in preferred direction"""
    _tf = temporal_frequency / optimal_frequency
    return 2 * _tf / (1 + _tf**2)  # the factor 2 normalizes the reponse

MR_PARAMETERS = {"optimal_frequency": 1.0}

# Experimental data and model response against temporal frequency
ax_tfreq = plt.subplot2grid((1, 2), (0, 1))
ax_tfreq.set_xlabel("Temporal Frequency (Hz)")
ax_tfreq.set_xlim(0, 5.)
ax_tfreq.set_ylim(0, 1.1)
ax_tfreq.set_yticklabels([])
ax_tfreq.plot(MR2BA_tf, MR2BA_n, color="black", marker="s", lw=2, markersize=3.)
ax_tfreq.plot(MR2BB_tf, MR2BB_n, color="grey", marker="o", lw=2, markersize=3.)

# the analytical motion response
tf = np.linspace(0, 5)
ax_tfreq.plot(tf, MotionResponsePD(tf, **MR_PARAMETERS), color="red", lw=2)
# </figureS2A>
plt.savefig("figureS2A.svg")


print "Generating Figure S2B"
# <figureS2B>
#plt.figure(figsize=(10, 4.5))
plt.figure()
plt.title("Asymmetric Motion Response")
plt.xlabel("Temporal Frequency (Hz)")
plt.ylabel("Response (a.u.)")
plt.xlim(-5., 5.)
plt.ylim(-0.5, 1.1)


class motion_response_pattern_wavelength(object):
    def __init__(self, optimal_frequency, asymmetry,
                 pattern_wavelength, preferred_direction):
        assert preferred_direction in ["CW", "CCW"]
        if preferred_direction == "CW":
            self._sign = 1
        else:
            self._sign = -1
        self._fopt = optimal_frequency
        self._asym = asymmetry
        self._lambdap = pattern_wavelength

    def _mrpd(self, tf):
        return MotionResponsePD(tf, self._fopt)

    def __call__(self, omega):
        tf = self._sign * np.atleast_1d(omega) / self._lambdap
        sign = np.clip(np.sign(tf), -self._asym, 1)
        return sign * self._mrpd(np.abs(tf))

# Constant taken from Fig 1D
MR_PARAMETERS["asymmetry"] = 0.4
# Motion responses of both cells:
test_pattern_wavelength = np.pi / 4.

mr_l = motion_response_pattern_wavelength(preferred_direction="CW",
        pattern_wavelength=test_pattern_wavelength, **MR_PARAMETERS)
mr_r = motion_response_pattern_wavelength(preferred_direction="CCW",
        pattern_wavelength=test_pattern_wavelength, **MR_PARAMETERS)

# Plot motion response
omega = np.linspace(-5, 5, 1001)
plt.plot(omega, mr_l(omega), "r", lw=2)  # left cell
plt.plot(omega, mr_r(omega), "b", lw=2)  # right cell
# </figureS2B>
plt.savefig("figureS2B.svg")


print "Generating Figure S2C"
# <figureS2C>
plt.figure()
plt.title("Receptive Field")
plt.xlabel("Angular Position (deg)")
plt.ylabel("Response (a.u.)")
plt.xlim(-180., 180.)
plt.xticks([-180, -120, -60, 0, 60, 120, 180])
plt.ylim(0., 0.6)
plt.yticks([8 / (5 * np.pi)], [u"8/5\u03C0"])

# Receptive Field from Fig 3C
SCHNELL_RF3C2Y = np.array([0.12, 0.32, 0.38, 0.35, 0.28, 0.32, 0.33, 0.38,
                           0.40, 0.44, 0.50, 0.63, 0.72, 0.76, 0.83, 0.84,
                           0.87, 0.95, 0.95, 0.92, 0.90, 0.90, 0.86, 0.83,
                           0.72])
SCHNELL_RF3C2X = np.linspace(-42, 94, 25)

RF_X = SCHNELL_RF3C2X
RF_Y = SCHNELL_RF3C2Y / SCHNELL_RF3C2Y.max() * 8 / (5 * np.pi)

RF_PARAMETERS = {"phi_max": np.radians(60.)}
phi = np.linspace(-np.pi, np.pi, 101)


def receptive_field(phi, phi_max):
    return 8 / (5 * np.pi) * np.cos((phi - phi_max) / 2.) ** 6

# Plot experimental data
plt.plot(RF_X, RF_Y, color="black", marker="s", ls="", markersize=3.)
plt.plot(np.degrees(phi), receptive_field(phi, **RF_PARAMETERS), "r-", lw=2)
# </figureS2C>
plt.savefig("figureS2C.svg")


print "Generating Figure S2D"
# <figureS2D>
plt.figure()
plt.title("Receptive Fields")
plt.xlabel("Angular Position (deg)")
plt.ylabel("Response (a.u.)")
plt.xlim(-180., 180.)
plt.xticks([-180, -120, -60, 0, 60, 120, 180])
plt.ylim(0., 0.6)
plt.yticks([8 / (5 * np.pi)], [u"8/5\u03C0"])
# left cell
plt.plot(np.degrees(phi), receptive_field(phi, **RF_PARAMETERS), "r-", lw=2)
# right cell
plt.plot(-np.degrees(phi), receptive_field(phi, **RF_PARAMETERS), "b-", lw=2)
# </figureS2D>
plt.savefig("figureS2D.svg")


print "Generating Figure S2F"
# <figureS2F>
plt.figure()

class receptive_field_figure_width(object):
    def __init__(self, figure_width, phi_max, position):
        assert 0 <= figure_width <= 2 * np.pi
        assert -np.pi <= phi_max <= np.pi
        assert position in ["left", "right"]
        self._fwhalf = figure_width / 2.
        self._phi_max = abs(phi_max) if position is "left" else -abs(phi_max)
        _root = sympy.integrate(
                    8 / (5 * sympy.pi) * sympy.cos(sympy.abc.x / 2) ** 6,
                    sympy.abc.x)
        self._rootfunc = sympy.lambdify(sympy.abc.x, _root, np)

    def __call__(self, phi, background=False):
        phi = phi - self._phi_max
        if not background:
            return (self._rootfunc(phi + self._fwhalf)
                   - self._rootfunc(phi - self._fwhalf))
        else:
            return (self._rootfunc(phi - self._fwhalf + 2 * np.pi)
                   - self._rootfunc(phi + self._fwhalf))


class cell_response(object):
    def __init__(self, position, figure_width,
                 pattern_wavelength_fg, pattern_wavelength_bg=np.inf):
        assert position in ["left", "right"]
        if position == "left":
            PD = "CW"
        else:
            PD = "CCW"
        self._mrf = motion_response_pattern_wavelength(preferred_direction=PD,
                    pattern_wavelength=pattern_wavelength_fg, **MR_PARAMETERS)
        self._mrg = motion_response_pattern_wavelength(preferred_direction=PD,
                    pattern_wavelength=pattern_wavelength_bg, **MR_PARAMETERS)
        self._rf = receptive_field_figure_width(position=position,
                                figure_width=figure_width, **RF_PARAMETERS)

    def __call__(self, phi, omega_F, omega_G=0):
        return (self._rf(phi) * self._mrf(omega_F) +
                self._rf(phi, background=True) * self._mrg(omega_G))


class torque_response(object):
    def __init__(self, figure_width, pattern_wavelength_fg,
                       pattern_wavelength_bg=np.inf, ct=1.0):
        self._hsl = cell_response("left", figure_width=figure_width,
                            pattern_wavelength_fg=pattern_wavelength_fg,
                            pattern_wavelength_bg=pattern_wavelength_bg)
        self._hsr = cell_response("right", figure_width=figure_width,
                            pattern_wavelength_fg=pattern_wavelength_fg,
                            pattern_wavelength_bg=pattern_wavelength_bg)
        self._ct = ct

    def __call__(self, phi, omega_F, omega_G=0):
        return self._ct * (self._hsl(phi, omega_F, omega_G)
                         - self._hsr(phi, omega_F, omega_G))

phi, tf_sf = np.meshgrid(np.linspace(-np.pi, np.pi, 101),
             np.radians(22.5)*np.linspace(-6, 6, 101))
extents1f = [-180, 180, -6, 6]

torque = torque_response(np.radians(22.5), np.radians(22.5))

response = torque(phi, tf_sf)

plt.title("Predicted Small-Field Behavioural Response")
plt.xlabel("Angular Position (deg)")
plt.ylabel("Temporal Frequency (Hz)")
im = plt.imshow(response, origin="lower", cmap=plt.get_cmap("RdBu_r"),
           extent=extents1f, aspect='auto')
plt.colorbar(label="Response (a.u.)")
# </figureS2F>
plt.savefig("figureS2F.svg")


print "Generating Figure S2E"
# <figureS2E>
fig = plt.figure()
grid = mpl_toolkits.axes_grid1.ImageGrid(fig, 111, (5, 2), axes_pad=0.1)

phi, omega = np.meshgrid(np.linspace(-np.pi, np.pi, 101),
                         np.linspace(-np.pi, np.pi, 101))

lambdap = np.radians(22.5)
widths = np.radians([360, 270, 180, 90, 45])
positions = ["left", "right"]

griddata = {}
for i, (w, p) in enumerate(itertools.product(widths, positions)):
    HS = cell_response(p, figure_width=w, pattern_wavelength_fg=lambdap)
    DATA = HS(phi, omega)
    griddata[i] = DATA

VMAX, VMIN = 1, -1

for i, data in griddata.items():
    grid[i].imshow(data, origin="lower", vmin=VMIN, vmax=VMAX,
                   cmap=plt.get_cmap("RdBu_r"))
    grid[i].set_xticks([1, 51, 101])
    grid[i].set_xticklabels([-180, 0, 180 if i == 9 else ""])
    grid[i].set_yticks([1, 51, 101])
    grid[i].set_yticklabels([-8, 0, 8 if i < 2 else ""])

grid[4].set_ylabel("Temporal Frequency (Hz)")
grid[8].set_xlabel("Angular Position (deg)", x=1.)
grid[0].set_title("left\ncell")
grid[1].set_title("right\ncell")
grid[1].text(115, 110, "figure\nwidth")
grid[1].text(110, 40, "360 deg")
grid[3].text(110, 40, "270 deg")
grid[5].text(110, 40, "180 deg")
grid[7].text(110, 40, "90 deg")
grid[9].text(110, 40, "45 deg")
# </figureS2E>
plt.savefig("figureS2E.svg")


print "Generating Figure 3A top"
# <figure3Atop>
plt.figure()

FIGURE_WIDTH_DEG = 22.5
BACKGROUND_WAVELENGTH_DEG = 22.5
M0 = 0.1


def figure_potential(phi, omega_bg, figure_width, pwlbg, m0):

    rfl = receptive_field_figure_width(position="left",
                        figure_width=figure_width, **RF_PARAMETERS)
    rfr = receptive_field_figure_width(position="right",
                        figure_width=figure_width, **RF_PARAMETERS)
    bgl = cell_response("left", figure_width=figure_width,
                                pattern_wavelength_fg=np.inf,
                                pattern_wavelength_bg=pwlbg)
    bgr = cell_response("right", figure_width=figure_width,
                                 pattern_wavelength_fg=np.inf,
                                 pattern_wavelength_bg=pwlbg)

    return ((rfl(phi) - rfr(phi)) * m0 +
            (bgl(phi, 0, omega_G=omega_bg) - bgr(phi, 0, omega_G=omega_bg)))

force = lambda phi, rho: -figure_potential(phi, rho,
                            np.radians(FIGURE_WIDTH_DEG),
                            np.radians(BACKGROUND_WAVELENGTH_DEG), M0)

rho_deg = np.linspace(0, -0.5, 71)  # The upper bound should be > 2*OMEGAC
phi_deg = np.linspace(0, 180, 51)
rho_rad = np.radians(rho_deg)
phi_rad = np.radians(phi_deg)
RHO_RAD, PHI_RAD = np.meshgrid(rho_rad, phi_rad)
RHO_DEG, PHI_DEG = np.meshgrid(rho_deg, phi_deg)
Z = force(PHI_RAD, RHO_RAD)
extent = (rho_deg[0], rho_deg[-1], phi_deg[-1], phi_deg[0])
clim = 0.5 * max(abs(Z.max()), abs(Z.min()))

plt.imshow(Z, aspect='auto', extent=extent,
              clim=(-clim, clim), cmap=plt.get_cmap("PiYG"))
im = plt.colorbar(label='Torque (a.u.)')
quadcontour = plt.contour(RHO_DEG, PHI_DEG, Z, colors='k',
                          linestyles='dashed', linewidths=2, levels=[0.0])

OMEGAC = None
# color stable branches in green
for linecollection in quadcontour.collections:
    lines = linecollection.get_segments()
    for line in lines:
        brho, bphi = line.T
        bx, by = np.radians(line.T)
        OMEGAC = np.degrees(bx.min())
        m2 = stability = (force(by - 0.01, bx) - force(by + 0.01, bx)) >= 0
        plt.plot(brho[m2], bphi[m2], color='k', ls='-', lw=2)

plt.xlim(rho_deg[0], 2 * OMEGAC)
plt.xticks([0, OMEGAC, 2 * OMEGAC], ["0", "$\\omega_{c}$", "2$\\omega_{c}$"])
plt.xlabel("Background Velocity")
plt.ylabel("Figure Position (deg)")
plt.ylim(phi_deg[0], phi_deg[-1])
# </figure3Atop>
plt.savefig("figure3Atop.svg")


print "Generating Figure S2G / 3A bottom"
# <figureS2G-3Abottom>
FIG_WIDTH_RAD = np.radians(FIGURE_WIDTH_DEG)
GND_WL_RAD = np.radians(BACKGROUND_WAVELENGTH_DEG)


def NEWpotential(figure_width, m0):
    N = 128
    rfl = receptive_field_figure_width(position="left",
                        figure_width=figure_width, **RF_PARAMETERS)
    rfr = receptive_field_figure_width(position="right",
                        figure_width=figure_width, **RF_PARAMETERS)

    def _integral_part(phi):
        _rf = lambda x: (rfl(x) - rfr(x)) * m0
        return scipy.integrate.fixed_quad(_rf, 0., phi)[0]

    DATA = []
    for pos in np.linspace(0, 2 * np.pi, N, endpoint=False):
        DATA.append(_integral_part(pos))
    res = np.array(DATA + DATA + DATA)
    phi = np.linspace(-2 * np.pi, 4 * np.pi, 3 * N, endpoint=False)
    func = scipy.interpolate.interp1d(phi, res)

    def _potential(phi):
        arg = (phi + np.pi) % (2 * np.pi) - np.pi
        return func(arg)

    return _potential


def LINpotential(omega_bg, figure_width, pwlbg):

    N = 256
    bgl = cell_response("left", figure_width=figure_width,
                                pattern_wavelength_fg=np.inf,
                                pattern_wavelength_bg=pwlbg)
    bgr = cell_response("right", figure_width=figure_width,
                                 pattern_wavelength_fg=np.inf,
                                 pattern_wavelength_bg=pwlbg)

    def Vstar(phi):
        _cmr = lambda x: -(bgl(x, 0, omega_G=omega_bg)
                           - bgr(x, 0, omega_G=omega_bg))
        return scipy.integrate.fixed_quad(_cmr, 0, phi)[0]

    DATA = []
    phis = np.linspace(-2 * np.pi, 2 * np.pi, N + 1, endpoint=True)
    for phi in phis:
        DATA.append(Vstar(phi))

    res = np.array(DATA)
    func = scipy.interpolate.interp1d(phis, res)
    return func


def Ustar(m0, omega_bg=0.0, fw=FIG_WIDTH_RAD, bwl=GND_WL_RAD):

    func0 = NEWpotential(fw, m0)
    func1 = LINpotential(omega_bg, fw, bwl)

    def _potential(x):
        ret = func1(x) - func0(x)
        if np.isscalar(ret):
            try:
                return ret.asscalar()
            except AttributeError:
                return ret
        else:
            return ret

    return _potential


def probability(omega_bg, m0, d, c):

    N = 91
    phis, dphi = np.linspace(0, 2 * np.pi, N, endpoint=False, retstep=True)

    eta = d / float(c)
    ustar = Ustar(m0, omega_bg)

    _f1 = lambda x: np.exp(-eta * ustar(x))

    P = []
    for phi in phis:
        f0 = np.exp(eta * ustar(phi))
        f1, _ = scipy.integrate.quad(_f1, phi - 2 * np.pi, phi,
                                     epsabs=1e-3, epsrel=1e-3)
        P.append(f0 * f1)

    arr = np.array(P)
    arr = arr / arr.sum()
    return arr, ustar(phis)


plt.figure()
ax0 = plt.subplot2grid((1, 2), (0, 0))
ax1 = plt.subplot2grid((1, 2), (0, 1))

ax0.set_xlabel("Figure Position (deg)")
ax0.set_ylabel("Probability (%/deg)")

DATA = []

rho_rad = rho_rad[::4]
M = 1. * len(rho_rad)

LC = matplotlib.colors.LinearSegmentedColormap.from_list(
                                      "FigureGroundPaper",
                                      [(0x05/255., 0x71/255., 0xB0/255., 1.0),
                                      (0xCA/255., 0x00/255., 0x20/255., 1.0)],
                                      N=len(rho_rad))

for i, rho in enumerate(rho_rad):
    PROB, POT = probability(rho, M0, 4., 0.01)
    N = len(PROB)/2
    PROB = PROB/np.sum(PROB) * len(PROB) / 360. * 100.
    ax0.plot(np.linspace(-180, 180, 91, endpoint=False),
                         np.roll(PROB, N), color=LC(i/M))
    DATA.append(np.roll(PROB, N))

ax0.set_ylim(bottom=0)
ax0.set_xlim(-180, 180)

arr = np.vstack(DATA).T[::-1, :]
ax1.imshow(arr,
           cmap=plt.get_cmap("Greys"), aspect='auto',
           interpolation="None", extent=[np.degrees(rho_rad[0]),
           np.degrees(rho_rad[-1]), -180, 180])

# color stable branches in green
for linecollection in quadcontour.collections:
    lines = linecollection.get_segments()
    for line in lines:
        brho, bphi = line.T
        bx, by = np.radians(line.T)
        m2 = stability = (force(by - 0.01, bx) - force(by + 0.01, bx)) <= 0
        ax1.plot(brho[m2], bphi[m2], color='w', ls='--', lw=4)
        ax1.plot(brho[~m2], bphi[~m2], color='w', ls='-', lw=4)
        ax1.plot(brho[m2], bphi[m2], color='k', ls='--', lw=2)
        ax1.plot(brho[~m2], bphi[~m2], color='k', ls='-', lw=2)

ax1.set_ylim(-20, 80)
ax1.set_xlim(rho_deg[0], 2*OMEGAC)
ax1.set_xticks([0, OMEGAC, 2*OMEGAC])
ax1.set_xticklabels(["0", "$\\omega_{c}$", "2$\\omega_{c}$"])
ax1.set_xlabel("Background Velocity")
ax1.set_ylabel("Figure Position (deg)", labelpad=-8)
# </figureS2G-3Abottom>
plt.savefig("figureS2G-3Abottom.svg")

plt.show()

