---
title: Supplemental Experimental Procedures
subtitle: Asymmetric processing of visual motion for simultaneous object and background responses
author:
- name: Lisa Fenk
  affiliation: IMP
- name: Andreas Poehlmann
  affiliation: IMP
- name: Andrew Straw
  affiliation: IMP
...

# The phenomenological model #

Here we describe a phenomenological model of fly visual behaviors based on
a few well-known properties of wide-field motion integrators such as the lobula
plate tangential cells (LPTCs) of flies. A primary goal
was to formulate a very simple model to facilitate rigorous
mathematical analysis. It predicts turning responses of the
fly to piecewise-defined open- and closed-loop stimuli. The model abstracts fly
behavior to a one degree of freedom system (azimuth) and reduces stimuli to
horizontal piecewise-defined patterns at zero elevation. We provide
a theoretical framework to analytically predict turning responses here, as well
as a software implementation in Supplemental Code S1 (also available at [http://strawlab.org/asymmetric-motion/](http://strawlab.org/asymmetric-motion/)). All retinal positions are
specified as azimuthal angle from the visual midline (right handed coordinate
system).

<!--@ .INCLUDECODESNIPPET python phenomenological_model_tutorial_code.py imports @-->

The internal parameters of the model are, with a single exception,
based exclusively on measurements from electrophysiological
experiments of well investigated LPTCs, _Drosophila_ Horizontal System (HS) cells [Schnell
et al. (2010)](http://dx.doi.org/10.1152/jn.00950.2009). The
remaining free parameter, $C_\text{T}$, is used to account for experimental coupling conditions and
predict behavioral responses under a variety of those conditions. The parameter is a constant specifying the cell output to torque coupling,
and we discuss it below.

## The motion response ##

A fundamental property of an HS-cell is its change in membrane
potential when a wide-field stimulus, such as horizontal sinusoidal
grating covering the whole visual field, is rotated azimuthally around
the fly. We call the function that describes the dependence of
membrane potential responses to movement of the visual
stimulus the motion response. Its dependent variable is the temporal
frequency, $f_{\text{t}}$, which is described by $f_{\text{t}} = \omega / \lambda_{\text{p}}$ given a 
stimulus with a pattern wavelength, $\lambda_{\text{p}}$, and
an angular velocity, $\omega$. We use the measured relative motion
response from Fig. 2B of [Schnell et
al. (2010)](http://dx.doi.org/10.1152/jn.00950.2009), to fit
the parameters of our function. The mathematical form of our function
is given by the Reichardt correlator as in [Borst, Bahde (1986)](http://dx.doi.org/10.1007/BF00363978) after
[Reichardt, Varju (1959)](http://www.ncbi.nlm.nih.gov/pubmed/13854486). Like most LPTCs, HS-cells respond with
depolarizing and hyperpolarizing responses in the so-called
preferred and null directions, respectively. The relative motion
response to preferred direction stimuli with constant pattern wavelength is:

<!-- As in other fly species, _Drosophila_  -->

\begin{equation}
M_{\text{PD}}(f_{\text{t}}) =  \frac{
    \tau_\text{opt} f_{\text{t}}
    }{
    1 + \left( \tau_\text{opt} f_{\text{t}}\right)^2
    }
\end{equation}

Adjusting the parameter, $f_{\text{opt}} = \tau_{\text{opt}}^{-1}$, to
$1\,\mathrm{Hz}$, approximates the data reasonably well<!--@
.OPTIONALTEXT " (Supplemental Fig. S2A)" @-->.

<!--@ .INCLUDECODESNIPPET python phenomenological_model_tutorial_code.py figureS2A @-->
<!--@ .OPTIONALFIGURE figureS2A "Experimentally measured motion response in the preferred direction by [Schnell et al. 2010](http://dx.doi.org/10.1152/jn.00950.2009) for angular velocity (left) and temporal frequency (right). Model motion response (red line) was manually fit to the temporal frequency data." @-->

A critical feature of our model is the asymmetry of the response to
the direction of motion of the stimulus. In particular, the motion
response to null direction stimuli is not only inverted in sign but
also reduced in amplitude. Basing our estimate on the values published in Fig. 1D of [Schnell et
al. (2010)](http://dx.doi.org/10.1152/jn.00950.2009), we
account for this by a factor $C_{\text{ND}} = 0.4$. The
relative motion response for stimuli moving in both directions can be
written as:

$$M(f_{\text{t}}) = \left\{ \begin{array}{rl}
M_{\text{PD}}(f_{\text{t}}) &\mbox{ if $f_{\text{t}}
\ge 0$} \\ -C_{\text{ND}} M_{\text{PD}}(-f_{\text{t}})
&\mbox{ otherwise} \end{array} \right.$$


The normalized response to a moving sinusoidal pattern in
both the null (amplitude is _negative_) and the preferred (amplitude is _positive_) direction in
temporal frequency space, $f_{\text{t}} = \omega/\lambda_{\text{p}}$,
for a pair of neurons with opposite direction selectivities
(e.g. contralateral cells)<!--@ .OPTIONALTEXT " is shown in Supplemental Fig. S2B." @-->

<!--@ .INCLUDECODESNIPPET python phenomenological_model_tutorial_code.py figureS2B @-->
<!--@ .OPTIONALFIGURE figureS2B "Motion response in preferred and null direction for left (red) and right (blue) model cells." @-->

## The receptive field ##

For a small stimulus, the response of an LPTC is also dependent on
retinal position.  We use the local motion sensitivity measured
at zero degrees elevation from Fig. 3C of [Schnell et
al. (2010)](http://dx.doi.org/10.1152/jn.00950.2009) to
define the spatial receptive field for our model. Because we simulated visual
stimuli along the unit circle, the receptive field is required to be
$2\pi$ periodic, and we therefore choose a trigonometric function.
However, its exact nature is not important and was chosen for
convenience. Our model receptive field <!--@ .OPTIONALTEXT " (see
Supplemental Fig. S2C)" @--> is given as:

\begin{equation}
R(\varphi) = \frac{8}{5\pi}\cos{\left(\frac{1}{2}(\varphi \pm \varphi_\text{max})\right)}^{6}
\end{equation}

In this model, $\varphi_{\text{max}}$ is the azimuthal angle at which the model cell
responds most strongly to a small field stimulus and the factor $\frac{8}{5\pi}$
normalizes the area under the curve to 1.

<!--@ .INCLUDECODESNIPPET python phenomenological_model_tutorial_code.py figureS2C @-->
<!--@ .OPTIONALFIGURE figureS2C "Model receptive field was manually fit (red line) to experimental data (black points) from [Schnell et al. 2010](http://dx.doi.org/10.1152/jn.00950.2009)." @-->

The receptive fields of the left (red) and right (blue) model cells
are symmetric around zero<!--@ .OPTIONALTEXT ", as shown in Supplemental Fig. S2D" @-->.

<!--@ .INCLUDECODESNIPPET python phenomenological_model_tutorial_code.py figureS2D @-->
<!--@ .OPTIONALFIGURE figureS2D "Receptive fields for left (red) and right (blue) model cells." @-->


## The cell's response ##

We model the cell's complete output as the integral over all angular
positions weighted with the receptive field:

\begin{equation}
W(f_{\text{t}}) = \int_{-\pi}^{\pi} R(\varphi) M(f_{\text{t}}) \,\mathrm{d}\varphi
\end{equation}

For clarity, we expand this to explicitly show all parameters
(enclosed in square brackets):

$$W(f_{\text{t}}, [\varphi_{\text{max}},
                               f_{\text{opt}},
                               C_{\text{ND}}]) =
        \int_{-\pi}^{\pi} R(\varphi, [\varphi_{\text{max}}])
                          M(f_{\text{t}},[
                               f_{\text{opt}},
                               C_{\text{ND}}]) \,\mathrm{d}\varphi
$$

As described above, each of these three parameters, $\varphi_{\text{max}},
f_{\text{opt}}, C_{\text{ND}}$, was adjusted by eye to the experimental
data of [Schnell et
al. (2010)](http://dx.doi.org/10.1152/jn.00950.2009). This
model describes the noise-free normalized graded potential of a
wide-field motion integrator cell. We do not adjust these parameters for
any predictions that we make later.


## Introducing the figure-ground stimulus #

So far, the model did not consider responses to spatially complex
stimuli. The motion response (equation 1) was defined for a panoramic
sinusoidal grating, which for a certain angular velocity, $\omega$, and
pattern wavelength, $\lambda_{\text{p}}$, generates a temporal frequency,
$f_{\text{t}}=\omega / \lambda_{\text{p}}$. The receptive field
(equation 2) was defined for a small stimulus. We now introduce a figure-ground
stimulus, where the background pattern is occluded by the figure, with the
following parameters:

  * figure width $\beta$
  * figure position $\varphi_{\text{F}}$
  * figure pattern wavelength ${\lambda_{\text{F}}}$
  * figure velocity $\omega_{\text{F}}$
  * ground pattern wavelength $\lambda_{\text{G}}$
  * ground velocity $\omega_{\text{G}}$

All the parameters above are inherent to the stimulus, not to the
model. This is important because the responses to the figure or to the
background --- which we derive below --- evolve from these parameters,
and are not properties of the model cell. Please note that in this
text, we use the word "figure" to describe a small-field visual
object, even when there is no background contrast.

To use these stimulus parameters in the model, we define $\omega_{\text{FG}}$
and $\lambda_{\text{FG}}$ dependent on the angular position, $\varphi$, the
figure position, $\varphi_{\text{F}}$, and width, $\beta$, for the figure-ground
stimulus.

$$
\lambda_{\text{FG}}(\varphi, [\varphi_{\text{F}}, \beta]) =  \left\{ \begin{array}{rl}
                    \lambda_{\text{F}} & \mbox{ if $\varphi_{\text{F}}-\frac{\beta}{2} \le \varphi \le \varphi_{\text{F}}+\frac{\beta}{2}$ } \\
                    \lambda_{\text{G}} & \mbox{ otherwise } \\
                       \end{array} \right.
$$

$$
\omega_{\text{FG}}(\varphi, [\varphi_{\text{F}}, \beta]) =  \left\{ \begin{array}{rl}
                    \omega_{\text{F}} & \mbox{ if $\varphi_{\text{F}}-\frac{\beta}{2} \le \varphi \le \varphi_{\text{F}}+\frac{\beta}{2}$ } \\
                    \omega_{\text{G}} & \mbox{ otherwise } \\
                       \end{array} \right.
$$

With this we can rewrite the cell response as:

$$
\begin{split}
W(f_{\text{t}}) &= W(\omega_{\text{FG}} / \lambda_{\text{FG}}) =
        \int_{-\pi}^{\pi} R(\varphi) M(\omega_{\text{FG}} / \lambda_{\text{FG}})\,\mathrm{d}\varphi \\
 &= \int_{-\frac{1}{2}\beta}^{+\frac{1}{2}\beta} R(\varphi_{\text{F}} + \varphi)M(\omega_{\text{F}} / \lambda_{\text{F}})\,\mathrm{d}\varphi
  + \int_{+\frac{1}{2}\beta}^{2\pi-\frac{1}{2}\beta} R(\varphi_{\text{F}} + \varphi) M(\omega_{\text{G}} / \lambda_{\text{G}})\, \mathrm{d}\varphi
\end{split}
$$

We can interpret these two summands as a response to the figure, and a response
to the background. But, again, they originate from the properties of the stimulus, not from
the properties of the wide-field motion integrator. We define the figure and ground responses
as:

$$
W_{\text{F}}(\varphi_{\text{F}}, \omega_{\text{F}}, [\lambda_{\text{F}}, \beta]) = \int_{-\frac{1}{2}\beta}^{+\frac{1}{2}\beta} R(\varphi_{\text{F}} + \varphi)M(\omega_{\text{F}} / \lambda_{\text{F}})\,\mathrm{d}\varphi
$$

$$
W_{\text{G}}(\varphi_{\text{F}}, \omega_{\text{G}}, [\lambda_{\text{G}}, \beta]) = \int_{+\frac{1}{2}\beta}^{2\pi-\frac{1}{2}\beta} R(\varphi_{\text{F}} + \varphi) M(\omega_{\text{G}} / \lambda_{\text{G}})\, \mathrm{d}\varphi
$$

Since the basic property of our figure-ground stimulus is that the
figure occludes the background, the background response, $W_{\text{G}}$,
is of course dependent on the figure position, $\varphi_{\text{F}}$, as
a variable and the figure width, $\beta$, as a parameter.

If we use a figure of width, $\beta$, and pattern wavelength, $\lambda_{\text{F}}$,
in front of a contrast-free background ($\lambda_{\text{G}} = \infty$), we can
plot the response of the cell in phase space<!--@ .OPTIONALTEXT " (see Supplemental Fig. S2F)" @-->.

<!--@ .INCLUDECODESNIPPET python phenomenological_model_tutorial_code.py figureS2F @-->
<!--@ .OPTIONALFIGURE figureS2F "Torque resulting from the model in response to a (22.5 deg width and wavelength) figure stimulus." @-->

## From cell output to torque #

We model the torque output of the fly to be proportional to the
difference of the left and the right model cell responses. We introduce a
factor, $C_{\text{T}}$, which scales the normalized output to torque. This is the only
free parameter of our model.

$$
T(f_{\text{t}}) = C_{\text{T}} \Big[ W^{\text{(left)}}(f_{\text{t}}) - W^{\text{(right)}}(f_{\text{t}})  \Big]
$$

To simulate a noisy cell response, we add noise to the torque response
of our system at this point to get:

$$
T^{*}(f_{\text{t}}, t) = T(f_{\text{t}}) + N(t) = C_{\text{T}} \Big[ W^{\text{(left)}}(f_{\text{t}}) - W^{\text{(right)}}(f_{\text{t}})  \Big] + N(t)
$$

We rewrite this torque response for a figure-ground stimulus, using
the following abbreviations:

$$
\begin{split}
R_{\text{F}}^{(\text{left/right})}(\varphi_{\text{F}}) &= \int_{-\frac{1}{2}\beta}^{+\frac{1}{2}\beta} R^{(\text{left/right})}(\varphi_{\text{F}} + \varphi)\,\mathrm{d}\varphi \\
R_{\text{G}}^{(\text{left/right})}(\varphi_{\text{F}}) &= \int_{+\frac{1}{2}\beta}^{2\pi - \frac{1}{2}\beta} R^{(\text{left/right})}(\varphi_{\text{F}} + \varphi)\,\mathrm{d}\varphi \\
X^{(-)} &= X^{\text{(left)}} - X^{\text{(right)}} \\
X^{(+)} &= X^{\text{(left)}} + X^{\text{(right)}}
\end{split}
$$

Putting this and the figure-ground stimulus parameters into $T(f_{\text{t}})$ gives:

$$
\begin{split}
T(f_{\text{t}}) &= T(\varphi_{\text{F}}, \omega_{\text{F}}, \omega_{\text{G}}, [\lambda_{\text{F}}, \lambda_{\text{G}}, \beta]) = C_{\text{T}} W^{(-)}(\varphi_{\text{F}}, \omega_{\text{F}}, \omega_{\text{G}}, [\lambda_{\text{F}}, \lambda_{\text{G}}, \beta]) \\
 &= C_{\text{T}} \Big( W^{(-)}_{\text{F}}(\varphi_{\text{F}}, \omega_{\text{F}}, [\lambda_{\text{F}}, \beta]) + W^{(-)}_{\text{G}}(\varphi_{\text{F}}, \omega_{\text{G}}, [\lambda_{\text{G}}, \beta]) \Big) \\
\end{split}
$$

For the sake of simplicity, we set the
background contrast to zero for now ($\lambda_{\text{G}} = \infty
\rightarrow W_{\text{G}}^{(-)} = 0$).

$$
T(\varphi_{\text{F}}, \omega_{\text{F}}) = \frac{C_{\text{T}}}{2} \Big[
    R_{\text{F}}^{(-)}(\varphi_{\text{F}}) M^{(+)}(\omega_{\text{F}} / \lambda_{\text{F}})
  + R_{\text{F}}^{(+)}(\varphi_{\text{F}}) M^{(-)}(\omega_{\text{F}} / \lambda_{\text{F}})
  \Big]
$$

This gives the response properties of a figure in front of a homogeneous background.
For different figure widths, the phase space representation of this output is<!--@ .OPTIONALTEXT " shown in Supplemental Fig. S2E" @-->. 
The receptive fields of both model cells make the initially homogeneous and optomotor like torque response ($360\,^{\circ}$) more and more position specific with decreasing figure width.

<!--@ .INCLUDECODESNIPPET python phenomenological_model_tutorial_code.py figureS2E @-->
<!--@ .OPTIONALFIGURE figureS2E "Responses of left and right cells (left and right panels, respectively) to small field stimuli of varying width $\beta$." @-->


Note: The mathematical conversion on the figure terms can be applied to the ground
terms in the same manner. For the figure-ground stimulus with $\lambda_{\text{G}} < \infty$ we would get:

$$
\begin{split}
T(\varphi_{\text{F}}, \omega_{\text{F}}, \omega_{\text{G}}) &= \frac{C_{\text{T}}}{2} \Big[
    R_{\text{F}}^{(-)}(\varphi_{\text{F}}) M^{(+)}(\omega_{\text{F}} / \lambda_{\text{F}})
  + R_{\text{F}}^{(+)}(\varphi_{\text{F}}) M^{(-)}(\omega_{\text{F}} / \lambda_{\text{F}})
  \Big] \\
 &+ \frac{C_{\text{T}}}{2} \Big[
    R_{\text{G}}^{(-)}(\varphi_{\text{F}}) M^{(+)}(\omega_{\text{G}} / \lambda_{\text{G}})
  + R_{\text{G}}^{(+)}(\varphi_{\text{F}}) M^{(-)}(\omega_{\text{G}} / \lambda_{\text{G}})
  \Big]
\end{split}
$$

## The dynamic equation #

Above we derived an analytical model for the torque output of a fly based on
the fundamental properties of motion integrator cells. In the following
sections, we use this model to predict fly behavior during different
figure-ground discrimination tasks, where two coupling conditions will be important.
These are
closed-loop and open-loop. For an open-loop stimulus, the torque
 produced does not influence the stimulus. The stimulus moves along
 a predefined trajectory in space time, and we measure torque output over time.
In a closed-loop scenario, we apply the torque produced by the model on
the part of the stimulus in closed-loop. For example, for a closed-loop
figure (position: $\varphi_{\text{F}}$, velocity: $\dot{\varphi}_{\text{F}}$)
with no or stationary background, we obtain a Langevin
type dynamic equation of the form:

$$\begin{split}
\Theta \ddot{\varphi}_{\text{F}} &= T^{*}(\varphi_{\text{F}}, \dot{\varphi}_{\text{F}}) \\
 &= \frac{C_{\text{T}}}{2} \Big[
    R_{\text{F}}^{(-)}(\varphi_{\text{F}}) M^{(+)}(\dot{\varphi}_{\text{F}} / \lambda_{\text{F}})
  + R_{\text{F}}^{(+)}(\varphi_{\text{F}}) M^{(-)}(\dot{\varphi}_{\text{F}} / \lambda_{\text{F}})
  \Big] + N(t)
\end{split}
$$

Here we introduced an additional constant, $\Theta$, which describes the
effective moment of inertia in the model. In our tethered flight
experiments, this corresponds to the inverse of the coupling gain, $g$, described in the manuscript. As can be
seen easily $\Theta$ and $C_\text{T}$ are not linearly independent, so for this
case only the sign of $\Theta$ will be important. The case where $\Theta$ is
set to a negative value corresponds to the natural condition in which the
generated torque moves the figure in the opposing direction, referred to as
_normal gain_.

## Equivalence with Poggio and Reichardt 1973 model for specific conditions ##

As mentioned in the main text, because this analytical model is based
explicitly on a pair of visual neurons, it can predict responses to arbitrary
visual stimuli in both open- and closed-loop. We found that for a particular
stimulus configuration (closed-loop figure fixation without input from a visual
background), it is formally equivalent to the classical model proposed by
Poggio and Reichardt based on torque measurements [Poggio, Reichardt (1973)](http://www.ncbi.nlm.nih.gov/pubmed/4718020), [Reichardt, Poggio (1975)](http://www.ncbi.nlm.nih.gov/pubmed/1138981)
(Eq. 2). For the specific notation refer to [Reichardt, Poggio (1975)](http://www.ncbi.nlm.nih.gov/pubmed/1138981). One can identify
$$
\begin{split}
R_{\text{F}}^{(-)}(\varphi_{\text{F}}) M^{(+)}(\dot{\varphi}_{\text{F}} / \lambda_{\text{F}})\; &\widehat{=} \; D^{\dagger}(\varphi(t), \dot{\varphi}(t)) \\
R_{\text{F}}^{(+)}(\varphi_{\text{F}}) M^{(-)}(\dot{\varphi}_{\text{F}} / \lambda_{\text{F}})\; &\widehat{=} \;\rho(\varphi(t), \dot{\varphi}(t))
\end{split}
$$ as the even ($f(x)=f(-x)$) and odd ($f(x)=-f(-x)$) symmetric parts of the torque response in $\dot{\varphi}$.

We now linearize the
Langevin equation around $\dot{\varphi}_{\text{F}} = 0$. Since our system is
described by stochastic processes we use [Kazakov (1961)](http://www.worldcat.org/oclc/561912994) to justify
the linearization. This assumes that the temporal average of the squared velocity does not vanish $\langle\dot{\varphi}^2\rangle > 0$ and
the figure is almost stationary $\dot{\varphi} \approx 0$, which is a valid assumption since we are only interested in the
quasi stationary behavior of the system.

For the even term, the linearization around $\dot{\varphi} \approx 0$
returns a constant and terms of the order of 2 or higher.
$$
\langle M^{(+)}(\dot{\varphi}) \rangle = M^{(+)}_0 + O(\dot{\varphi}^2)
$$
which depends on the initial slope of the motion response and on the
asymmetry, $C_{\text{ND}}$, of the motion responses. For a symmetric motion
response, this constant would be equal to zero. The asymmetry is therefore
required for the system to have a restoring force, which is essential
for figure fixation. Additionally, it is important to understand that
$M^{(+)}_0 = 0$ if $\langle\dot{\varphi}^2\rangle = 0$, i.e. in the absence of
noise. It is a common misunderstanding that the linearization would throw
away the asymmetric response to a moving figure.

After linearizing the even term, we can write:
$$
\frac{C_{\text{T}}}{2} R_{\text{F}}^{(-)}(\varphi_{\text{F}}) M^{(+)}_0 = \frac{\partial}{\partial \varphi} U(\varphi)
$$
Which describes the restoring force of the system. Since we assume the
potential energy of the dynamic system to be path independent, we can define a scalar potential, $U(\varphi)$, as written above.
Due to our circular boundary conditions, $U(\varphi) = U(\varphi + n 2\pi)$ is the periodic
potential in which the figure is moving. Please note that this potential
depends on the shape of the receptive field and, more importantly,
on the properties of the presented stimulus. It results from the stimulus properties and
a spatially integrating (wide-field) motion responsive cell having an anisotropic receptive field and
an asymmetric motion response.

For the odd term, the linearization around $\dot{\varphi} \approx 0$
returns only a linear term (see [Kazakov (1961)](http://www.worldcat.org/oclc/561912994)) and terms of the order of 3 or higher.
$$
\langle M^{(-)}(\dot{\varphi}) \rangle = M^{(-)}_0 \dot{\varphi} + O(\dot{\varphi}^3)
$$
The odd function $M^{(-)}(\dot{\varphi})$ is by definition already point
symmetric and shows no asymmetry. For small velocities this linearization
is valid.

After linearizing the odd term, we get:
$$
\frac{C_{\text{T}}}{2} R_{\text{F}}^{(+)}(\varphi_{\text{F}}) M^{(-)}_0 \dot{\varphi} \approx r \dot{\varphi}
$$
in [Reichardt, Poggio (1975)](http://www.ncbi.nlm.nih.gov/pubmed/1138981), they assume $R^{(+)}(\varphi) = r$ to be constant.
Since the odd term describes the dampening force of the system for quasi
stationary processes, this assumption is completely reasonable from
a modelling perspective (note that $R_{\text{F}}^{(+)}(\varphi) > 0$ is true for all
$\varphi$).

With all this, we can write down the linearized torque response to a figure as the Langevin equation:

$$\begin{split}
\Theta \ddot{\varphi}_{\text{F}} &= T^{*}(\varphi_{\text{F}}, \dot{\varphi}_{\text{F}}) \\
 &= \frac{C_{\text{T}}}{2} \Big[
    R_{\text{F}}^{(-)}(\varphi_{\text{F}}) M^{(+)}(\dot{\varphi}_{\text{F}} / \lambda_{\text{F}})
  + R_{\text{F}}^{(+)}(\varphi_{\text{F}}) M^{(-)}(\dot{\varphi}_{\text{F}} / \lambda_{\text{F}})
  \Big] + N(t) \\
  &\approx \frac{\partial}{\partial \varphi_{\text{F}}} U(\varphi_{\text{F}}) + r \dot{\varphi}_{\text{F}} + N(t)
\end{split}
$$

using standard techniques, this Langevin equation can be transformed into
a Fokker-Planck type equation, which describes the distribution arising 
from the time evolution of this stochastic differential equation. 
This representation allowed us to calculate probability distributions for
the figure fixation process.

$$
\frac{\partial}{\partial t} P = -\omega \frac{\partial}{\partial \varphi} P - \frac{\partial}{\partial \omega} \Big(\frac{\partial}{\partial \varphi} U + r \omega \Big) P + \zeta \frac{\partial^2}{\partial \omega^2} P
$$

For this case, closed-loop figure fixation, this equation is
equivalent to the theoretical framework presented in [Poggio, Reichardt
(1973)](http://www.ncbi.nlm.nih.gov/pubmed/4718020), [Reichardt, Poggio (1975)](http://www.ncbi.nlm.nih.gov/pubmed/1138981). Again, note that we did not derive this equation from
observed torque responses of tethered flies, but from the measured
response properties of the wide-field motion integrator cells --- the HS-cells --- in _Drosophila_. Furthermore, our spatially explicit
representation of visual stimuli and our explicit treatment of open- or
closed-loop dynamics allow use of our framework in other scenarios, as
performed below.

To summarize the above, the subtracted output of two opposing
_HS-like_ model cells predicts figure fixation in
_Drosophila_.

## Theoretical predictions #

In the following experiments, the fly is in closed-loop control of a figure
displayed in front of a moving background. It was previously shown
that flies fixate a figure in front of a moving background at a
position shifted from zero in the direction opposite to the background movement direction, and it was suggested
that this results from a summation of a motion independent response
and a motion response [(Bahl et al.,
2013)](http://dx.doi.org/10.1038/nn.3386).

However, we now show that our model, with only a single pathway, also
predicts figure displacement. The quasistationary dynamic equation for
our system looks like:

$$
\begin{split}
\Theta \ddot{\varphi}_{\text{F}} &= T(\varphi_{\text{F}}, \dot{\varphi}_{\text{F}}, \omega_{\text{G}}) \\
 &= \frac{C_{\text{T}}}{2} \Big[
    R_{\text{F}}^{(-)}(\varphi_{\text{F}}) M^{(+)}(\dot{\varphi}_{\text{F}} / \lambda_{\text{F}})
  + R_{\text{F}}^{(+)}(\varphi_{\text{F}}) M^{(-)}(\dot{\varphi}_{\text{F}} / \lambda_{\text{F}})
  \Big] \\
 &+ \frac{C_{\text{T}}}{2} \Big[
    R_{\text{G}}^{(-)}(\varphi_{\text{F}}) M^{(+)}(\omega_{\text{G}} / \lambda_{\text{G}})
  + R_{\text{G}}^{(+)}(\varphi_{\text{F}}) M^{(-)}(\omega_{\text{G}} / \lambda_{\text{G}})
  \Big]
\end{split}
$$

We would now like to perform a stability analysis of figure fixation. For
this, we calculate the points in space for which the torque
response vanishes and then check each of these fixed points to determine if it is a stable or an
unstable equilibrium. 

Since we are only interested in the fixed points of the system, the forces
in this equation should vanish, and we can again linearize around a stable
figure position and rewrite the dynamic equation to:

$$
\Theta \ddot{\varphi}_{\text{F}} = \frac{C_{\text{T}}}{2} \Big[
    \frac{\partial}{\partial \varphi_{\text{F}}} V(\varphi_{\text{F}}, \omega_{\text{G}}) + r \dot{\varphi}_{\text{F}}
\Big]
$$

with

$$\frac{\partial}{\partial \varphi_{\text{F}}} V(\varphi_{\text{F}}, \omega_{\text{G}}) =
\frac{\partial}{\partial \varphi_{\text{F}}} U(\varphi_{\text{F}})
  + R_{\text{G}}^{(-)}(\varphi_{\text{F}}) M^{(+)}(\omega_{\text{G}} / \lambda_{\text{G}})
  + R_{\text{G}}^{(+)}(\varphi_{\text{F}}) M^{(-)}(\omega_{\text{G}} / \lambda_{\text{G}})
$$

At the fixed points of the system, the force,
$\frac{\partial}{\partial \varphi_{\text{F}}} V(\varphi_{\text{F}}, \omega_{\text{G}})$,
vanishes. By solving for the positions where torque was zero, we calculated the
fixed points as we altered background velocity<!--@ .OPTIONALTEXT ", with the results shown
in Fig. 3A (top) of the main text" @-->.

<!--@ .INCLUDECODESNIPPET python phenomenological_model_tutorial_code.py figure3Atop @-->
<!--@ .OPTIONALFIGURE figure3Atop "Fixation behavior predicted by the analytical model for Drosophila, showing model responses to the closed-loop coupled figure. Dynamics of figure fixation dependent on background velocity are color-coded (regions with negative figure velocity are colored in purple, regions of positive velocities are shown in green). Positions at which responses are zero (fixed points) are indicated by a line, with the solid line denoting stable equilibrium — predicting figure fixation — and dashed an unstable equilibrium. Velocity is given relative to critical velocity $\omega_{\text{c}}$." @-->

We furthermore evaluated whether
these were stable or unstable. In Fig. 3A (top), the solid black line
is the stable position of a quasistationary figure dependent on the
background velocity. At this stable point, the restoring force in the
system pushes the figure back to the stable branch (indicated by the
color coding). The dashed black line is the unstable position of the
figure and perturbation from this equilibrium would grow and push
the figure away.

With increasing background velocity, the stable and unstable points
approach each other. At a critical background velocity,
$\omega_{\text{c}}$, the fixed points meet and vanish (Fig. 3A, top)
and therefore the figure can no longer be fixated. This is a classical
saddle-node bifurcation. At higher background velocities, the total
torque is dominated by the response component due to the open-loop
background, much like an optomotor response, and exceeds the fraction
generated by the response to the figure.

To calculate the probability distribution, $P(\varphi)$, of figure positions dependent on background velocity, we
solve the Fokker-Planck equation derived above, but for the potential, $V(\varphi_{\text{F}}, \omega_{\text{G}})$.
The analytical solution for this type of Fokker-Planck equation for
arbitrary potentials can be found in [Reichardt, Poggio 1975](http://www.ncbi.nlm.nih.gov/pubmed/1138981).
With this, we can calculate the probability distribution for the figure in closed-loop, and
compare it to the stationary solution of the bifurcation<!--@ .OPTIONALTEXT " (see Supplemental Fig. S2G and Fig. 3A bottom)"@-->.
Note how the mode of the distribution shifts in position and how the distribution becomes flatter when getting closer to the critical velocity.

<!--@ .INCLUDECODESNIPPET python phenomenological_model_tutorial_code.py figureS2G-3Abottom @-->
<!--@ .OPTIONALFIGURE figureS2G-3Abottom "Left: Probability of figure position under open-loop background conditions. Colors correspond to open-loop background velocity ranging from 0 (blue) to twice the critical background velocity (red). Right: Probability of figure position as a function of background velocity together with the black line indicating stable equilibrium." @-->


## Summarizing the phenomenological model #

We derived an analytical model for an asymmetric motion detection
pathway from fundamental response properties of HS-cells. We showed
that this model is capable of performing figure fixation in a manner
equivalent to a classical model [Poggio, Reichardt
(1973)](http://www.ncbi.nlm.nih.gov/pubmed/4718020), [Reichardt, Poggio (1975)](http://www.ncbi.nlm.nih.gov/pubmed/1138981). Furthermore, utilizing a novel aspect of our model, its
explicit dependence on arbitrary visual input, we also showed that it
predicts the displacement of the fixation point as a function of
background motion. This behavior was previously argued to result from
a motion independent position pathway [(Bahl et al., 2013)](http://dx.doi.org/10.1038/nn.3386). Thus, from a theoretical
perspective, a motion pathway alone leads to responses that have been
argued to require a position pathway.
