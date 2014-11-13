
import numpy as np
import scipy.linalg
import scipy.signal

np.seterr(divide='ignore')


# helper function
def arange_endpoint( start, stop, inc, dtype=np.float64, eps=1e-16 ):
    """create an arrange from start to stop, including endpoints"""
    return np.arange( start, stop+eps, inc, dtype=dtype )


# filter recovery
def prony(h, nb, na):
    """Prony's method for time-domain IIR filter design.

    Description:

      Finds a filter with numerator order na, denominator order nb,
      and having the impulse response in array h.  The IIR filter
      coefficients are returned in length nb+1 and na+1 row vectors b
      and a, ordered in descending powers of Z.

    Inputs:

      h --- impulse response to fit
      nb, na -- number of filter coefficients

    Outputs:

      b,a -- Numerator and denominator of the iir filter.
    """
    h = np.asarray(h)
    K = len(h) - 1
    M = nb
    N = na
    if K < max(M,N):
        raise ValueError("Model order too large for data. Zero-pad data to fix?")
    c = h[0] if h[0] != 0 else 1 # avoid divide by zero
    row = np.zeros((K+1,))
    row[0] = (h/c)[0] # avoid scipy warning
    H = scipy.linalg.toeplitz(h/c,row)
    if K > N:
        H = H[:, :N+1]
    # Partition H matrix
    H1 = H[:(M+1),:]
    h1 = H[(M+1):,:1]
    H2 = H[(M+1):,1:]
    x,resids,rank,s = np.linalg.lstsq(-H2,h1)
    a = np.vstack(([1],x))[:,0]
    b = np.dot(H1, c*a[:, np.newaxis])[:,0]
    return b, a


class FilterMaker(object):

    def __init__(self, hz, impulse_response_cutoff=20, test_filter=False, normalize=True, incremental_increase_till_valid=False):
        """FilterMaker( framerate_hz )

        > iir_lowpass1

        > iir_highpass1

        > james_lmc

        """
        self.hz = float(hz)
        self.dt = 1.0 / self.hz
        self._impulse_response_cutoff = float(impulse_response_cutoff) # default at np.exp(-20) ~= 2e-09
        self._test_filter = bool(test_filter)
        self._filter_order_lowpass = (0, 1)
        self._filter_order_highpass = (1, 1)
        self._filter_order_lmc = (17, 13)
        #self._filter_order_bandpass = (30, 20)
        self._normalize = normalize
        self._iitv = incremental_increase_till_valid

    def iir_lowpass1_direct(self, tau=0.008):
        """first order low-pass IIR filter"""
        max_t = tau * self._impulse_response_cutoff
        t = arange_endpoint(0, max_t, self.dt)
        V = 1.0/tau * np.exp(-t/tau)
        # normalize to unity gain
        if self._normalize: V = V/abs(np.sum(V))
        return V

    def iir_lowpass1(self, tau=0.008):
        V = self.iir_lowpass1_direct(tau)
        b, a = prony(V, *self._filter_order_lowpass)
        if self._test_filter:
            if not self._iitv:
                assert self._test_filter_func(V, (b, a))
            else:
                while not self._test_filter_func(V, (b, a)):
                    self._filter_order_lowpass = tuple(map(lambda x: x+1, self._filter_order_lowpass))
                    b, a = prony(V, *self._filter_order_lowpass)
                    print "LOWPASS FILTER CORRECTION (b,a):", (b,a)
        return b, a

    def iir_highpass1_direct(self, tau=0.05):
        """first order high-pass IIR filter"""
        max_t = tau * self._impulse_response_cutoff
        t = arange_endpoint(0, max_t, self.dt)
        V = -1.0/tau * np.exp(-t/tau)
        # normalize to unity gain
        if self._normalize: V = V/abs(np.sum(V))
        V[0] = V[0]+1.0 # --> dirac delta at t=0
        return V

    def iir_highpass1(self, tau=0.05):
        V = self.iir_highpass1_direct(tau)
        b, a = prony(V, *self._filter_order_highpass)
        if self._test_filter:
            if not self._iitv:
                assert self._test_filter_func(V, (b, a))
            else:
                while not self._test_filter_func(V, (b, a)):
                    self._filter_order_highpass = tuple(map(lambda x: x+1, self._filter_order_highpass))
                    b, a = prony(V, *self._filter_order_highpass)
                    print "HIHGPASS FILTER CORRECTION (b,a):", (b,a)
        return b, a

    def james_lmc_direct(self, a1=-1.06,tau1=0.012,sigma1=0.197,
                               a2=0.167,tau2=0.021,sigma2=0.345):
        # see Lindemann, et al. 2005
        # NOTE: With default parameters, this is not a perfect high pass
        # and has a lowpass gain of approximately -0.003.
        max_t = max( tau1*np.exp(np.sqrt(2)*sigma1*np.sqrt(self._impulse_response_cutoff)),
                     tau2*np.exp(np.sqrt(2)*sigma1*np.sqrt(self._impulse_response_cutoff)) )
        t = arange_endpoint(0, max_t, self.dt)
        V= ( a1 * np.exp(-(np.log(t/tau1))**2/(2*sigma1**2)) 
           + a2 * np.exp(-(np.log(t/tau2))**2/(2*sigma2**2)) )
        # normalize to unity gain
        if self._normalize: V = V/abs(np.sum(V))
        return V

    def james_lmc(self, a1=-1.06,tau1=0.012,sigma1=0.197,
                        a2=0.167,tau2=0.021,sigma2=0.345):
        V = self.james_lmc_direct(a1, tau1, sigma1, a2, tau2, sigma2)
        b, a = prony(V, *self._filter_order_lmc)
        if self._test_filter:
            if not self._iitv:
                assert self._test_filter_func(V, (b, a))
            else:
                while not self._test_filter_func(V, (b, a)):
                    self._filter_order_lmc = tuple(map(lambda x: x+1, self._filter_order_lmc))
                    b, a = prony(V, *self._filter_order_lmc)
                    print "LMC FILTER CORRECTION (b,a):", (b,a)
        return b, a

    def _test_filter_func(self, direct_impulse_response, ba, rtol=1e-05, atol=1e-06):
        # ensure order of filter is high enough
        b, a = ba
        dirac = np.zeros_like(direct_impulse_response)
        dirac[0] = 1
        indirect_impulse_response = scipy.signal.lfilter(b, a, dirac)
        return np.allclose(direct_impulse_response, indirect_impulse_response, rtol, atol)



if __name__ == "__main__":

    import matplotlib.pyplot as plt

    sample_rate = 1000. # Hz
    filters = FilterMaker(sample_rate, test_filter=True)

    # Impulse responses
    lp1_ir = filters.iir_lowpass1_direct()
    hp1_ir = filters.iir_highpass1_direct()
    lmc_ir = filters.james_lmc_direct()
    IRs = [lp1_ir, hp1_ir, lmc_ir]

    # Numerator/Denominator representation
    lp1_ba = filters.iir_lowpass1()
    hp1_ba = filters.iir_highpass1()
    lmc_ba = filters.james_lmc()
    BAs = [lp1_ba, hp1_ba, lmc_ba]


    # plotting axes
    ax0 = plt.subplot2grid((3, 1), (0, 0))
    ax1 = plt.subplot2grid((3, 1), (1, 0))
    ax2 = plt.subplot2grid((3, 1), (2, 0))
    AXs = [ax0, ax1, ax2]

    for ir, (b, a), ax in zip(IRs, BAs, AXs):
        # time
        t = np.arange(len(ir)) * filters.dt # time
        # delta function arrays for impulse reponse recovery
        d = np.zeros_like(t)
        d[0] = 1
        # recovered impulse responses
        rir = scipy.signal.lfilter(b, a, d)

        ax.plot(t, ir, 'r-', linewidth=2)
        ax.plot(t, rir, 'b--', linewidth=1)

        ax.set_xlabel('time (s)')
        ax.set_ylabel('response')

    plt.show()

     # plotting axes
    ax0 = plt.subplot2grid((3, 1), (0, 0))
    ax1 = plt.subplot2grid((3, 1), (1, 0))
    ax2 = plt.subplot2grid((3, 1), (2, 0))
    AXs = [ax0, ax1, ax2]

    for ir, (b, a), ax in zip(IRs, BAs, AXs):
        # time
        t = np.arange(len(ir)) * filters.dt # time
        p = 20*np.log10(np.abs(np.fft.fft(ir)))
        freqs = np.fft.fftfreq(ir.size, 1/512.)
        freqs = freqs[freqs>0]
        p = p[freqs>0]
        idx = np.argsort(freqs)
        ax.semilogx(freqs[idx], p[idx])
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('Magnitude (dB)')

    plt.show()
