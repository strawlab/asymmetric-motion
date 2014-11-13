
import numpy as np
import scipy.sparse
import os
import sys

from .wrapping import wrap_to

import progressbar



class SimpleBlurMatrix(object):
    """SimpleBlurMatrix

    This class is equivalent to the FastBlurMatrix, but easier to understand.
 
    """
    def __init__(self, input_phi, input_theta, output_phi, output_theta, sigma):
        # calculate the weights
        weights_az = np.exp(-wrap_to(input_phi[:,np.newaxis] - output_phi, -np.pi, np.pi)**2 / (2*sigma**2))
        weights_el = np.exp(-wrap_to(input_theta[:,np.newaxis] - output_theta, -np.pi, np.pi)**2 / (2*sigma**2))
        # store in class
        self._A = weights_az.copy()
        self._B = weights_el.T.copy()

    def blur(self, INPUT):
        """gaussian blur on the INPUT array"""
        return np.dot(np.dot(self._B, INPUT), self._A)


class FastBlurMatrix(object):
    """FastBlurMatrix

    This class is equivalent to the SimpleBlurMatrix, but faster.

    """
    def __init__(self, input_phi, input_theta, output_phi, output_theta, sigma, cachefile='.cached_blurmatrix.npz'):
        # generate shapes.
        self.ishape = (input_phi.size*input_theta.size,)
        self.oshape = (output_theta.size, output_phi.size)
        self.shape = np.array((output_phi.size * output_theta.size,
                                input_phi.size * input_theta.size))

        self._cache = str(cachefile)
        # attempt to load cached blur matrix:
        try:
            cache = np.load(self._cache)
            assert np.all(cache['shape'] == self.shape)
        except IOError: # file not found don't load cache
            pass
        except AssertionError: # grid_changed. remove cache
            print >> sys.stderr, 'emd_blurmatrix: stale cache, deleting...'
            os.unlink(self._cache)
        else:
            data = cache['data']
            indices = cache['indices']
            indptr = cache['indptr']
            self._csr_bm = scipy.sparse.csr_matrix((data, indices, indptr), shape=self.shape)
            return

        # generate the blurmatrix
        self._csr_bm = self._generate_blur_matrix(input_phi, input_theta,
                                            output_phi, output_theta, sigma)

    def _generate_blur_matrix(self, ip, it, op, ot, sigma, eps=1e-7):
        """generate_blur_matrix

        Parameters
        ----------
        ip : ndarray
          Input phi angles 1d
        it : ndarray
          Input theta angles 1d
        op : ndarray
          Output phi angles 1d
        ot : ndarray
          Output theta angles 1d

        Returns
        -------
        csr : scipy.sparse.csr_matrix
          A sparse matrix to generate the blur on an
          input array by multiplication only.

        """
        ophi, otheta = np.meshgrid(op, ot)
        iphi, itheta = np.meshgrid(ip, it)
        ophi = ophi.reshape((-1,1))
        iphi = iphi.reshape((-1,1))
        otheta = otheta.reshape((-1,1))
        itheta = itheta.reshape((-1,1))

        indptr = [0]
        indices = []
        data = []
        cum_len = 0

        print "First run: generatig blurmatrix"
        pbar = progressbar.ProgressBar()
        for out_phi, out_theta in pbar(zip(ophi, otheta)):
            dist2 = wrap_to(out_phi - iphi, -np.pi, np.pi)**2 + (out_theta - itheta)**2
            weight = np.exp(-(dist2) / (2*sigma**2))
            weight /= np.sum(weight)  # unity gain
            cond = weight > eps
            nz = np.nonzero(cond)[0]
            cum_len += len(nz)
            indptr.append(cum_len)
            indices.append(nz)
            data.append(weight[cond])

        data = np.hstack(data)
        indices = np.hstack(indices)
        indptr = np.array(indptr)
        # save cached blur matrix:
        np.savez(self._cache, data=data, indices=indices, indptr=indptr, shape=self.shape)
        return scipy.sparse.csr_matrix((data, indices, indptr), shape=self.shape)


    def blur(self, INPUT):
        """blur input array onto new grid"""
        result = self._csr_bm * INPUT.reshape((-1,1))
        return result.reshape(self.oshape)




if __name__ == '__main__':

    import matplotlib.pyplot as plt

    ax0 = plt.subplot2grid((3,1),(0,0))
    ax1 = plt.subplot2grid((3,1),(1,0))
    ax2 = plt.subplot2grid((3,1),(2,0))

    X0 = np.linspace(-np.pi, np.pi, 512, endpoint=False)
    Y0 = np.linspace(-np.pi/3., np.pi/3., 256)
    extent = (X0[0], X0[-1], Y0[0], Y0[-1])
    stim = np.clip(np.cos(8*(np.outer(Y0, X0)))*1e10, 0, 1)
    ax0.imshow(stim, extent=extent, cmap=plt.get_cmap('gray'), interpolation='nearest')

    X1 = np.linspace(-np.pi, np.pi, 67, endpoint=False)
    Y1 = np.linspace(-np.pi/3., np.pi/3., 31)

    BM = SimpleBlurMatrix(X0, Y0, X1, Y1, np.radians(2))
    blurredstim = BM.blur(stim)
    ax1.imshow(blurredstim, extent=extent, cmap=plt.get_cmap('gray'), interpolation='nearest')

    BM2 = FastBlurMatrix(X0, Y0, X1, Y1, np.radians(2))
    blurredstim2 = BM2.blur(stim)
    ax2.imshow(blurredstim2, extent=extent, cmap=plt.get_cmap('gray'), interpolation='nearest')

    plt.show()

    import time

    s = time.time()
    for _ in range(100):
        blurredstim = BM.blur(stim)
    print "2 * np.dot takes on average %f sec" % ((time.time()-s)/100.)
    s = time.time()
    for _ in range(100):
        blurredstim2 = BM2.blur(stim)
    print "csr takes on average %f sec" % ((time.time()-s)/100.)




