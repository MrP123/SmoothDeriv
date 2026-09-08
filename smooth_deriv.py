# Matthias Panny - 2026-09-08

import numpy as np
from math import factorial, comb


class SmoothDeriv:
    """
    Smooth numerical differentiation via convolution kernels.

    Parameters
    ----------
    Ts : float
        Sampling period of the signal.
    T : float
        Kernel window length in seconds. Larger values give more smoothing. Must be designed based on the signal that is to be differentiated/filtered.
    n : int, optional
        Kernel order (default 4). Controls smoothness (kernel is 2n+1 times differentiable).
    r_max : int, optional
        Maximum derivative order to support (default 4).

    Examples
    --------
    >>> sd = SmoothDeriv(Ts=1e-3, T=5e-2)
    >>> y_smooth = sd(y, deriv_order=0)
    >>> y_dot = sd(y, deriv_order=1)
    >>> y_ddot = sd(y, deriv_order=2)
    """

    def __init__(self, Ts: float, T: float, n: int = 4, r_max: int = 4):
        self.Ts = Ts
        self.T = T
        self.n = n
        self.r_max = r_max

        self._phi_polys = self._build_polynomials()
        self._kernels = self._build_kernels()

    def _build_polynomials(self) -> list[np.ndarray]:
        """Build prototype polynomials for the convolution kernels and their derivatives.

        Returns
        -------
        list[np.ndarray]
            List of polynomials (i.e. array of coefficients) for the convolution kernels and their derivatives.
            Index 0 corresponds to the smoothing kernel, index 1 to the first derivative kernel, etc.
        """

        n, T = self.n, self.T
        max_order = 2 * n + 1
        Cn = factorial(max_order) / (factorial(n) ** 2)

        S_p = np.zeros(max_order + 1)  # prototype polynomial for the smoothing kernel
        for k in range(n + 1):
            idx = max_order - (n + k + 1)
            S_p[idx] = Cn * comb(n, k) * ((-1) ** k / (n + k + 1))

        phi_polys = [None] * (self.r_max + 1)
        phi_polys[0] = np.polyder(S_p)
        for r in range(self.r_max):
            phi_polys[r + 1] = np.polyder(phi_polys[r]) / T

        return phi_polys

    def _build_kernels(self) -> dict[int, np.ndarray]:
        """Build discrete convolution kernels from the prototype polynomials.

        Returns
        -------
        dict[int, np.ndarray]
            A dictionary mapping derivative orders to their corresponding discrete convolution kernels.
        """

        t_kernel = np.arange(0, self.T + self.Ts, self.Ts)
        tau_kernel = t_kernel / self.T

        norm_factor = np.sum(np.polyval(self._phi_polys[0], tau_kernel))

        kernels: dict[int, np.ndarray] = {}
        for r in range(self.r_max + 1):
            kernels[r] = np.polyval(self._phi_polys[r], tau_kernel) / norm_factor
        return kernels

    def kernel(self, deriv_order: int = 1) -> np.ndarray:
        """Return the discrete convolution kernel for the given derivative order."""
        if deriv_order not in self._kernels:
            raise ValueError(f"deriv_order={deriv_order} exceeds r_max={self.r_max}")
        return self._kernels[deriv_order]

    def __call__(self, y: np.ndarray, deriv_order: int = 1) -> np.ndarray:
        """
        Compute a smooth derivative of the signal.

        Parameters
        ----------
        y : np.ndarray
            Input signal.
        deriv_order : int, optional
            0 = smoothing, 1 = first derivative, etc. (default 1).

        Returns
        -------
        np.ndarray
            Result, same length as y.
        """
        return np.convolve(y, self.kernel(deriv_order), mode="same")

    def __repr__(self) -> str:
        return f"SmoothDeriv(Ts={self.Ts}, T={self.T}, n={self.n}, r_max={self.r_max})"
