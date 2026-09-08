# SmoothDeriv

Smooth numerical differentiation using polynomial convolution kernels.

`SmoothDeriv` implements a family of smooth differentiators based on a compactly supported polynomial kernel. Instead of differentiating the sampled signal directly, the signal is convolved with a kernel that performs smoothing and differentiation simultaneously.

This is particularly useful for experimental data, where ordinary finite differences can strongly amplify measurement noise.

## Mathematical Background

### Issue with numerical differentiation & the basic idea
A common engineering problem is computing the derivative of a measured signal

$$
y(t),
$$

for which the most straightforward approach is to use a finite difference approximation scheme.
With first order forward differences, this results in

$$
\dot{y}(t) \approx \frac{y(t+\Delta t)-y(t)}{\Delta t}.
$$

which has the disadvantage of amplifying high-frequency measurement noise, making it unsuitable for many real applications.

Instead this project makes use of a property of the convolution operation. For any two 
functions (or discrete signals) $f$ and $g$, it follows that

$$
(f * g)' = f' * g = f * g'.
$$

This means that instead of differentiating the signal directly, we can differentiate a kernel and convolve it with the signal.
This also holds for higher-order derivatives:

$$
f^{(r)} * g = f * g^{(r)}.
$$


As the kernel can be freely chosen, it can be constructed to have finite support  and a prescribed degree of smoothnes at its boundaries. The resulting derivative estimate is therefore smooth and less sensitive to noise.

This means, that we now have to find a suitable method to construct a family of kernels $\phi$ that can be used to compute derivatives of any order $r$.
The proposed nomenclature for this is to denote the kernel for the $r$-th derivative, for the kernel order $n$ as $\phi_{r,n}$.
Therefore, it holds that

$$
\dot{y}(t) \approx y(t) * \phi_{1,n}(t).
$$


### Construction of the kernel
The kernel is constructed from a smooth polynomial prototype function $p_n(\tau)$ defined on the normalized interval $0 \leq \tau \leq 1$, where $\tau = \frac{t}{T}$ for the kernel with width $T$. 

The polynomial is constructed such that it fulfills the following properties:
- $p_n(0) = 0, \quad p_n(1) = 1$
- $p_n'(0) = 0, \quad p_n'(1) = 0$
- up to the $n$-th derivative:
- $p_n^{(n)}(0) = 0, \quad p_n^{(n)}(1) = 0$
This gives $2(n+1) = 2n+2$ boundary conditions, in a similar fashion to Hermite interpolation. This gives the required coefficients for a polynomial of order $2n+1$.

This is best formulated in terms of the derivative of the polynomial

$$
\phi_{0,n}(\tau) = p_n'(\tau) = C_n\tau^n(1-\tau)^n
$$

with $C_n$ being chosen in a way that the kernel is normalized, i.e.:

$$
\int_0^1 \phi_{0,n}(\tau)\,d\tau=1.
$$

From this normalization it follows that

$$
\int_0^1\tau^n(1-\tau)^n\,d\tau = B(n+1,n+1) = \frac{(n!)^2}{(2n+1)!},
$$

where $B$ denotes the Beta function.

Consequently,

$$
C_n = \frac{1}{B(n+1,n+1)} = \frac{(2n+1)!}{(n!)^2}
$$

Thus $\phi_{0,n}$ is a normalized smoothing kernel.

#### Example for $n=4$
For the default value $n=4$,

$$
C_4 = \frac{9!}{(4!)^2} = 630.
$$

Therefore,

$$
\phi_{0,4}(\tau) = 630 \, \tau^4(1-\tau)^4.
$$

Expanded,

$$
\phi_{0,4}(\tau) = 630 \left( \tau^8 -4\tau^7 +6\tau^6 -4\tau^5 +\tau^4 \right). 
$$

The corresponding polynomial $p_4$ is

$$
p_4(\tau) = 70\tau^9 -315\tau^8 +540\tau^7 -420\tau^6 +126\tau^5.
$$

#### Example for $n = 1$
For $n = 1$, we find

$$
C_1 = \frac{3!}{(1!)^2} = 6.
$$

Therefore,

$$
\phi_{0,1}(\tau) = 6 \, \tau^1(1-\tau)^1 = 6\tau(1-\tau),
$$

with the corresponding polynomial $p_1$ being

$$
p_1(\tau) = \frac{6}{2}\tau^2 - \frac{6}{3}\tau^3 = 3\tau^2 - 2\tau^3.
$$

This is also known as the cubic Hermite spline or "smoothstep" function in computer graphics.

### Derivatives of the kernel
Differentiation the smoothing kernel $\phi_{0,n}(\tau)$ with respect to $t$ requires the chain rule, since $\tau = t/T$:

$$
\frac{d}{dt} \phi_{0,n}(\tau) =  \frac{d}{dt} \phi_{0,n}\left(\frac{t}{T}\right) = \frac{d}{d\tau} \phi_{0,n}(\tau) \, \frac{d\tau}{dt} = \frac{1}{T} \phi_{0,n}'(\tau).
$$

As $\phi_{0,n}(\tau)$ is the derivative of the prototype polynomial $p_n(\tau)$, we can also write

$$
\frac{d}{dt}
p_n\left(\frac{t}{T}\right)
=
\frac{1}{T}p_n'(\tau),
$$

and more generally,

$$
\frac{d^r}{dt^r} p_n\left(\frac{t}{T}\right) = \frac{1}{T^r} \, p_n^{(r)}(\tau).
$$

The resulting polynomials represent

$$
\phi_{0,n} = p_n',
$$

$$
\phi_{1,n} = \frac{1}{T}p_n'',
$$

$$
\phi_{2,n} = \frac{1}{T^2}p_n''',
$$

and so on.

The index therefore corresponds to the desired derivative order:

$$
\phi_{r,n} = \frac{1}{T^r} \, p_n^{(r+1)}.
$$

As $p_n$ is a polynomial this can be done analytically without any numerical differentiation.


### Discrete convolution kernels
For a sampled signal with sampling period $T_s$ the kernel is sampled at

$$
t_k = k \, T_s,
$$

as such the corresponding normalized coordinate is

$$
\tau_k=\frac{t_k}{T}.
$$

For derivative order $r$, the discrete convolution kernel $h$ is then given by

$$
h_{r,n}[k] = \phi_{r,n}(\tau_k) \frac{1}{\Phi_\Sigma}.
$$

The smoothing kernel is normalized using the discrete sum

$$
\Phi_\Sigma = \sum_k \phi_{0,n}(\tau_k),
$$

which is the discrete equivalent of the continuous integral.
This gives the final discrete convolution kernels.

### Applying the derivative kernel
Once the kernel has been constructed, differentiation is performed through discrete convolution.

For derivative order $r$,

$$
y^{(r)}[j] \approx \sum_k y[j-k] \, h_{r,n}[k].
$$

and for the first order derivative in particular,

$$
\dot{y}[j] \approx \sum_k y[j-k] \, h_{1,n}[k].
$$


## Parameter selection

### Kernel width $T$
The parameter $T$ determines the width of the convolution kernel in seconds.
A larger $T$ means that more samples contribute to each derivative estimate.
This means, that larger $T$ results in more smoothing and less sensitivity to noise, while smaller $T$ results in less smoothing and better preservation of rapid signal changes.

### Kernel order $n$
The parameter $n$ determines the order of the polynomial kernel and therefore its smoothness. Increasing $n$ changes the shape and smoothness of the kernel.
The kernel becomes increasingly concentrated toward the center of the interval while retaining zero values at the boundaries.


## Provided code
The provided code implements the kernel construction and convolution in Python as module `smooth_deriv.py` with test cases shown in `kernel_convolve_deriv.ipynb`.
There is also a MATLAB implementation in `kernel_convolve_deriv.m`.