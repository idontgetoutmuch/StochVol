% Stochastic Volatility
% Dominic Steinitz
% 28th December 2014

\newcommand{\condprob} [3] {#1 \left( #2 \,\vert\, #3 \right)}

---
bibliography: Bayesian.bib
---

Introduction
============

In their tutorial on particle filtering, @doucet2011tutorial give an
example of stochastic volatiltiy.

$$
X_1 \sim {\mathcal{N}}\bigg(0, \frac{\sigma^2}{1 - \alpha^2}\bigg)
$$

Stochastic volatility models treat the volatility (i.e., variance) of
a return on an asset, such as an option to buy a security, as
following a latent stochastic process in discrete time
\citep{KimShephardChib:1998}.  The data consist of mean corrected
(i.e., centered) returns $y_t$ on an underlying asset at $T$ equally
spaced time points.  Kim et al.\ formulate a typical stochastic
volatility model using the following regression-like equations, with a
latent parameter $h_t$ for the log volatility, along with parameters
$\mu$ for the mean log volatility, and $\phi$ for the persistence of
the volatility term.  The variable $W_t$ represents the
white-noise shock (i.e., multiplicative error) on the asset return at
time $t$, whereas $V_t$ represents the shock on volatility at
time $t$.

$$
\begin{aligned}
X_{n+1} &= \mu + \alpha (X_n - \mu) + V_n \sigma \\
Y_n     &= \beta \exp(X_t / 2) W_n \\
h_1     &\sim {\mathcal{N}}\left( \mu, \frac{\sigma}{\sqrt{1 - \phi^2}} \right)
\end{aligned}
$$

$$
W_t \sim {\mathcal{N}}(0,1); \ \ \ \ \  V_t \sim {\mathcal{N}}(0,1)
$$

Rearranging the first line, $W_t = y_t \exp(-h_t / 2)$,
allowing the sampling distribution for $y_t$ to be written as
\[ 
y_t \sim {\mathcal{N}}(0,\exp(h_t/2)).
\]
The recurrence equation for $h_{t+1}$ may be combined with the
scaling and sampling of $V_t$ to yield the sampling distribution
\[
h_t \sim {\mathcal{N}}(\mu + \phi(h_t - \mu), \sigma).
\]

TBD
===

From the state equation, we have

$$
\begin{align}
X_{n+1} &= \mu + \alpha(X_n - \mu)  + V_n \sigma \\
\alpha^2(X_n - \mu) &= \alpha(X_{n+1} - \mu) -  \alpha V_n \sigma \\
\end{align}
$$

We also have

$$
\begin{align}
X_n - \mu &=  \alpha(X_{n-1} - \mu)  + V_{n-1} \sigma \\
\end{align}
$$

Adding the two expressions together gives

$$
\begin{align}
(1+\alpha^2)(X_n - \mu) &= \alpha((X_{n-1} - \mu) + (X_{n+1} - \mu)) + \sigma (V_{n-1} - \alpha V_n) \\
(X_n - \mu) &= \frac{\alpha}{1+\alpha^2}((X_{n-1} - \mu) + (X_{n+1} - \mu)) + \frac{\sigma}{1+\alpha^2} (V_{n-1} - \alpha V_n) \\
X_n &= \mu + \frac{\alpha}{1+\alpha^2}((X_{n-1} - \mu) + (X_{n+1} - \mu)) + \frac{\sigma}{1+\alpha^2} (V_{n-1} - \alpha V_n)
\end{align}
$$

Since $\{V_n\}$ are standard normal, then $X_n$ conditional on
$X_{n-1}$ and $X_{n+1}$ is normally distributed, and

$$
\begin{align}
\mathbb{E}(X_n\mid X_{n-1}, X_{n+1}) &= \mu + \frac{\alpha}{1+\alpha^2}((X_{n-1} - \mu) + (X_{n+1} - \mu)) \\
\mathbb{V}(X_n\mid X_{n-1}, X_{n+1}) &= \frac{\sigma^2}{1+\alpha^2}
\end{align}
$$

Scribbles
=========

$$
f(h_t | h_{t^-}, \theta, y) = \frac{f(h, \theta, y)}{f(h_{t^-}, \theta, y)}
$$

$$
f(h_t | h_{t^-}, y) = \frac{f(h, y)}{f(h_{t^-}, y)}
$$

$$
f(h_t | h_{t^-})f(y_t | h_t) = \frac{f(h)}{f(h_{t^-})} \frac{f(y_t, h_t)}{f(h_t)}
$$

MCMC
====

To construct the chain we specify a non-informative prior for the parameters

$$
p(\alpha, \beta, \eta) \propto \eta^{-2}
$$

Standard [Bayesian analysis for
regression](http://en.wikipedia.org/wiki/Bayesian_linear_regression)
tells us that the conditional distribution is given by

$$
\condprob {p}{\alpha, \beta, \eta}{y^{(t)}, x^{(t)}} = {\cal{N}}((\alpha, \beta); \mu_n, \sigma^2\Lambda_n^{-1})\,{\cal{IG}}(a_n, b_n)
$$

with

$$
\Lambda_n = X^\top X  + \Lambda_0
$$

$$
\Lambda_n = \boldsymbol{x}_n^\top \boldsymbol{x}_n  + \Lambda_{n-1}
$$

In the case of our model we can specialise these as

$$
\Lambda_n = \begin{bmatrix} 1 & 1 & \ldots & 1 \\
                            x_1 & x_2 & \ldots & x_n
            \end{bmatrix}
            \begin{bmatrix} 1 & x_1 \\
                            1 & x_2 \\
                            \ldots & \ldots \\
                            1 & x_n
            \end{bmatrix}
+ \Lambda_0
$$

<span style="color: red">colour</span>

$$
\color{blue}{
B_t^{-1} = B_{t-1}^{-1} + \begin{bmatrix} 1 \\ x_{t-1} \end{bmatrix}
                          \begin{bmatrix} 1 &  x_{t-1} \end{bmatrix}
}
$$

$$
\color{blue}
{y_i = \boldsymbol{x}_i^{T}\boldsymbol{\beta} + \epsilon_i}
$$

Bibliography
============