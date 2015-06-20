% A Short Note on Hidden Markov Models
% Dominic Steinitz
% 5th May 2015

\newcommand{\condprob} [3] {#1 \left( #2 \,\vert\, #3 \right)}

---
bibliography: Bayesian.bib
---

Markov Process and Chains
=========================

If you look at the [wikipedia
article](http://en.wikipedia.org/wiki/Hidden_Markov_model) on Hidden
Markov Models then you might be forgiven for concluding that these
deal only with discrete time and finite state spaces.

Recall that a **transition kernel** is a mapping $\mu : X \times {\cal{Y}}
\rightarrow \overline{\mathbb{R}}_{+}$ where $(X, {\cal{X}})$ and $(Y,
{\cal{Y}})$ are two measurable spaces such that $\mu(s, \cdot)$ is a
probability measure on ${\cal{Y}}$ for all $s \in X$ and such that
$\mu(\cdot, A)$ is a measurable function on $X$ for all $A \in
{\cal{Y}}$.

For example, we could have $X = Y = \{a, b\}$ and ${\cal{X}} =
{\cal{Y}} = \{\emptyset, \{a\}, \{b\}, \{a,b\}\}$ and $\mu(a,\{a\}) =
0.4, \mu(a,\{b\}) = 0.6, \mu(b,\{a\}) = 0.6, \mu(b,\{b\}) =
0.4$. Hopefully this should remind you of the transition matrix of a
Markov chain.

Recall further that a family of such transitions $\{\mu_t\}_{t \in T}$
where $T$ is some index set satisfying

$$
\begin{eqnarray}
\mu_{t+s}(x, A) = \int_{Y} \mu_s(x, {\mathrm{d}y})\mu_t(y, A)
\end{eqnarray}
$$

gives rise to a Markov process (under some mild conditions --- see
@rogerswilliams1 and @kallenberg2002foundations for much more detail),
that is, a process in which what happens next only depends on where the
process is now and not how it got there.

Let us carry on with our example and take $T = \mathbb{N}$. With a slight
abuse of notation and since $Y$ is finite we can re-write the integral
as a sum

$$
\begin{eqnarray}
\mu_{n+m}(x, z) = \sum_{y \in Y} \mu_m(x, y)\mu_n(y, z)
\end{eqnarray}
$$

which we recognise as a restatement of how Markov transition matrices
combine.

Dynamical Systems
=================

A dynamical system can be formulated as a dynamical system with a
particularly simple transition kernel given by

$$
\begin{equation}
\mu_t(x_s, A) = \delta(f_t(x_s), A) \triangleq
\begin{cases}
1           & \text{if } f_t(x_s) \in    A \\
0           & \text{if } f_t(x_s) \notin A
\end{cases}
\end{equation}
$$

where $f_t$ is the deterministic state update function (the flow) and
$\delta$ is the Dirac delta function.

Parameters
==========

Now let's add parameters to the dynamical system. Particle filters
work badly if we assume that the parameters are fixed; at some point
we end up with just one value for all the particles. Thus we pretend
that each parameter is not quite fixed but undergoes a small Brownian
motion or diffusion. This is called regularization in some of the
literature.

Using Greek letters for the parameters (and Roman letters for state),
we can write this more precisely as a transition kernel.

$$
\begin{equation}
\mu_t(\theta_s, {\mathrm{d}\phi}) =
\condprob{{\cal{N}}}{{\mathrm{d}\phi}}{\theta_s, \sigma^2(t-s)}
\end{equation}
$$

where we use e.g. ${\mathrm{d}\phi}$ to indicate probability densities.

Now suppose we believe that our parameters actually change over time
then as well as regularisation (which we use to make the particle
filter behave nicely) we could write

$$
\begin{equation}
\mu_t(\theta_s, {\mathrm{d}\phi}) =
\condprob{{\cal{N}}}{{\mathrm{d}\phi}}{\theta_s, \tau^2(t-s)}
\end{equation}
$$

This is *identical* to the regularisation equation but now we are
using it to describe our view on how parameters evolve rather than as
a particle filter trick.

Of course Brownian motion or diffusion may not be a good model for our
parameters; with Brownian motion, the parameters could drift off to
$\pm\infty$. We might believe that our parameters tend to stay close
to some given value and use the Ornstein-Uhlenbeck kernel.

$$
\mu_t(\theta_s, {\mathrm{d}\phi}) =
\condprob{{\cal{N}}}{{\mathrm{d}\phi}}{\alpha + (\theta_s - \alpha)e^{-\beta t},\frac{\sigma^2}{2\beta}\big(1 - e^{-2\beta t}\big)}
$$

where $\beta$ expresses how strongly we expect the parameter to
respond to perturbations, $\alpha$ is the mean to which the process
wants to revert (aka the asymptotic mean) and $\sigma^2$ expresses how
noisy the process is.

It is sometimes easier to view these transition kernels in terms of
stochastic differential equations. Brownian motion can be expressed as

$$
\mathrm{d}X_t = \sigma\mathrm{d}W_t
$$

and Ornstein-Uhlenbeck can be expressed as

$$
\mathrm{d}X_t = -\beta(X_t - \alpha)\mathrm{d}t + \sigma\mathrm{d}W_t
$$

where $W_t$ is the Wiener process.

Other
=====
From wikipedia

$$
\mathrm{d}x_t = \theta(\mu - x_t)\mathrm{d}t + \sigma\mathrm{d}W_t
$$

From Jeff's email

The OU update equations themselves are very straightforward, and are
given in terms of simple parameters $\lambda$ and $\sigma$.  The OU solves
the Ito SDE:

$$
\mathrm{d}X_t = -\lambda X_t\mathrm{d}t + \sqrt{2\lambda}\sigma\mathrm{d}W_t
\quad
X_0 \sim {\cal{N}}(0,\sigma^2)
$$

From Susanne

$$
\mathrm{d}X_t = -\beta(X_t - \alpha)\mathrm{d}t + \sigma\mathrm{d}W_t
$$

In integral form

$$
X_t = \alpha + (x_0 - \alpha)e^{-\beta t} + \sigma\int_0^t e^{-\beta(t - s)}\mathrm{d}W_s
$$

Assume without loss of generality that $s = 0$ then

$$
\mathbb{E}[X_t \,\vert\, X_0 = x_0] =
\mathbb{E}\Bigg[\alpha + (x_0 - \alpha)e^{-\beta t} + \sigma\int_0^t e^{-\beta(t - s)}\mathrm{d}W_s \,\vert\, X_0 = x_0\Bigg] =
\alpha + (x_0 - \alpha)e^{-\beta t}
$$

and

$$
\mathbb{V}[X_t \,\vert\, X_0 = x_0] =
\mathbb{E}\Bigg[\Bigg( \sigma\int_0^t e^{-\beta(t - s)}\mathrm{d}W_s  \Bigg)^2\Bigg]
$$

Using the [[http://en.wikipedia.org/wiki/It%C5%8D_isometry][Ito Isometry]]

$$
\mathbb{V}[X_t \,\vert\, X_0 = x_0] =
\sigma^2 \mathbb{E}\Bigg[ \int_0^t e^{-2\beta(t-s)}\mathrm{d}s  \Bigg] =
\frac{\sigma^2}{2\beta}\big(1 - e^{-2\beta t}\big)
$$

Summarising

$$
X_t \,\vert\, X_0 = x_0 \sim {\cal{N}}\bigg(\alpha + (x_0 - \alpha)e^{-\beta t},\frac{\sigma^2}{2\beta}\big(1 - e^{-2\beta t}\big)\bigg)
$$

Perhaps we can generalise slightly? Without loss of generality assume $t \leq u$

$$
\begin{aligned}
\mathbb{C}[X_u, X_t \,\vert\, X_0 = x_0] &=
\mathbb{E}\Bigg[\Bigg( \sigma\int_0^u e^{-\beta(u - s)}\mathrm{d}W_s  \Bigg)
                \Bigg( \sigma\int_0^t e^{-\beta(t - s)}\mathrm{d}W_s  \Bigg)\Bigg] \\
&=
\sigma^2e^{-\beta(u + t)}
\mathbb{E}\Bigg[\Bigg(\int_0^u e^{\beta s}\mathrm{d}W_s\Bigg)
          \Bigg(\int_0^t e^{\beta s}\mathrm{d}W_s\Bigg)\Bigg] \\
&=
\sigma^2e^{-\beta(u + t)}
\mathbb{E}\Bigg[\Bigg(\int_0^t e^{\beta s}\mathrm{d}W_s + \int_t^u e^{\beta s}\mathrm{d}W_s\Bigg)
          \Bigg(\int_0^t e^{\beta s}\mathrm{d}W_s\Bigg)\Bigg]
\end{aligned}
$$

And now we can use Ito and independence

$$
\begin{aligned}
\mathbb{C}[X_u, X_t \,\vert\, X_0 = x_0] &=
\sigma^2e^{-\beta(u + t)}\mathbb{E}\Bigg[ \int_0^t e^{2\beta s}\mathrm{d}s  \Bigg] \\
&=
\frac{\sigma^2e^{-\beta(u + t)}}{2\beta}\big(e^{2\beta t} - 1\big)
\end{aligned}
$$

Bibliography
============