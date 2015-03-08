% Stochastic Volatility
% Dominic Steinitz
% 28th December 2014

\newcommand{\condprob} [3] {#1 \left( #2 \,\vert\, #3 \right)}

---
bibliography: Bayesian.bib
---

Introduction
============

Simple models for e.g. financial option pricing assume that the
volatility of an index or a stock is constant, see
[here](https://idontgetoutmuch.wordpress.com/2013/02/10/parallelising-path-dependent-options-in-haskell-2/)
for example. However, simple observation of time series show that this
is not the case.

One approach which addresses this, GARCH (Generalised AutoRegressive
Conditional Heteroskedasticity), models the evolution of volatility
deterministically.

Stochastic volatility models treat the volatility of a return on an
asset, such as an option to buy a security, as a Hidden Markov Model
(HMM). Typically, the observable data consist of noisy mean-corrected
returns on an underlying asset at equally spaced time points.

There is evidence that Stochastic Volatility models
(@RePEc:bla:restud:v:65:y:1998:i:3:p:361-93) offer increased
flexibility over the GARCH family, e.g. see
@Geweke94bayesiancomparison, @citeulike:2576658 and
@Jacquier94bayesiananalysis. Despite this and judging by the numbers
of questions on the R Special Interest Group on Finance [mailing
list](https://stat.ethz.ch/pipermail/r-sig-finance/), the use of GARCH
in practice far outweighs that of Stochastic Volatility. Reasons cited
are the multiplicity of estimation methods for the latter and the lack
of packages (but see
[here](http://cran.r-project.org/web/packages/stochvol/vignettes/article.pdf)
for a recent improvement to the paucity of packages).

In their tutorial on particle filtering, @doucet2011tutorial give an
example of stochastic volatility. We save this approach for future
blog posts and follow [Lopes and
Polson](http://faculty.chicagobooth.edu/nicholas.polson/research/papers/lopes-polson-2010.pdf)
and the excellent [lecture
notes](http://hedibert.org/wp-content/uploads/2013/12/UPCcourse-handouts.pdf)
by [Hedibert Lopes](http://hedibert.org).

Here's the model.

$$
\begin{aligned}
H_0     &\sim {\mathcal{N}}\left( m_0, C_0\right) \\
H_t     &= \mu + \phi H_{t-1} + \tau \eta_t \\ 
Y_n     &= \beta \exp(H_t / 2) \epsilon_n \\
\end{aligned}
$$

We wish to estimate $\mu, \phi, \tau$ and $\boldsymbol{h}$. To do this
via a Gibbs sampler we need to sample from

$$
\condprob{p}{\mu, \phi, \tau}{\boldsymbol{h}, \boldsymbol{y}} \quad \text{and} \quad
\condprob{p}{\boldsymbol{h}}{\mu, \phi, \tau, \boldsymbol{y}}
$$

Marginal Distribution for Parameters
====================================

Let us take a prior that is standard for linear regression

$$
(\boldsymbol{\theta}, \tau^2) \sim {\mathcal{NIG}}(\boldsymbol{\theta}_0, V_0, \nu_0, s_0^2)
$$

where

$$
\boldsymbol{\theta} = \begin{bmatrix} \mu \\ \phi \end{bmatrix}
$$

and use standard results for linear regression to obtain the required
marginal distribution.

That the prior is Normal Inverse Gamma (${\cal{NIG}}$) means

$$
\begin{aligned}
\boldsymbol{\theta} \, | \, \tau^2 & \sim {\cal{N}}(\boldsymbol{\theta}_0, V_0) \\
\tau^2                             & \sim {\cal{IG}}(\nu_0 / 2, \nu_0 s_0^2 / 2)
\end{aligned}
$$

Standard [Bayesian analysis for
regression](http://en.wikipedia.org/wiki/Bayesian_linear_regression)
tells us that the (conditional) posterior distribution for

$$
y_i = \beta + \alpha x_i + \epsilon_i
$$

where the $\{\epsilon_i\}$ are IID normal with variance $\sigma^2$ is
given by

$$
\condprob {p}{\alpha, \beta, \eta}{\boldsymbol{y}, \boldsymbol{x}} =
{\cal{N}}((\alpha, \beta); \mu_n, \sigma^2\Lambda_n^{-1})\,{\cal{IG}}(a_n, b_n)
$$

with

$$
\Lambda_n = X_n^\top X_n  + \Lambda_0
$$

$$
\begin{matrix}
\mu_n = \Lambda_n^{-1}({X_n}^{\top}{X_n}\hat{\boldsymbol{\beta}}_n + \Lambda_0\mu_0) &
\textrm{where} &
\hat{\boldsymbol\beta}_n = (\boldsymbol{X}_n^{\rm T}\boldsymbol{X}_n)^{-1}\boldsymbol{X}_n^{\rm T}\boldsymbol{y}_n
\end{matrix}
$$

$$
\begin{matrix}
a_n = \frac{n}{2} + a_0 & \quad &
b_n = b_0 +
      \frac{1}{2}(\boldsymbol{y}^\top\boldsymbol{y} +
                  \boldsymbol{\mu}_0^\top\Lambda_0\boldsymbol{\mu}_0 -
                  \boldsymbol{\mu}_n^\top\Lambda_n\boldsymbol{\mu}_n)
\end{matrix}
$$

Recursive Form
--------------

We can re-write the above recursively. We do not need to for this blog
article but it will be required in any future blog article which uses
Sequential Monte Carlo techniques.

$$
\Lambda_n = \boldsymbol{x}_n^\top \boldsymbol{x}_n  + \Lambda_{n-1}
$$

Furthermore

$$
\Lambda_{n}\mu_{n} =
\boldsymbol{X}_{n}^{\rm T}\boldsymbol{y}_{n} + \Lambda_0\mu_0 =
\boldsymbol{X}_{n-1}^{\rm T}\boldsymbol{y}_{n-1} + \boldsymbol{x}_n^\top y_n + \Lambda_0\mu_0 =
\Lambda_{n-1}\mu_{n-1} + \boldsymbol{x}_n^\top y_n
$$

so we can write

$$
\boldsymbol{\mu}_n = \Lambda_n^{-1}(\Lambda_{n-1}\mu_{n-1} + \boldsymbol{x}_n^\top y_n)
$$

and

$$
\begin{matrix}
a_n = a_{n-1} + \frac{1}{2} & \quad &
b_n = b_{n-1} + \frac{1}{2}\big[(y_n - \boldsymbol{\mu}_n^\top \boldsymbol{x}_n)y_n + (\boldsymbol{\mu}_{n-1} - \boldsymbol{\mu}_{n})^\top \Lambda_{n-1}\boldsymbol{\mu}_{n-1}\big]
\end{matrix}
$$

Specialising
------------

In the case of our model we can specialise the non-recursive equations
as

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

$$
\begin{matrix}
\mu_n = (\Lambda_n)^{-1}({X_n}^{\top}{X_n}\hat{\boldsymbol{\beta}}_n + \Lambda_0\mu_0) &
\textrm{where} &
\hat{\boldsymbol\beta}_n = (\boldsymbol{X}_n^{\rm T}\boldsymbol{X}_n)^{-1}\boldsymbol{X}_n^{\rm T}\boldsymbol{x}_{2:n+1}
\end{matrix}
$$

$$
\begin{matrix}
a_n = \frac{n}{2} + a_0 & \quad &
b_n = b_0 +
      \frac{1}{2}(\boldsymbol{x}_{2:n+1}^\top\boldsymbol{x}_{2:n+1} +
                  \boldsymbol{\mu}_0^\top\Lambda_0\boldsymbol{\mu}_0 -
                  \boldsymbol{\mu}_n^\top\Lambda_n\boldsymbol{\mu}_n)
\end{matrix}
$$

Let's re-write the notation to fit our model.

$$
\Lambda_n = \begin{bmatrix} 1 & 1 & \ldots & 1 \\
                            h_1 & h_2 & \ldots & h_n
            \end{bmatrix}
            \begin{bmatrix} 1 & h_1 \\
                            1 & h_2 \\
                            \ldots & \ldots \\
                            1 & h_n
            \end{bmatrix}
+ \Lambda_0
$$

$$
\begin{matrix}
\mu_n = (\Lambda_n)^{-1}({H_n}^{\top}{H_n}\hat{\boldsymbol{\theta}}_n + \Lambda_0\mu_0) &
\textrm{where} &
\hat{\boldsymbol\theta}_n = (\boldsymbol{X}_n^{\rm T}\boldsymbol{X}_n)^{-1}\boldsymbol{X}_n^{\rm T}\boldsymbol{x}_{2:n+1}
\end{matrix}
$$

$$
\begin{matrix}
a_n = \frac{n}{2} + a_0 & \quad &
b_n = b_0 +
      \frac{1}{2}(\boldsymbol{x}_{2:n+1}^\top\boldsymbol{x}_{2:n+1} +
                  \boldsymbol{\mu}_0^\top\Lambda_0\boldsymbol{\mu}_0 -
                  \boldsymbol{\mu}_n^\top\Lambda_n\boldsymbol{\mu}_n)
\end{matrix}
$$

Marginal Distribution for State
===============================

From the state equation, we have

$$
\begin{align}
H_{t+1} &=  \mu + \phi H_{t} + \tau \eta_t \\ 
\phi^2 H_{t} &=  -\phi\mu + \phi H_{t+1} - \phi \tau \eta_t \\ 
\end{align}
$$

We also have

$$
\begin{align}
H_{t} &=  \mu + \phi H_{t-1} + \tau \eta_{t-1} \\ 
\end{align}
$$

Adding the two expressions together gives

$$
\begin{align}
(1 + \phi^2)H_{t} &= \phi (H_{t-1} + H_{t+1}) + \mu (1 - \phi) + \tau(\eta_{t-1} - \phi\eta_t) \\
H_{t} &= \frac{\phi}{1 + \phi^2} (H_{t-1} + H_{t+1}) + \mu \frac{1 - \phi}{1 + \phi^2} + \frac{\tau}{1 + \phi^2}(\eta_{t-1} - \phi\eta_t) \\
\end{align}
$$

Since $\{\eta_t\}$ are standard normal, then $H_t$ conditional on
$H_{t-1}$ and $H_{t+1}$ is normally distributed, and

$$
\begin{align}
\mathbb{E}(H_n\mid H_{n-1}, H_{n+1}) &= \frac{1 - \phi}{1+\phi^2}\mu +
                                        \frac{\phi}{1+\phi^2}(H_{n-1} + H_{n+1}) \\
\mathbb{V}(H_n\mid H_{n-1}, H_{n+1}) &= \frac{\tau^2}{1+\phi^2}
\end{align}
$$

Writing

$$
\boldsymbol{h}_{-t} = \begin{bmatrix}
                      h_0, &
                      h_1, &
                      \ldots, &
                      h_{t-1}, &
                      h_{t+1}, &
                      \ldots, &
                      h_{n-1}, &
                      h_n
                      \end{bmatrix}
$$

by Bayes' Theorem we have

$$
\condprob{p}{h_t}{\boldsymbol{h}_{-t}, \theta, \boldsymbol{y}} \propto
\condprob{p}{y_t}{h_t} \condprob{p}{h_t}{\boldsymbol{h}_{-t}, \theta, y_{-t}} =
f_{\cal{N}}(y_t;0,e^{h_t}) f_{\cal{N}}(h_t;\mu_t,\nu_t^2)
$$

where $f_{\cal{N}}(x;\mu,\sigma^2)$ is the probability density
function of a normal distribution.

Markov Chain Monte Carlo
========================

Now we can write down a single step for our Gibbs sampler. Steps 1 - 3
specify the Metropolis step within this.

1. For each $t$, sample $h_t^\flat$ from ${\cal{N}}(h_t, \gamma^2)$ where $\gamma^2$
is the tuning variance.

2. For each $t$, compute the acceptance probability

$$
p_t = \min{\Bigg(\frac{f_{\cal{N}}(h^\flat_t;\mu_t,\nu_t^2) f_{\cal{N}}(y_t;0,e^{h^\flat_t})}{f_{\cal{N}}(h_t;\mu_t,\nu_t^2) f_{\cal{N}}(y_t;0,e^{h_t})}, 1 \Bigg)}
$$

3. For each $t$, compute a new value of $h_t$

$$
h^\sharp_t =
\begin{cases}
h^\flat_t \text{with probability } p_t \\
h_t \text{with probability } 1 - p_t
\end{cases}
$$

4. Sample from $\condprob{p}{\boldsymbol{\theta},
\tau^2}{\boldsymbol{h}, \boldsymbol{y}} \sim
{\mathcal{NIG}}(\boldsymbol{\theta}_1, V_1, \nu_1, s_1^2)$

Haskell Preamble
================

> {-# OPTIONS_GHC -Wall                      #-}
> {-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
> {-# OPTIONS_GHC -fno-warn-type-defaults    #-}
> {-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
> {-# OPTIONS_GHC -fno-warn-missing-methods  #-}
> {-# OPTIONS_GHC -fno-warn-orphans          #-}

> {-# LANGUAGE RecursiveDo                   #-}
> {-# LANGUAGE ExplicitForAll                #-}
> {-# LANGUAGE TypeOperators                 #-}
> {-# LANGUAGE TypeFamilies                  #-}
> {-# LANGUAGE ScopedTypeVariables           #-}
> {-# LANGUAGE DataKinds                     #-}
> {-# LANGUAGE FlexibleContexts              #-}

> module StochVol (
>     bigM
>   , bigM0
>   , vols
>   , ys
>   , runMC
>   ) where

> import Numeric.LinearAlgebra.HMatrix hiding ( (===), (|||), Element )
> import qualified Numeric.LinearAlgebra.Static as S
> import GHC.TypeLits
> import Data.Proxy
> import Data.Maybe ( fromJust )

> import Data.Random
> import Data.Random.Source.PureMT
> import Control.Monad.Fix
> import Control.Monad.State.Lazy
> import Control.Monad.Writer
> import Control.Monad.Loops
> import Control.Applicative

> import qualified Data.Vector as V


Let's create some test data.

> mu', phi', tau2', tau' :: Double
> mu'   = -0.00645
> phi'  =  0.99
> tau2' =  0.15^2
> tau'  = sqrt tau2'

> n :: Int
> n = 5000

Arbitrarily let us start the process at

> h0 :: Double
> h0 = 0.0

We define the process as a stream (aka co-recursively) using the
Haskell *recursive do* construct. It is not necessary to do this but
streams are a natural way to think of stochastic processes.

> hs, vols, sds, ys :: V.Vector Double
> hs = V.fromList $ take n $ fst $ runState hsAux (pureMT 1)
>   where
>     hsAux :: (MonadFix m, MonadRandom m) => m [Double]
>     hsAux = mdo { x0 <- sample (Normal (mu' + phi' * h0) tau')
>                 ; xs <- mapM (\x -> sample (Normal (mu' + phi' * x) tau')) (x0:xs)
>                 ; return xs
>                 }

> vols = V.map exp hs

We can plot the volatility (which we cannot observe directly).

```{.dia height='500'}
import Data.Vector ( toList )

import StochVol
import StochVolChart

dia = diag 1.0 "Volatility" (zip (map fromIntegral [0..]) (toList vols))
```

And we can plot the log returns.

> sds = V.map sqrt vols

> ys = fst $ runState ysAux (pureMT 2)
>   where
>     ysAux = V.mapM (\sd -> sample (Normal 0.0 sd)) sds

```{.dia height='500'}
import Data.Vector ( toList )

import StochVol
import StochVolChart

dia = diag 1.0 "Log Return" (zip (map fromIntegral [0..]) (toList ys))
```

> m0, c0 :: Double
> m0 = 0.0
> c0 = 100.0

FIXME: Why start the sampler with the values of the simulated process?

> mu0, phi0, tau20 :: Double
> mu0   = -0.1 -- -0.00645
> phi0  =  0.95 -- 0.99
> tau20 =  0.5 -- 0.15^2

> bigV0, invBigV0 :: S.Sq 2
> bigV0 = S.diag $ S.fromList [100.0, 100.0]
> invBigV0 = sinv bigV0

> nu0, s02 :: Double
> nu0    = 10.0
> s02    = (nu0 - 2) / nu0 * tau20

Tuning parameter

> vh :: Double
> vh = 0.1

General MCMC setup
------------------

> bigM, bigM0 :: Int
> bigM0 = 2000
> bigM  = 3000

> iC0 :: Double
> iC0 = recip c0

> iC0m0 :: Double
> iC0m0  = iC0 * m0

> type StatsM a = RVarT (Writer [((Double, Double), Double)]) a

> (|||) :: (KnownNat ((+) r1 r2), KnownNat r2, KnownNat c, KnownNat r1) =>
>          S.L c r1 -> S.L c r2 -> S.L c ((+) r1 r2)
> (|||) = (S.¦)

> runMC :: [((Double, Double), Double)]
> runMC = take bigM $ drop bigM0 $
>         execWriter (evalStateT (sample multiStep) (pureMT 42))

> multiStep :: StatsM (Double, Double, Double, Double, V.Vector Double)
> multiStep = iterateM_ (singleStep vh ys) (mu0, phi0, tau20, h0, vols)

> singleStep :: Double -> V.Vector Double ->
>               (Double, Double, Double, Double, V.Vector Double) ->
>               StatsM (Double, Double, Double, Double, V.Vector Double)
> singleStep vh y (mu, phi, tau2, h0, h) = do
>   lift $ tell [((mu, phi),tau2)]
>   hNew <- randomWalkMetropolis y mu phi tau2 h0 h vh
>   let var = recip $ (iC0 + phi^2 / tau2)
>       mean = var * (iC0m0 + phi * ((hNew V.! 0) - mu) / tau2)
>   h0New <- rvarT (Normal mean (sqrt var))
>   let bigX :: S.L 5000 2 -- FIXME
>       bigX = (S.col $ S.vector $ replicate n 1.0)
>              |||
>              (S.col $ S.vector $ V.toList $ h0 `V.cons` V.init hNew)
>   newParms <- sampleParms (S.vector $ V.toList h) bigX (S.vector [mu0, phi0]) invBigV0 nu0 s02
>   return ( (S.extract (fst newParms))!0
>          , (S.extract (fst newParms))!1
>          , snd newParms
>          , h0New
>          , hNew
>          )

> randomWalkMetropolis :: V.Vector Double ->
>                         Double ->
>                         Double ->
>                         Double ->
>                         Double ->
>                         V.Vector Double ->
>                         Double ->
>                         RVarT m (V.Vector Double)
> randomWalkMetropolis ys mu phi tau2 h0 hs vh = do
>   let eta2s = V.replicate (n-1) (tau2 / (1 + phi^2)) `V.snoc` tau2
>       etas  = V.map sqrt eta2s
>       coef1 = (1 - phi) / (1 + phi^2) * mu
>       coef2 = phi / (1 + phi^2)
>       mu_n  = mu + phi * (hs V.! (n-1))
>       mu_1  = coef1 + coef2 * ((hs V.! 1) + h0)
>       innerMus = V.zipWith (\hp1 hm1 -> coef1 + coef2 * (hp1 + hm1)) (V.tail (V.tail hs)) hs
>       mus = mu_1 `V.cons` innerMus `V.snoc` mu_n
>   hs' <- V.mapM (\mu -> rvarT (Normal mu vh)) hs
>   let num1s = V.zipWith3 (\mu eta h -> logPdf (Normal mu eta) h) mus etas hs'
>       num2s = V.zipWith (\y h -> logPdf (Normal 0.0 (exp (0.5 * h))) y) ys hs'
>       nums  = V.zipWith (+) num1s num2s
>       den1s = V.zipWith3 (\mu eta h -> logPdf (Normal mu eta) h) mus etas hs
>       den2s = V.zipWith (\y h -> logPdf (Normal 0.0 (exp (0.5 * h))) y) ys hs
>       dens = V.zipWith (+) den1s den2s
>   us <- V.replicate n <$> rvarT StdUniform
>   let ls   = V.zipWith (\n d -> min 0.0 (n - d)) nums dens
>   return $ V.zipWith4 (\u l h h' -> if log u < l then h' else h) us ls hs hs'

> sampleParms ::
>   forall n m .
>   (KnownNat n, (1 <=? n) ~ 'True) =>
>   S.R n -> S.L n 2 -> S.R 2 -> S.Sq 2 -> Double -> Double ->
>   RVarT m (S.R 2, Double)
> sampleParms y bigX theta_0 invBigV_0 a_0 s_02 = do
>   let n = natVal (Proxy :: Proxy n)
>       a_n = 0.5 * (a_0 + fromIntegral n)
>       invBigV_n = invBigV_0 + (tr bigX) S.<> bigX
>       bigV_n = sinv invBigV_n
>       theta_n = bigV_n S.#> ((tr bigX) S.#> y + (tr invBigV_0) S.#> theta_0)
>       r = y - bigX S.#> theta_n
>       s = r `S.dot` r
>       b_0 = 0.5 * a_0 * s_02
>       p21 = a_0 * s_02 + s
>       p22 = d `S.dot` (invBigV_0 S.#> d) where d = theta_n - theta_0
>       b_n = 0.5 * (p21 + p22)
>       b_m = b_0 +
>             0.5 * (S.extract (S.row y S.<> S.col y)!0!0) +
>             0.5 * (S.extract (S.row theta_0 S.<> invBigV_0 S.<> S.col theta_0)!0!0) -
>             0.5 * (S.extract (S.row theta_n S.<> invBigV_n S.<> S.col theta_n)!0!0)
>   error ("\n\n" ++
>          show b_0 ++ " " ++ show b_n ++ " " ++ show b_m)
>   g <- rvarT (Gamma a_n (recip b_n))
>   let s2 = recip g
>   let bigV_n' = m S.<> bigV_n
>         where
>           m = S.diag $ S.vector (replicate 2 s2)
>   m1 <- rvarT StdNormal
>   m2 <- rvarT StdNormal
>   let theta_n' :: S.R 2
>       theta_n' = theta_n + S.chol (S.sym bigV_n') S.#> (S.vector [m1, m2])
>   return (theta_n', s2)

> sinv :: (KnownNat n, (1 <=? n) ~ 'True) => S.Sq n -> S.Sq n
> sinv m = fromJust $ S.linSolve m S.eye

```{.dia height='600'}
dia = image (DImage (ImageRef "mus.png") 600 600 (translationX 0.0))
```

```{.dia height='400'}
dia = image (DImage (ImageRef "phis.png") 400 400 (translationX 0.0))
```

Bibliography
============
