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
X_0     &\sim {\mathcal{N}}\bigg(0, \frac{\sigma^2}{1 - \alpha^2}\bigg) \\
X_{n+1} &= \mu + \alpha (X_n - \mu) + V_n \sigma \\
H_0     &\sim {\mathcal{N}}\left( \mu, \frac{\sigma}{\sqrt{1 - \phi^2}} \right) \\
H_t     &= \mu + \phi H_{t-1} + \tau \eta_t \\ 
Y_n     &= \beta \exp(X_t / 2) W_n \\
\end{aligned}
$$

Temporarily

$$
m_0 = \mu \quad C_0 = \frac{\sigma}{\sqrt{1 - \phi^2}}
$$

Let

$$
\boldsymbol{\theta} = \begin{bmatrix} \mu \\ \phi \end{bmatrix}
$$

Let us take a prior that is standard for linear regression

$$
(\boldsymbol{\theta}, \tau^2) \sim {\mathcal{NIG}}(\boldsymbol{\theta}_0, V_0, \nu_0, s_0^2)
$$

In other words

$$
\begin{aligned}
\boldsymbol{\theta} \, | \, \tau^2 & \sim {\cal{N}}(\boldsymbol{\theta}_0, V_0) \\
\tau^2                             & \sim {\cal{IG}}(\nu_0 / 2, \nu_0 s_0^2 / 2)
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
H_{t+1} &=  \mu + \phi H_{t} + \tau \eta_t \\ 
\phi^2 H_{t} &=  -\phi\mu + \phi H_{t+1} - \phi \tau \eta_t \\ 
\end{align}
$$

We also have

$$
\begin{align}
X_n - \mu &=  \alpha(X_{n-1} - \mu)  + V_{n-1} \sigma \\
H_{t} &=  \mu + \phi H_{t-1} + \tau \eta_{t-1} \\ 
\end{align}
$$

Adding the two expressions together gives

$$
\begin{align}
(1+\alpha^2)(X_n - \mu) &= \alpha((X_{n-1} - \mu) + (X_{n+1} - \mu)) + \sigma (V_{n-1} - \alpha V_n) \\
(1 + \phi^2)H_{t} &= \phi (H_{t-1} + H_{t+1}) + \mu (1 - \phi) + \tau(\eta_{t-1} - \phi\eta_t) \\
(X_n - \mu) &= \frac{\alpha}{1+\alpha^2}((X_{n-1} - \mu) + (X_{n+1} - \mu)) + \frac{\sigma}{1+\alpha^2} (V_{n-1} - \alpha V_n) \\
X_n &= \mu + \frac{\alpha}{1+\alpha^2}((X_{n-1} - \mu) + (X_{n+1} - \mu)) + \frac{\sigma}{1+\alpha^2} (V_{n-1} - \alpha V_n) \\
H_{t} &= \frac{\phi}{1 + \phi^2} (H_{t-1} + H_{t+1}) + \mu \frac{1 - \phi}{1 + \phi^2} + \frac{\tau}{1 + \phi^2}(\eta_{t-1} - \phi\eta_t) \\
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

Algorithm
---------

$$
\condprob{p}{x_t}{x_{t-1}, x_{t+1}, \boldsymbol{y}} \propto
\exp{\bigg\{-\frac{1}{2}\bigg[ x_n + \frac{y_n^2}{e^{x_n}} + (x_t - m_t)^2/b^2 \bigg]  \bigg\}}
$$

where

$$
b^2 = \frac{\sigma^2}{1 + \beta^2}
$$

$$
m_t = \frac{\alpha(1- \beta) + \beta(x_{t-1} + x_{t+1})}{1 + \beta^2}
$$

$$
\alpha = \mu - \alpha\mu
$$

$$
\beta = \alpha
$$

$$
b^2 = \frac{\sigma^2}{1 + \alpha^2}
$$

$$
m_t = \frac{(\mu - \alpha\mu)(1- \alpha) + \alpha(x_{t-1} + x_{t+1})}{1 + \alpha^2}
$$



MCMC
====

Standard Analysis
-----------------

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

We also have

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

Hyperparameters
===============

Choose $X_0 \sim {\cal{N}}(m_0, C_0)$

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


> mu, phi, tau2, tau :: Double
> mu   = -0.00645
> phi  =  0.99
> tau2 =  0.15^2
> tau  = sqrt(tau2)

> n :: Int
> n = 5000

> h0 :: Double
> h0 = 0.0

> hs, vols, sds, ys :: V.Vector Double
> hs = V.fromList $ take n $ fst $ runState hsAux (pureMT 1)
>   where
>     hsAux :: (MonadFix m, MonadRandom m) => m [Double]
>     hsAux = mdo { x0 <- sample (Normal (mu + phi * h0) tau)
>                 ; xs <- mapM (\x -> sample (Normal (mu + phi * x) tau)) (x0:xs)
>                 ; return xs
>                 }

> vols = V.map exp hs

```{.dia height='500'}
import Data.Vector ( toList )

import StochVol
import StochVolChart

dia = diag 1.0 "Volatility" (zip (map fromIntegral [0..]) (toList vols))
```

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

$\theta = (\mu, \phi)^\top$ then

$$
\begin{matrix}
\theta \,|\, \tau^2 & \sim & {\cal{N}}(\theta_0, \tau^2V_0) \\
\tau^2              & \sim & {\cal{IG}}(\nu_0/2, \nu_0 s_0^2/2)
\end{matrix}
$$

> mu0, phi0 :: Double
> mu0   = -0.00645
> phi0  =  0.99

> bigV0, invBigV0 :: S.Sq 2
> bigV0 = S.diag $ S.fromList [100.0, 100.0]
> invBigV0 = sinv bigV0

> nu0, s02 :: Double
> nu0    = 10.0
> s02    = (nu0 - 2) / nu0 * tau2

Tuning parameter

> vh :: Double
> vh = 0.1

General MCMC setup
------------------

> bigM, bigM0 :: Int
> bigM0 = 1000
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
> multiStep = iterateM_ (singleStep vh ys) (mu0, phi0, tau2, h0, vols)

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
> sampleParms y bigX b bigA v lam = do
>   let n = natVal (Proxy :: Proxy n)
>       p1 = 0.5 * (v + fromIntegral n)
>       var = sinv $ bigA + (tr bigX) S.<> bigX
>       mean = var S.#> ((tr bigX) S.#> y + (tr bigA) S.#> b)
>       r = y - bigX S.#> mean
>       s = r `S.dot` r
>       p21 = v * lam + s
>       p22 = d `S.dot` (bigA S.#> d) where d = mean - b
>       p2 = 0.5 * (p21 + p22)
>   g <- rvarT (Gamma p1 (recip p2))
>   let s2 = recip g
>   let var' = m S.<> var
>         where
>           m = S.diag $ S.vector (replicate 2 s2)
>   m1 <- rvarT StdNormal
>   m2 <- rvarT StdNormal
>   let mean' :: S.R 2
>       mean' = mean + S.chol (S.sym var') S.#> (S.vector [m1, m2])
>   return (mean', s2)

> sinv :: (KnownNat n, (1 <=? n) ~ 'True) => S.Sq n -> S.Sq n
> sinv m = fromJust $ S.linSolve m S.eye

```{.dia width='800'}
import Data.Histogram ( asList )
import StochVolMain

dia = barDiag (zip (map fst $ asList (hist mus)) (map snd $ asList (hist mus)))
  where
    result = runMC
    mus = map (fst. fst) result
```


Bibliography
============
