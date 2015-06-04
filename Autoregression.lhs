% Predicting the Flow of the Thames
% Dominic Steinitz
% 29th May 2015

\newcommand{\condprob} [3] {#1 \left( #2 \,\vert\, #3 \right)}

---
bibliography: Bayesian.bib
---

```{.dia width='400'}
dia = image (DImage (ImageRef "diagrams/AutoregressionVary.png") 600 600 (translationX 0.0))
```

Thames Flux
===========

It is roughly 150 miles from the source of the Thames to [Kingston
Bridge](http://www.nationaltrail.co.uk/sites/default/files/thames_path_distances_source_to_teddington_jan_11.pdf). If
we assume that it flows at about 2 miles per second then the water
at Thames Head will have reached Kingston very roughly at
$\frac{150}{24\times 2} \approxeq 3$ days.

The Environmental Agency measure the flux at Kingston Bridge on a
twice daily basis. Can we predict this? In the first instance without
any other data and using our observation that Thames flushes itself
every 3 days, let us try

$$
X_t = \theta_1 X_{t-1} + \theta_2 X_{t-2} + \theta_3 X_{t-3} + \epsilon_t
$$

where $X_t$ is the flux on day $t$ and $\{\epsilon_t\}_{t \in
\mathbb{N}}$ are independent normal errors with mean $0$ and variance
some given value $\sigma^2$.

Kalman
======

As it stands, our model is not Markov so we cannot directly apply
techniques such as Kalman filtering or particle filtering to estimate
the parameters. However we can re-write the model as

$$
\begin{bmatrix}
\theta_1^{(t)} \\
\theta_2^{(t)} \\
\theta_3^{(t)}
\end{bmatrix} =
\begin{bmatrix}
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1
\end{bmatrix}
\begin{bmatrix}
\theta_1^{(t-1)} \\
\theta_2^{(t-1)} \\
\theta_3^{(t-1)}
\end{bmatrix} +
\begin{bmatrix}
\eta_{t} \\
\eta_{t} \\
\eta_{t}
\end{bmatrix}
$$

$$
y_t = \begin{bmatrix}
      x_{t-1} & x_{t-2} & x_{t-3}
      \end{bmatrix}
\begin{bmatrix}
\theta_1^{(t)} \\
\theta_2^{(t)} \\
\theta_3^{(t)}
\end{bmatrix} +
\epsilon_{t}
$$

Note that the observation map now varies over time so we have modify
our [Kalman filter
implementation](https://idontgetoutmuch.wordpress.com/2014/08/06/fun-with-kalman-filters-part-ii/)
to accept a different matrix at each step.

> {-# OPTIONS_GHC -Wall                     #-}
> {-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
> {-# OPTIONS_GHC -fno-warn-type-defaults   #-}
> {-# OPTIONS_GHC -fno-warn-unused-do-bind  #-}
> {-# OPTIONS_GHC -fno-warn-missing-methods #-}
> {-# OPTIONS_GHC -fno-warn-orphans         #-}

> {-# LANGUAGE DataKinds                    #-}
> {-# LANGUAGE ScopedTypeVariables          #-}
> {-# LANGUAGE RankNTypes                   #-}
> {-# LANGUAGE TypeOperators                #-}
> {-# LANGUAGE TypeFamilies                 #-}

> module Autoregression (
>     predictions
>   ) where

> import GHC.TypeLits
> import Numeric.LinearAlgebra.Static
> import Data.Maybe ( fromJust )

> import qualified Data.Vector as V

> inv :: (KnownNat n, (1 <=? n) ~ 'True) => Sq n -> Sq n
> inv m = fromJust $ linSolve m eye

> outer ::  forall m n . (KnownNat m, KnownNat n,
>                         (1 <=? n) ~ 'True, (1 <=? m) ~ 'True) =>
>           R n -> Sq n -> [L m n] -> Sq m -> Sq n -> Sq n -> [R m] -> [(R n, Sq n)]
> outer muPrior sigmaPrior bigHs bigSigmaY bigA bigSigmaX ys = result
>   where
>     result = scanl update (muPrior, sigmaPrior) (zip ys bigHs)
>
>     update :: (R n, Sq n) -> (R m, L m n) -> (R n, Sq n)
>     update (xHatFlat, bigSigmaHatFlat) (y, bigH) =
>       (xHatFlatNew, bigSigmaHatFlatNew)
>       where
>         v :: R m
>         v = y - bigH #> xHatFlat
>         bigS :: Sq m
>         bigS = bigH <> bigSigmaHatFlat <> (tr bigH) + bigSigmaY
>         bigK :: L n m
>         bigK = bigSigmaHatFlat <> (tr bigH) <> (inv bigS)
>         xHat :: R n
>         xHat = xHatFlat + bigK #> v
>         bigSigmaHat :: Sq n
>         bigSigmaHat = bigSigmaHatFlat - bigK <> bigS <> (tr bigK)
>         xHatFlatNew :: R n
>         xHatFlatNew = bigA #> xHat
>         bigSigmaHatFlatNew :: Sq n
>         bigSigmaHatFlatNew = bigA <> bigSigmaHat <> (tr bigA) + bigSigmaX

We can now set up the parameters to run the filter.

> stateVariance :: Double
> stateVariance = 1e-2

> bigSigmaX :: Sq 3
> bigSigmaX = fromList [ stateVariance, 0.0,           0.0
>                      , 0.0,           stateVariance, 0.0
>                      , 0.0,           0.0,           stateVariance
>                      ]

> bigA :: Sq 3
> bigA = eye

> muPrior :: R 3
> muPrior = fromList [0.0, 0.0, 0.0]

> sigmaPrior :: Sq 3
> sigmaPrior = fromList [ 1e1, 0.0, 0.0
>                       , 0.0, 1e1, 0.0
>                       , 0.0, 0.0, 1e1
>                       ]

> bigHsBuilder :: V.Vector Double -> [L 1 3]
> bigHsBuilder flows =
>   V.toList $
>   V.zipWith3 (\x0 x1 x2 -> fromList [x0, x1, x2])
>   (V.tail flows) (V.tail $ V.tail flows) (V.tail $ V.tail $ V.tail flows)

> obsVariance :: Double
> obsVariance = 1.0e-6

> bigSigmaY :: Sq 1
> bigSigmaY = fromList [ obsVariance ]

> predict :: R 3 -> Double -> Double -> Double -> Double
> predict theta f1 f2 f3 = h1 * f1 + h2 * f2 + h3 * f3
>   where
>     (h1, t1) = headTail theta
>     (h2, t2) = headTail t1
>     (h3, _)  = headTail t2

> thetas :: V.Vector Double -> [(R 3, Sq 3)]
> thetas flows = outer muPrior sigmaPrior (bigHsBuilder flows)
>                bigSigmaY bigA bigSigmaX (map (fromList . return) (V.toList flows))

> predictions :: V.Vector Double -> V.Vector Double
> predictions flows =
>   V.zipWith4 predict
>   (V.fromList $ map fst (thetas flows))
>   flows (V.tail flows) (V.tail $ V.tail flows)

How Good is Our Model?
======================

If we assume that parameters are essentially fixed by taking the state
variance to be e.g. $10^{-6}$ then the fit is not good.

```{.dia height='600'}
dia = image (DImage (ImageRef "diagrams/AutoregressionFix.png") 600 600 (translationX 0.0))
```

However, if we assume the parameters to undergo Brownian motion by
taking the state variance to be e.g. $10^{-2}$ then we get a much
better fit. Of course, Brownian motion is probably not a good way of
modelling the parameters; we hardly expect that these could wander off
to infinity.

```{.dia height='600'}
dia = image (DImage (ImageRef "diagrams/AutoregressionVary.png") 600 600 (translationX 0.0))
```
