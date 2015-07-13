% Conditional Expectation under Change of Measure
% Dominic Steinitz
% 12th July 2015

\newcommand{\condprob} [3] {#1 \left( #2 \,\vert\, #3 \right)}

---
bibliography: Bayesian.bib
---

Theorem
-------

Let $\mathbb{P}$ and $\mathbb{Q}$ be measures on $(\Omega, {\mathcal{F}})$
with $\mathbb{Q} \ll \mathbb{P}$, ${\mathcal{G}} \subset {\mathcal{F}}$ a sub
$\sigma$-algebra and $X$ an integrable random variable
($\mathbb{P}\lvert{X}\rvert < \infty$) then

$$
\mathbb{P}(X\vert {\mathcal{G}}) =
\frac
{\mathbb{Q}\bigg(X\frac{\mathrm{d}\mathbb{P}}{\mathrm{d}\mathbb{Q}}\bigg\vert {\mathcal{G}}\bigg)}
{\mathbb{Q}\bigg(\frac{\mathrm{d}\mathbb{P}}{\mathrm{d}\mathbb{Q}}\bigg\vert {\mathcal{G}}\bigg)}
$$

Proof
-----

$$
\begin{aligned}
\mathbb{Q}\bigg(\mathbb{I}_A \mathbb{Q}\bigg(X \frac{\mathrm{d}\mathbb{P}}{\mathrm{d}\mathbb{P}}\bigg\vert {\mathcal{G}}\bigg)\bigg)
&=
\mathbb{Q}\bigg(\mathbb{I}_A X \frac{\mathrm{d}\mathbb{P}}{\mathrm{d}\mathbb{P}}\bigg) \\
&=
\mathbb{P}\big(\mathbb{I}_A X \big) \\
&=
\mathbb{P}\big(\mathbb{I}_A \mathbb{P}(X \vert {\mathcal{G}})\big) \\
&=
\mathbb{Q}\bigg(\mathbb{I}_A \frac{\mathrm{d}\mathbb{P}}{\mathrm{d}\mathbb{Q}}\mathbb{P}(X \vert {\mathcal{G}})\bigg) \\
&=
\mathbb{Q}\bigg(\mathbb{I}_A \mathbb{Q}\bigg(\frac{\mathrm{d}\mathbb{P}}{\mathrm{d}\mathbb{Q}}\bigg\vert {\mathcal{G}}\bigg)\mathbb{P}(X \vert {\mathcal{G}})\bigg) \\
\end{aligned}
$$

Thus

$$
\mathbb{Q}\bigg(\mathbb{I}_A\mathbb{P}(X\vert {\mathcal{G}}){\mathbb{Q}\bigg(\frac{\mathrm{d}\mathbb{P}}{\mathrm{d}\mathbb{Q}}\bigg\vert {\mathcal{G}}\bigg)}\bigg) =
\mathbb{Q}\bigg(\mathbb{I}_A \mathbb{Q}\bigg(X\frac{\mathrm{d}\mathbb{P}}{\mathrm{d}\mathbb{Q}}\bigg\vert {\mathcal{G}}\bigg)\bigg)\quad \mathrm{for\,all}\, A \in {\mathcal{G}}
$$

Hence

$$
\mathbb{P}(X\vert {\mathcal{G}}){\mathbb{Q}\bigg(\frac{\mathrm{d}\mathbb{P}}{\mathrm{d}\mathbb{Q}}\bigg\vert {\mathcal{G}}\bigg)} =
\mathbb{Q}\bigg(X\frac{\mathrm{d}\mathbb{P}}{\mathrm{d}\mathbb{Q}}\bigg\vert {\mathcal{G}}\bigg)\quad {\mathbb{Q}-\mathrm{a.s.}}
$$

We note that

$$
A = \bigg\{\omega \,\bigg\vert\, \mathbb{Q}\bigg(\frac{\mathrm{d}\mathbb{P}}{\mathrm{d}\mathbb{Q}}\bigg)\bigg\}
$$

is ${\mathcal{G}}$-measurable (it is the result of a projection) and that

$$
0
=
\mathbb{Q}\bigg(\mathbb{I}_A\mathbb{Q}\bigg(
\frac{\mathrm{d}\mathbb{P}}{\mathrm{d}\mathbb{Q}}
\bigg\vert {\mathcal{G}}\bigg)\bigg)
=
\mathbb{Q}\bigg(\mathbb{I}_A
\frac{\mathrm{d}\mathbb{P}}{\mathrm{d}\mathbb{Q}}
\bigg)
=
\mathbb{P}(\mathbb{I}_A)
$$

Hence

$$
\mathbb{P}(X\vert {\mathcal{G}}) =
\frac
{\mathbb{Q}\bigg(X\frac{\mathrm{d}\mathbb{P}}{\mathrm{d}\mathbb{Q}}\bigg\vert {\mathcal{G}}\bigg)}
{\mathbb{Q}\bigg(\frac{\mathrm{d}\mathbb{P}}{\mathrm{d}\mathbb{Q}}\bigg\vert {\mathcal{G}}\bigg)}\quad {\mathbb{P}-\mathrm{a.s.}}
$$

as required.