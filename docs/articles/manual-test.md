# Manual test

## Purpose

This validation reconstructs the fitted parameters of
[`gipslda()`](https://antonikingston.github.io/gipsDA/reference/gipslda.md),
[`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md),
and
[`gipsmultqda()`](https://antonikingston.github.io/gipsDA/reference/gipsmultqda.md)
independently. The reference calculation uses only base R, `stats`, and
the public `gips` API. It does not call any `gipsDA` helper or internal
function.

The executable version is
`tests/testthat/manual/test-manual-parameter-reconstruction.R`. The
formulas and implementation notes are kept beside the calculations in
that test so that differences between the reference and package
implementations remain reviewable.

## Mock data in $`\mathbb{R}^3`$

Three Gaussian classes are generated with sample sizes 9, 10, and 11.
Their means are

``` math

\mu_A=(-1.5,0.2,0.8),\quad
\mu_B=(0.6,1.7,-0.5),\quad
\mu_C=(1.4,-1.0,1.2).
```

Each class has a different positive-definite covariance matrix:

``` math

\Sigma_A =
\begin{pmatrix}
1.20&0.35&0.10\\
0.35&0.80&-0.15\\
0.10&-0.15&0.65
\end{pmatrix},\quad
\Sigma_B =
\begin{pmatrix}
0.75&-0.20&0.18\\
-0.20&1.35&0.30\\
0.18&0.30&0.90
\end{pmatrix},
```

``` math

\Sigma_C =
\begin{pmatrix}
1.05&0.25&-0.22\\
0.25&0.70&0.12\\
-0.22&0.12&1.25
\end{pmatrix}.
```

With a fixed seed, standard-normal matrices $`Z_g`$ are transformed as

``` math

X_g = Z_g\,\operatorname{chol}(\Sigma_g) + \mu_g.
```

The unequal class sizes test the empirical priors and the sample-size
vector used by the joint projection. Three dimensions keep brute-force
permutation search small and deterministic.

## Covariance projection

For every reference projection, the test constructs a `gips` object
directly:

``` r
search <- gips::gips(
  empirical_covariances,
  sample_sizes,
  was_mean_estimated = TRUE
)
search <- gips::find_MAP(
  search,
  optimizer = "BF",
  show_progress_bar = FALSE
)
permutation <- search[[1L]]
projected <- lapply(
  empirical_covariances,
  gips::project_matrix,
  permutation
)
```

When only one covariance matrix is projected, that matrix is passed
directly to
[`gips::gips()`](https://przechoj.github.io/gips/reference/gips.html).
The brute-force optimizer enumerates the permutations in
$`\mathbb{R}^3`$, avoiding Monte Carlo variation.

The projected covariance $`S`$ is regularized only if its eigenvalue
$`\lambda`$ nearest zero satisfies $`|\lambda| < 0.05`$. In that case,

``` math

s=\frac{0.05-\lambda}{1-0.05},
\qquad
S^*=\frac{S+sI}{1+s}.
```

This transformation is implemented explicitly in the test.

## LDA reconstruction

Let $`M_g`$ denote the sample mean of class $`g`$, and let
$`E_i=x_i-M_{g_i}`$ be the within-class residuals. Define

``` math

D=\operatorname{diag}\left(
\frac{1}{\sqrt{\operatorname{diag}(\operatorname{var}(E))}}
\right).
```

For $`n`$ observations and $`G`$ classes, the covariance supplied to
`gips` is calculated manually as

``` math

S_W=\frac{n}{n-G}\operatorname{cov}(ED).
```

After MAP projection and regularization, write
$`S_W^*=V\operatorname{diag}(d)V^\mathsf{T}`$. The whitening transform
is

``` math

W=DV\operatorname{diag}(d^{-1/2}).
```

Using the empirical priors $`\pi_g=n_g/n`$ and
$`\bar{x}=\sum_g\pi_gM_g`$, the between-class matrix is

``` math

B_g=
\sqrt{\frac{n\pi_g}{G-1}}\,(M_g-\bar{x})W.
```

If $`B=U_B\operatorname{diag}(d_B)V_B^\mathsf{T}`$, the stored
discriminant scaling is $`WV_B`$, truncated using the same numerical
rank rule as the public model. The test independently compares:

- empirical priors, class counts, means, labels, and sample size;
- the MAP permutation;
- between-class singular values;
- the discriminant scaling after resolving arbitrary SVD column signs.

## Separate QDA reconstruction

For each class, the sample covariance

``` math

S_g=\operatorname{cov}(X_g)
```

is projected in a separate brute-force `gips` search. The total sample
size $`n`$ is supplied to each search, matching the estimator’s
statistical parameterization. If

``` math

S_g^*=V_g\operatorname{diag}(d_g)V_g^\mathsf{T},
```

the reference parameters are

``` math

L_g=V_g\operatorname{diag}(d_g^{-1/2}),
\qquad
\log|S_g^*|=\sum_j\log d_{gj}.
```

Because singular vectors can change sign without changing the model, the
test compares $`L_gL_g^\mathsf{T}`$, the precision matrix used by QDA,
rather than the signs of individual columns.

## Joint QDA reconstruction

The joint model starts from the same class covariances but passes
$`(S_A,S_B,S_C)`$ and the sample-size vector $`(9,10,11)`$ to one `gips`
search. Its single MAP permutation is then used to project every
covariance. The scaling matrices and log determinants are calculated
from the resulting matrices using the QDA equations above.

This distinguishes
[`gipsmultqda()`](https://antonikingston.github.io/gipsDA/reference/gipsmultqda.md)
from
[`gipsqda()`](https://antonikingston.github.io/gipsDA/reference/gipsqda.md):
the former estimates one shared symmetry, while the latter estimates
each class symmetry independently.

## Assertions

The three public fitting functions are called only after all reference
parameters have been computed. Numerical quantities are compared with
tolerance $`10^{-10}`$. Discrete quantities, dimensions, labels, sample
size, and MAP permutations are compared exactly. The QDA scaling checks
use precision matrices, and the LDA scaling check aligns SVD signs
before applying the tolerance.

The test therefore validates the complete set of fitted numerical
parameters, apart from the recorded call, against a calculation that is
independent of the package implementation.
