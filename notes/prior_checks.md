# Checks and error messages for misspecified priors in `MCMCglmm`

This note surveys how `MCMCglmm()` responds to common mistakes in the `prior` argument. Line numbers refer to commit `db07ed7`, which added structure labels and expected/observed dimensions to the error messages in `priorformat()`.

All the examples use the built-in `PlodiaPO` data, with this setup:

```r
library(MCMCglmm)
data(PlodiaPO)
f <- function(prior, random = ~ FSfamily + plate) {
  MCMCglmm(PO ~ 1, random = random, data = PlodiaPO, prior = prior,
           verbose = FALSE, nitt = 20, burnin = 0, thin = 1)
}
R <- list(V = 1, nu = 0.002)
G <- list(G1 = list(V = 1, nu = 1), G2 = list(V = 1, nu = 1))
```

## Summary

| Category | Mistake | Result |
|---|---|---|
| Silent | Misspelled element in a G/R prior (e.g. `alpha.v`) | No error or warning; the element is ignored |
| Silent | G prior given without an R prior | No error; the R prior falls back to the default `V = 1, nu = 0` (improper) |
| Cryptic | Only `mu`, or only `V`, given in `prior$B` | `'data' must be of a vector type, was 'NULL'` |
| Cryptic | `fix` larger than the dimension of `V` | `subscript out of bounds` |
| Cryptic | Negative `nu` | `ill-conditioned cross-product: can't form Cholesky factor` |
| Cryptic | `nu = NA` | `NA/NaN/Inf in foreign function call (arg 24)` |
| Cryptic | `V` given as a string | `infinite or missing values in 'x'` |
| Could say more | Wrong number of G structures, or `G` not given as a list of lists | `prior$G has the wrong number of structures` |
| Could say more | `fix` out of range | `fix term in priorG must be at least one less than the dimension of V` |
| Could say more | Unknown top-level element of `prior` (e.g. `b`) | `prior list should contain elements R, G, and/or B and/or S only` |
| Could say more | Fixed-effect `V` not positive definite | `fixed effect V prior is not positive definite` |
| Informative | Missing `V` or `nu`; `V` not positive definite; dimension or length mismatches | Names the structure and term, and gives expected/observed dimensions |

## Silent failures

These are the most serious problems, because the model runs with a prior the user didn't intend.

### Misspelled element names inside a G/R prior

```r
f(list(R = R, G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.v = 1000),
                       G2 = list(V = 1, nu = 1))))
## runs without error or warning; no parameter expansion is used
```

`priorformat()` checks the element names against a list of valid ones and builds a message, but never passes it to `stop()` or `warning()` (`R/priorformat.R:71-74`):

```r
if(any(names(prior)%in%pnames==FALSE)){
   paste(paste(names(prior)[which(names(prior)%in%pnames==FALSE)], sep=" "), " are not valid prior specifications for G/R-structures")
}
```

A misspelled `nu` (e.g. `Nu`) is caught indirectly by "nu not specified". Misspellings of optional elements (`alpha.mu`, `alpha.V`, `fix`, `beta.mu`, `beta.V`, `covu`) are not caught.

**Suggested fix:** use `warning()` rather than `stop()`, so that existing code that passes extra elements doesn't start failing. The message should also name the structure and term, as the other `priorformat()` messages now do.

### G prior without an R prior

```r
m <- MCMCglmm(PO ~ 1, random = ~ FSfamily, data = PlodiaPO, verbose = FALSE,
              nitt = 200, burnin = 0, thin = 1,
              prior = list(G = list(G1 = list(V = 0.5, fix = 1))))
range(m$VCV[, "FSfamily"])
## [1] 0.5 0.5    (the G prior is used)
```

There is a check meant to catch this (`R/MCMCglmm.R:103`, "either both or neither R and G structures need a prior"). But the code that converts old-style priors runs first (`R/MCMCglmm.R:72-76`). It turns a missing `prior$R` into `list(R1 = NULL)`, so `is.null(prior$R)` is `FALSE` and the check never fires. `priorformat()` then gets `NULL` for the R structure and uses its default, `V = 1, nu = 0`.

**Suggested fix:** skip the old-style conversion when `prior$R` is `NULL`, or run the both-or-neither check before the conversion.

## Cryptic errors

These stop the fit, but the message comes from deep inside the code and doesn't point to the problem.

| Mistake | Example | Message |
|---|---|---|
| `prior$B` has `V` but no `mu` | `B = list(V = 1e3)` | `Error in matrix(prior$B$mu, length(prior$B$mu), 1): 'data' must be of a vector type, was 'NULL'` |
| `prior$B` has `mu` but no `V` | `B = list(mu = 0)` | `Error in array(x, ...): 'data' must be of a vector type, was 'NULL'` |
| `fix` larger than the dimension of `V` | `G1 = list(V = 1, nu = 1, fix = 2)` | `Error in prior$V[prior$fix:nrow(prior$V), prior$fix:nrow(prior$V)]: subscript out of bounds` |
| Negative `nu` | `G1 = list(V = 1, nu = -5)` | `ill-conditioned cross-product: can't form Cholesky factor` |
| `nu = NA` | `G1 = list(V = 1, nu = NA)` | `NA/NaN/Inf in foreign function call (arg 24)` |
| Non-numeric `V` | `G1 = list(V = "1", nu = 1)` | `Error in eigen(m, ...): infinite or missing values in 'x'` |

**Suggested fixes:**

- **`prior$B`:** check that both `mu` and `V` are present (`R/MCMCglmm.R:1104-1113`).
- **`fix`:** check its range in `priorformat()` before `prior$V` is subset (`R/priorformat.R:141`), naming the structure and giving the allowed range.
- **`nu`:** check in `priorformat()` that it is numeric, not `NA`, and not negative.
- **`V`:** check in `priorformat()` that it is numeric.

## Messages that could say more

- **`prior$G has the wrong number of structures`** (and the same for `start$G`; `R/MCMCglmm.R:658-662`). This doesn't give the expected and actual counts. It is also the message for `G = list(V = 1, nu = 1)` instead of `G = list(G1 = list(V = 1, nu = 1))`, which is probably the most common mistake. The message could give both counts and, when `prior$G` has an element named `V` or `nu`, suggest the nested form.
- **`fix term in priorG must be at least one less than the dimension of V`** (`R/MCMCglmm.R:1392-1396`). This doesn't say which structure. The wording also doesn't match the test, which allows `fix` to equal the dimension of `V`. In practice the "subscript out of bounds" error above fires first, so this check is rarely reached.
- **`prior list should contain elements R, G, and/or B and/or S only`** (`R/MCMCglmm.R:51`). This could name the offending element.
- **`fixed effect V prior is not positive definite`** (`R/MCMCglmm.R:1114`). This could give the dimension of `V`. The dimension-mismatch messages for `prior$B` already report expected and observed dimensions.

## Messages that are already informative

After commit `db07ed7`, the checks in `priorformat()` name the structure and the model term, and the dimension checks report expected and observed values. For example:

```r
f(list(R = R, G = list(G1 = list(V = 1, nu = 1), G2 = list(V = diag(2), nu = 1))))
## V is the wrong dimension for prior$G[[2]] (term 'plate') (expected 1 x 1, observed 2 x 2)
f(list(R = list(V = 1, Nu = 0.002), G = G))
## nu not specified for prior$R[[1]] (term 'units')
f(list(R = R, G = list(G1 = list(V = -1, nu = 1), G2 = list(V = 1, nu = 1))))
## V is not positive definite for prior$G[[1]] (term 'FSfamily')
```

This covers:

- missing `V` or `nu`;
- `V`, `alpha.V` or `beta.V` not positive definite;
- wrong dimension or length of `V`, `alpha.V`, `alpha.mu`, `beta.mu` or `beta.V`;
- `nu` not a scalar;
- the same checks on starting values in `start$G`/`start$R`;
- unsupported combinations of parameter expansion with residual, `antev` or `idvm` structures.

## Suggested order of work

1. Warn about misspelled element names in G/R priors.
2. Make the both-or-neither check for G and R priors work.
3. Check that `prior$B` has both `mu` and `V`.
4. Check the range of `fix` in `priorformat()`.
5. Check that `nu` is numeric, not `NA`, and not negative.
6. Give expected and actual counts, plus the hint about the nested form, in "wrong number of structures".
