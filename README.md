# HMM_Discretization

Codes accompanying *"Finite-State Markov Chain Approximations: A Hidden Markov Approach"*

This folder contains MATLAB functions that can be used to discretize continuous Markov processes as in:  
**“Finite-State Markov-Chain Approximations: A Hidden Markov Approach” by Eva F. Janssens and Sean McCrary**  
https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4137592.  
**Please cite appropriately.**

---

## Standard Processes

This folder contains functions that can be used to discretize the four standard processes from the paper:

- **AR(1) process**
- **AR(1)-Mixture process**: this is an AR(1) where the innovations are drawn from a Gaussian mixture.
- **AR(1)-SV process**: note there are multiple formulations for the AR(1)-SV. We have tailored this function to the AR(1)-SV as in the paper. If you want a different formulation, you will have to make some changes to the code or use `StationaryEM`.
- **VAR(1) process**: we have made these codes for a VAR(1) with TWO variables. If you want to discretize a higher-dimensional process, you need to use `StationaryEM` and provide your own simulation and initialization.

The HMM discretization method takes as input a simulated sequence from the process to be discretized, an initial grid, and an initial transition probability matrix. For these four processes, we have built in a good initial guess, and we optimize the length of the simulated sequence. Hence, all you need to provide to the functions is the parameters of the stochastic process, and the number of grid points of your desired discretization.

---

### DiscretizeAR1.m

**Inputs**:
- `m` – number of grid points  
- `rho` – persistence of the AR(1)  
- `sig` – standard deviation of the shock innovations

**Outputs**:
- `grid` – the grid of the discretization  
- `tpm` – the transition probability matrix  
- `ll` – the average log likelihood: this is a useful statistic when determining how many grid points to use. You can make a scree plot as in the paper, and choose `m` such that the likelihood has flattened.  
- `sigs` – the standard deviation of the HMM errors (interpreted as approximation error standard deviation)  
- `Y` – the simulated sequence

---

### DiscretizeAR1M.m

**Inputs**:
- `m` – number of grid points  
- `rho` – persistence  
- `mixp` – vector of mixture weights  
- `mixmu` – vector of mixture means  
- `mixsig` – vector of mixture standard deviations

**Outputs**:
- Same as for `DiscretizeAR1.m`

---

### DiscretizeAR1SV.m

**Inputs**:
- `m` – total number of grid points (note: we do not use tensor grids, so you tell the function how many grid points to use in total)  
- `rhoM` – persistence of levels  
- `rhoV` – persistence of volatility  
- `omega` – standard deviation of volatility shocks  
- `eta` – mean volatility

**Outputs**:
- `gridz` – the grid of the volatility dimension  
- `gridy` – the grid for the levels  
- `tpm` – the transition probability matrix  
- `ll` – the average log likelihood  
- `sigs` – a vector with the standard deviation of the HMM errors (interpreted as approximation error standard deviation): `sigs = [sigY sigZ]`. Given the process is multivariate, we have two of these.  
- `Y` – the simulated sequence used for the discretization method for the levels  
- `Z` – the simulated sequence for the volatility

---

### DiscretizeVAR2.m

**Inputs**:
- `m` – total number of grid points (note: we do not use tensor grids, so you tell the function how many grid points to use in total)  
- `A` – coefficient of lag matrix  
- `Sigma` – variance-covariance matrix of errors

**Outputs**:
- `grid` – the grid (two columns with `m` rows)  
- `tpm` – the transition probability matrix  
- `ll` – the average log likelihood  
- `sigs` – a vector with the standard deviation of the HMM errors: `sigs = [sigY1 sigY2]`  
- `Y1`, `Y2` – the simulated sequences for the first and second dimension

---

## Lifecycle Processes

This folder contains functions that can be used to discretize the two lifecycle processes from Section 4 of the paper:  
**“Finite-State Markov-Chain Approximations: A Hidden Markov Approach”**

Simulated data used as input to these functions can be found here:  
https://www.dropbox.com/scl/fo/tuw1535ckoz1u8m539g3b/AEcW9KQzMuHCUFUe4dn5Y00?rlkey=u815lzmy1g8jl8pmfmiy8dshh&st=kh3l6omh&dl=0

### DiscretizeGKOS.m

Returns an **age-dependent** discretization for the process of Guvenen, Karahan, Ozkan and Song (2021), as discussed in the paper.

> Note: the GKOS process is actually two-dimensional because of the non-employment shocks in the process. Hence we allow users to pick both `m1`, the number of zeros, and `m2`, the number of non-zero earnings states.

### DiscretizeABB.m

Returns an **age-dependent** discretization for the process of Arellano, Blundell and Bonhomme (2017).

**Inputs**:
- `m` (ABB), and `m1, m2` for GKOS

**Outputs**:
Both functions have the same structure. The input is the number of grid points, and the output is a structure `out` containing:

- `out.Grid` – grid: T × m  
- `out.Pi` – transition matrix: (T−1) cells of m × m  
- `out.initial` – the initial distribution over the states at time 0 (vector of length m)  
- Simulated data as additional fields in `out`

Everything else happens inside the function: it loads in a simulated data set, builds an initial guess, and does the discretization.

---

## General Process Discretization: StationaryEM.m

This folder contains a general-purpose function, `StationaryEM.m`, which can be used to discretize **any stationary stochastic process**.

The HMM discretization method takes as input a simulated sequence from the process to be discretized, an initial grid, and an initial transition probability matrix.

### Requirements:

For **stationary processes**, this means the following:

- You need to simulate a (long T) sequence from the process you want to discretize. For better performance, you may want to generate multiple (N) sequences.  
- If the process is multidimensional (k > 1), the function expects the following format: `Y` is T × N × k.  
  **We strongly recommend T × N ≥ 10,000.**

- The initial grid `Grid0` should be of size m × k:
  - For **univariate processes**, we recommend an equally spaced grid over ±3 standard deviations.
  - For **multivariate processes**, either:
    - Construct tensor-product grids over each dimension, or  
    - Apply `kmeans` to the simulated dataset and use the centroids as initial grid points (recommended).

- The initial transition matrix `Pi0` should be size m × m:
  - If using `kmeans`, use `pdist2` to assign bins, and compute a bin-based initialization from this.

### Function Call:

```matlab
[Pi, Grid, ll, sige] = StationaryEM(Ysim, Pi0, Grid0)




