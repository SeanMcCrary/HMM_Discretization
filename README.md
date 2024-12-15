# HMM_Discretization
Codes accompanying "Finite-State Markov Chain Approximations: A Hidden Markov Approach"

This folder contains Matlab functions that can be used to discretize continuous Markov processes as in:
“Finite-State Markov-Chain Approximations: A Hidden Markov Approach” by Eva F. Janssens and Sean McCrary.
https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4137592. Please cite appropriately. 

______________________________________________________________________________________________________________________

This folder contains functions that can be used to discretize the four standard processes from the paper:

* AR(1) process
* AR(1)-Mixture process: this is an AR(1) where the innovations are drawn from a Gaussian mixture.
* AR(1)-SV process: note there are multiple formulations for the AR(1)-SV, we have tailored this function to the AR(1)-SV as in the paper. If you want a different formulation, you will have to make some changes to the code or use StationaryEM.
* VAR(1) process: we have made these code for a VAR(1) with TWO variables. If you want to discretize a higher-dimensional process, you need to use StationaryEM and provide your own simulation and initialization.

The HMM discretization method takes as input a simulated sequence from the process to be discretized, an initial grid and an initial transition probability matrix. For these four processes, we have built in a good initial guess, and we optimize the length of the simulated sequence. Hence, all you need to provide to the functions is the parameters of the stochastic process, and the number of grid points of your desired discretization. 

For the AR(1) process, DiscretizeAR1.m takes as input:
* m   - number of grid points
* rho - persistence of the AR(1)
* sig - standard deviation of the shock innovations
  
The function outputs:
* grid - the grid of the discretization
* tpm - the transition probability matrix
* ll     - the average log likelihood: this is a useful statistic when determining how many grid points to use, you can make a scree plot as in the paper, and choose m such that the likelihood has flattened.
* sigs - the standard deviation of the HMM errors (interpret as approximation error standard deviation)
* Y     - the simulated sequence 

For the AR(1)-M process (Gaussian mixture), DiscretizeAR1M.m takes as input:
* m         - number of grid points 
* rho       - persistence 
* mixp     - vector of mixture weights 
* mixmu  - vector of mixture means 
* mixsig  - vector of mixture standard deviations 

 The function outputs: 
 * grid - the grid of the discretization
 * tpm - the transition probability matrix
 * ll     - the average log likelihood: this is a useful statistic when determining how many grid points to use, you can make a scree plot as in the paper, and choose m such that the likelihood has flattened.
 * sigs - the standard deviation of the HMM errors (interpret as approximation error standard deviation)
 * Y     - the simulated sequence used for the discretization method

For the AR(1)-SV process, DiscretizeAR1SV takes as input:
* m        - total number of grid points (note: we do not use tensor grids, so you tell the function how many grid points to use in total)
* rhoM  - persistence of levels 
* rhoV  - persistence of volatility
* omega - s.d. of volatility shocks
* eta   - mean volatility

The function outputs: 
* gridz - the grid of the volatility dimension
* gridy - the grid for the levels
* tpm - the transition probability matrix
* ll     - the average log likelihood: this is a useful statistic when determining how many grid points to use, you can make a scree plot as in the paper, and choose m such that the likelihood has flattened.
* sigs - a vector with the standard deviation of the HMM errors (interpret as approximation error standard deviation): sigs  = [sigY sigZ] given the process is multivariate we have two of these.
* Y     - the simulated sequence used for the discretization method for the levels
* Z      the simulated sequence for the volatility

For the bivariate VAR(1) process, DiscretizeVAR2.m takes as input:
* m        - total number of grid points (note: we do not use tensor grids, so you tell the function how many grid points to use in total)
* A         - coefficient of lag matrix 
* Sigma - variance-covariance matrix of errors   

The function outputs: 
* grid 	- the grid (two columns with m rows)
* tpm 	- the transition probability matrix
* ll     	- the average log likelihood: this is a useful statistic when determining how many grid points to use, you can make a scree plot as in the paper, and choose m such that the likelihood has flattened.
* sigs 	- a vector with the standard deviation of the HMM errors (interpret as approximation error standard deviation): sigs  = [sigY1 sigY2] given the process is multivariate we have two of these.
* Y1     - the simulated sequence for the first dimension
* Y2     - the simulated sequence for the second dimension

______________________________________________________________________________________________________________________

This folder contains functions that can be used to discretize the two lifecycle processes from Section 4 
“Finite-State Markov-Chain Approximations: A Hidden Markov Approach”: 

Simulated data used as an input to these functions can be found here: 
https://www.dropbox.com/scl/fo/tuw1535ckoz1u8m539g3b/AEcW9KQzMuHCUFUe4dn5Y00?rlkey=u815lzmy1g8jl8pmfmiy8dshh&st=kh3l6omh&dl=0

* DiscretizeGKOS.m returns an age-dependent discretization for the process of Guvenen, Karahan, Ozkan and Song (2021) as discussed in the paper.
* DiscretizeABB.m returns an age-dependent discretization for the process of Arellano, Blundell and Bonhomme (2017).

The function only requires you to input “m”, the number of grid points. Everything else happens inside the function, that is, it loads in a simulated data set, builds an initial guess and does the discretization. Note that the GKOS process is actually two-dimensional, because of the non-employment shocks in the process, hence we allow users to both pick m1, the numbers of zeros, and m2, the number of non-zero earnings states. 

The functions both have the same structure, where the input is “m”, the number of grid points, for ABB, and m1,m2 for GKOS, and the output is “out”, a structure containing:

* out.Grid   - grid : T x m grid
* out.Pi       - the transition probability matrix: (T-1) cells of m x m 
* out.initial  - the initial distribution over the states at time 0: vector size m

We also provide the simulated data as output in the structure. 


______________________________________________________________________________________________________________________

This folder contains functions that can be used to discretize a general procces with StationaryEM.m 

The HMM discretization method takes as input a simulated sequence from the process to be discretized, an initial grid and an initial transition probability matrix. 

For stationary processes, this means the following:
* You need to simulate a (long T) sequence from the process you want to discretize. For better performance, you may want to generate multiple (N) sequences. If the process is multidimensional (k>1), the function expects the following format: Y is T x N x k. We strongly recommend NxT to be at least 10,000. 
* The initial grid mu will be of size m x k, where m is the number of grid points, and k the dimension of your process. For univariate processes, we recommend an equal-spaced grid of 3 standard deviations as the initial grid. For multivariate processes, you can either use equal-spaced grids in each dimension, and take a tensor grid as initial, or, better, would be to apply k-means to your simulated data set and use the centroids of kmeans as an initial guess. For an example, you can look at the VAR(1) codes in the “Standard Processes” folder. 
* The initial transition probability matrix is size m x m, which we recommend you to compute using a simple binning method. In the case that you use kmeans to obtain an initial grid, you can use pdist2 to compute a bin-assignment, and you can then compute an initialization from this. 

The function to call is:

[Pi,Grid,ll,sige] = StationaryEM(Ysim,Pi0,Grid0)

Arguments
* Y:     an T-N-k array of data
* Pi0:   an m-m initial guess of transitions
* Grid0: an m-k initial guess of grid points

Output
* Pi:   transition probability matrix
* Grid: matrix of estimated grid points
* ll:   log-likelihood
* sige: vector of standard deviations of the error terms



