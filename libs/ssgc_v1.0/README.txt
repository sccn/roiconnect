This small toolbox provides an efficient Matlab implementation for state space
Granger causality computation, as presented in ref. [1] (referred to simply as
"the reference article" in the in-file documentation). See also ref. [2].

The code is provided free of charge and is neither exhaustively tested nor par-
ticularly well documented. The authors accept no liability for its use. You may
use this code in any way you see fit; we ask only that you acknowledge author-
ship and cite ref. [1]. This functionality will also feature as the new default
GC calculation method in the next major release of the MVGC Multivariate Granger
Causality Matlab Toolbox (ref. [3]: http://www.sussex.ac.uk/sackler/mvgc/).

To get started, we recommend that you run and work through the demonstration
script ssgc_demo.m.

Lionel Barnett and Anil K. Seth, May 2015.


Demonstration script
--------------------

ssgc_demo - Generate random state space model and calculate time- and
            frequency-domain causalities and causal graphs

Main computational routines
---------------------------

iss_GC    - Compute time-domain (conditional) GC for a state space model from
            innovations form parameters

iss_SGC   - Compute frequency-domain (conditional) GC for a state space model
            from innovations form parameters

iss_PWGC  - Compute pairwise-conditional time-domain GCs for a state space model
            from innovations form parameters

iss_SPWGC - Compute pairwise-conditional frequency-domain GCs for a state space
            model from innovations form parameters

Helper routines
---------------

iss_MA    - Compute moving-average representation (transfer function) for a
            state space model in innovations form

iss_AR    - Compute autoregressive representation (inverse transfer function)
            for a state space model in innovations form

iss_cpsd  - Calculate cross-power spectral density from transfer function and
            innovations covariance matrix

ss2iss    - Compute innovations form parameters for a general state space model
            by solution of a discrete algebraic Riccati equation (DARE)

parcovar  - Compute partial covariance matrix

Utility routines
----------------

specnorm  - Calculate spectral norm for a vector autoregressive model

cov_rand  - Generate random positive-definite covariance matrix (for testing)

iss_rand  - Generate random innovations form parameters (for testing)

iss_gen   - Generate a realization of a state space process with Gaussian
            innovations from innovations form parameters (for testing)

iss_ds    - Calculate innovations form state space parameters for a downsampled
            state space model

s4sid_CCA - Estimate an innovations form state space model from an empirical
            observation time series using Larimore's Canonical Correlations
            Analysis (CCA) state space-subspace algorithm

ar_rand   - Generate random autoregression coefficients (for testing)

ar_gen    - Generate a realization of a vector autoregressive process with
            Gaussian residuals

ar_IC     - Compute Akaike (AIC) and Schwarz' Bayesian (BIC) information
            criteria for autoregressive model order estimation

arsid_OLS - Estimate a vector autoregressive model from an empirical observation
            time series using Ordinary Least Squares

ar2iss    - Return innovations form state space model parameters corresponding
            to a vector autoregressive model


Notes
-----

1) No additional Matlab toolboxes are required to run this code.

2) Spectral (frequency-domain) quantities are computed for a user-supplied
   vector of points z on the unit circle in the complex plane. Thus, e.g., to
   obtain evenly-spaced z at a frequency resolution fres in the range [0,pi],
   where pi represents the Nyqvist frequency in the normalised frequency range
   [0,2*pi), the following code may be used:

   omega = pi*((0:fres)/fres); % normalised frequencies up to Nyqvist...
   z = exp(-1i*omega);         % ...on unit circle

3) The reason for the proliferation of seemingly pointless Cholesky decompo-
   sitions ("matrix square roots") is to ensure that covariance (resp. spectral)
   matrices are symmetric (resp. Hermitian), and positive-definite.

4) "Pairwise-conditional" Granger causalities for a set of n variables, refers
   to the "causal graph" matrix of all causalities F(i,j) from source variable
   j -> target variable i (i ~= j), conditioned on the remaining n-2 variables.
   This matrix may be significance-tested to obtain a directed functional
   connectivity graph between the n variables. It is also defined in the
   frequency domain. Since directed functional connectivity is frequently a main
   motivation for Granger-causal analysis, especially in the neurosciences, we
   have provided the optimised routines iss_PWGC and iss_SPWGC to calculate
   these quantities more efficiently than by separate computations for each pair
   of variables.

5) We have included an implementation in s4sid_CCA of Larimore's Canonical
   Correlations Analysis (CCA) state space-subspace system identification
   algorithm [4] to estimate a state space model in innovations form from
   empirical time series data. State space model order may optionally be
   estimated using Bauer's Singular Value Criterion (SVC) [5]. There are, of
   course, many other toolboxes and libraries out there for state space model
   identification. The utility iss_gen may be used to generate test data.

6) We also include the utility function ar2iss and a simple OLS vector auto-
   regression system identification routine arsid_OLS, so that this code may be
   used with estimated autoregressive models; but note that if the data has a
   moving average component - i.e. is ARMA - then a pure autoregressive model
   will not be parsimonious, and state space model estimation is preferable. The
   utility ar_gen may be used to generate test data.

7) For statistical inference, eq. 17 of ref. [1] and the discussion following
   suggest that, from the standard large sample theory, maximum likelihood
   estimators for (time domain) state space Granger causalities should be asymp-
   totically chi^2 with m*n1*n2 degrees of freedom, where m is the state space
   dimension, n1 the dimension of the target variable and n2 the dimension of
   the source variable. However, preliminary tests by the authors have not fully
   confirmed this, and surrogate methods are advisable (in the frequency domain
   there are in any case no known asymptotic distributions for GC estimators,
   even in the pure AR case).


References
----------

[1] L. Barnett and A. K. Seth, Granger causality for state-space models, Phys.
    Rev. E 91(4) Rapid Communication, 2015.

[2] V. Solo, State space methods for Granger-Geweke causality measures,
    arXiv:1501.04663.

[3] L. Barnett and A. K. Seth, The MVGC Multivariate Granger Causality Toolbox:
    A new approach to Granger-causal inference, J. Neurosci. Methods 223, 2014.

[4] W. E. Larimore, System identification, reduced-order filtering and modeling
    via canonical variate analysis, in Proceedings of the 1983 American Control
    Conference, IEEE Service Center, Piscataway, NJ, USA, vol. 2, 1983.

[5] D. Bauer, Order estimation for subspace methods, Automatica 37, 1561, 2001.


Acknowledgements
----------------

The authors are grateful to the Dr. Mortimer and Theresa Sackler Foundation,
which supports the Sackler Centre for Consciousness Science.


Contact
-------

Lionel Barnett (corresponding author) and Anil K. Seth
Sackler Centre for Consciousness Science
Department of Informatics
University of Sussex
Falmer
Brighton
BN1 9QJ
UK

http://www.sussex.ac.uk/sackler/

http://users.sussex.ac.uk/~lionelb/
http://users.sussex.ac.uk/~anils/

mailto:l.c.barnett@sussex.ac.uk
