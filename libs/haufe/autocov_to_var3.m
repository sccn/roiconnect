%% autocov_to_var
%
% Calculate VAR parameters from autocovariance sequence
%
% <matlab:open('autocov_to_var.m') code>
%
%% Syntax
%
%     [A,SIG] = autocov_to_var(G)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     G          autocovariance sequence
%
% _output_
%
%     A          VAR coefficients matrix
%     SIG        residuals covariance matrix
%
%% Description
%
% Calculates regression coefficients |A| and residuals covariance matrix
% |SIG| from the autocovariance sequence |G| defined as [[ii_acseq.png]]
% by solving the Yule-Walker equations
%
% <<eq_yweqs.png>>
%
% (where  [[ii_Sigma.png]] = |SIG|). For a |q|-lag autocovariance sequence,
% this routine corresponds to an autoregression of |q| lags. It also
% effects an efficient spectral factorisation if called with the
% autocovariance sequence derived from the cross-power spectral density
% (_e.g._ as calculated by <cpsd_to_autocov.html |cpsd_to_autocov|>).
%
% This routine implements Whittle's recursive LWR algorithm [2] which, for |n|
% variables, performs |2q| separate |n x n| matrix inversions as compared with a
% single |nq x nq| matrix inversion for the conventional "OLS" solution of the
% Yule-Walker equations (see [1]). The LWR algorithm also (unlike OLS)
% guarantees that if the "true" regression model is stable, then the estimated
% model is also stable, even if not of the correct order.
%
% *_Note_*: If the regressions are rank-deficient or ill-conditioned then A may
% be "bad" (i.e. will contain a |NaN| or |Inf|; see <isbad.html |isbad|>) and/or
% warnings will may be issued. The caller should test for both these
% possibilities, by calls to <isbad.html |isbad|> and <warn_supp.html
% |warn_supp|> ... <warn_test.html |warn_test|> respectively. The likely cause
% is that something went wrong in <var_to_autocov.html |var_to_autocov|>, which
% is typically called prior to this function; check the results of the latter
% with a <var_info.html |var_info|> call.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <https://urldefense.proofpoint.com/v2/url?u=http-3A__www.sciencedirect.com_science_article_pii_S0165027013003701&d=DwIGaQ&c=-35OiAkTchMrZOngvJPOeA&r=q2b2jk6iMqHvTzFUUrmWpFcLBpvSfdmimXkepsdyNwg&m=YCA0SzyGKYQTdDj_ywpYQmW8voUEJ-6jnGyIBmJvsrBHYyZFviLcf4kg1HT0q5_V&s=9HqgpUA5vX3_GdsIWjPg-B6EcSil2XSCYodzROIUowM&e=  The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] P. Whittle, "On the fitting of multivariate autoregressions, and the
% approximate canonical factorization of a spectral density matrix",
% _Biometrika_, 50, 1963.
%
%% See also
%
% <var_to_autocov.html |var_to_autocov|> |
% <cpsd_to_autocov.html |cpsd_to_autocov|> |
% <warn_supp.html |warn_supp|> |
% <warn_test.html |warn_test|> |
% <isbad.html |isbad|> |
% <var_info.html |var_info|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [AF,SIG] = autocov_to_var3(G)

[m, ~, n] = size(G);

[AF,SIG] = autocov_to_var(G);

sing = 0;
n_ = n - 1;
while ~isposdef(SIG)  
  sing = 1; 
  [~,SIG] = autocov_to_var(G(:, :, 1:n_));
%   AF = cat(3, AF, zeros(m, m, n-n_));
  n_ = n_ - 1;
  if n_ == 1
    break;
  end
end

if sing
  disp('krass')
end



