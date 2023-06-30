function out = cell_prior(theta,draw)

% Returns the product of independent priors for parameters of a g-and-k distribution
% Input:  - theta, the vector of parameters to be estimated
% Output: - out, the product of the prior distributions set for each parameter.
%                possibily multiplied with a jacobian for transformations
%                from log-parameter to parameter

if draw == 0


Pm = theta(1);
Pp = theta(2);

Pm_prior = unifpdf(Pm,0,1);
Pp_prior = unifpdf(Pp,0,1);

  out = Pm_prior * Pp_prior;
  
elseif draw == 1
    
    Pm_prior_rnd = unifrnd(0,1);
    Pp_prior_rnd = unifrnd(0,1);
    
    out = [Pm_prior_rnd,Pp_prior_rnd];
    
end
