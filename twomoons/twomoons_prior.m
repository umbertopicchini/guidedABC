function out = twomoons_prior(theta,draw)

% Returns the product of independent priors 
% Input:  - theta, the vector of parameters to be estimated
% Output: - out, the product of the prior distributions set for each parameter.
%                possibily multiplied with a jacobian for transformations
%                from log-parameter to parameter

if draw == 0

t1 = theta(1);
t2 = theta(2);

t1_prior = unifpdf(t1,-1,1);
t2_prior = unifpdf(t2,-1,1);

  out = t1_prior*t2_prior ;
  
elseif draw == 1
    
    t1_prior_rnd = unifrnd(-1,1);
    t2_prior_rnd = unifrnd(-1,1);
    
    out = [t1_prior_rnd,t2_prior_rnd];
    
end
