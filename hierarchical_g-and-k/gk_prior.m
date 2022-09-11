function out = gk_prior(theta,draw)

% Returns the product of independent priors 
% Input:  - theta, the vector of parameters to be estimated
% Output: - out, the product of the prior distributions set for each parameter.
%                possibily multiplied with a jacobian for transformations
%                from log-parameter to parameter

if draw == 0

alpha = theta(1);
As = theta(2:length(theta));

alpha_prior = unifpdf(alpha,-10,10);
As_prior = normpdf(As,alpha,1);

out = alpha_prior*prod(As_prior) ;
  
elseif draw == 1
    
    numdatasets = 20;
    
    alpha = unifrnd(-10,10);
    As = alpha + randn(numdatasets,1);
    
    out = [alpha,As'];
    
end