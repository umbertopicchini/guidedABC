function out = twisted_prior(theta,draw)

% Input:  - theta, the vector of parameters to be estimated. If draw=0 no sampling. If draw=1 returns sampels from the prior
% Output: - out, the joint prior or samples from the prior

if draw == 0
  % see Li, J., Nott, D. J., Fan, Y., & Sisson, S. A. (2017). Extending approximate Bayesian computation methods to high dimensions via a Gaussian copula model. Computational Statistics & Data Analysis, 106, 77-89.
  b=0.1;
  
  if length(theta)==2
     addterm = 0;
  else
     addterm = sum(theta(3:end).^2);
  end
  out = exp(-theta(1)^2/200 - (theta(2)-b*theta(1)^2+100*b)^2 / 2 - addterm) ;
  
elseif draw == 1
    
    % see Li, J., Nott, D. J., Fan, Y., & Sisson, S. A. (2017). Extending approximate Bayesian computation methods to high dimensions via a Gaussian copula model. Computational Statistics & Data Analysis, 106, 77-89.
    diag_cov = ones(1,5);
    diag_cov(1) = 100;
    theta = mvnrnd(zeros(1,length(diag_cov)),diag_cov);
    b = 0.1;
    theta(2) = theta(2) + b*theta(1)^2 -100*b;
    
    out = theta;
    
end
