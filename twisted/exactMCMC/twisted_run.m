
rng(100)

covariates = [];  % no covariates needed for this model

% set data
dimensionality = 5;
yobs = zeros(dimensionality,1);  %!!! MUST BE A VERTICAL VECTOR !!!!
yobs(1) = 10;




rng(200)
                  % t1   t2 
parbase        = [ zeros(1,dimensionality)];    % just some "baseline parameter values", e.g. prior means
parmask        = [ ones(1,dimensionality)];  % 1= this is a parameter to infer; 0 = this is a known constant
% that is sum(parmask) is the number of parameters to infer

thetastart = zeros(1,dimensionality);
thetastart(1) = 10;
MCMCiter = 300000;
MCMC = zeros(MCMCiter,dimensionality);
MCMC(1,:) = thetastart;

theta_old = thetastart;
%::::::::: PRIOR COMPUTATION ::::::::::::::::::::::::::::::
  % see Li, J., Nott, D. J., Fan, Y., & Sisson, S. A. (2017). Extending approximate Bayesian computation methods to high dimensions via a Gaussian copula model. Computational Statistics & Data Analysis, 106, 77-89.
  b=0.1;
  if length(theta_old)==2
     addterm = 0;
  else
     addterm = sum(theta_old(3:end).^2);
  end
  prior_old = exp(-theta_old(1)^2/200 - (theta_old(2)-b*theta_old(1)^2+100*b)^2 / 2 - addterm) ;
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
like_old = mvnpdf(yobs',theta_old,ones(1,dimensionality));

for ii=2:MCMCiter
    theta = mvnrnd(theta_old,0.1^2*ones(1,dimensionality));
    %::::::::: PRIOR COMPUTATION ::::::::::::::::::::::::::::::
  % see Li, J., Nott, D. J., Fan, Y., & Sisson, S. A. (2017). Extending approximate Bayesian computation methods to high dimensions via a Gaussian copula model. Computational Statistics & Data Analysis, 106, 77-89.
     b=0.1;
    if length(theta)==2
       addterm = 0;
    else
       addterm = sum(theta(3:end).^2);
    end
    prior = exp(-theta(1)^2/200 - (theta(2)-b*theta(1)^2+100*b)^2 / 2 - addterm) ;
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    like = mvnpdf(yobs',theta,ones(1,dimensionality));
    ratio = (like/like_old) * (prior/prior_old);
    if rand < min(1,ratio)
        theta_old = theta;
        prior_old = prior;
        like_old = like;
        MCMC(ii,:) = theta;
    else
        MCMC(ii,:) = theta_old;
    end
end
filename=sprintf('MCMC_dimensionality_%d.dat',dimensionality);
save(filename,'MCMC','-ascii')
