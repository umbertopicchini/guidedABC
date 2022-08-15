function x = twisted_model(bigtheta,~,numsim)

if nargin == 1
    numsim=1;
end

if nargin < 2
    numsim=1;
end

% see Li, J., Nott, D. J., Fan, Y., & Sisson, S. A. (2017). Extending approximate Bayesian computation methods to high dimensions via a Gaussian copula model. Computational Statistics & Data Analysis, 106, 77-89.
x = mvnrnd(bigtheta,ones(1,length(bigtheta)));


x = x(:);

end

