function out = lv_prior(theta,draw)

% Returns the product of independent priors 
% Input:  - theta, the vector of parameters to be estimated
% Output: - out, the product of the prior distributions set for each
% parameter.

if draw == 0
    
  logc1 = theta(1);
  logc2 = theta(2);
  logc3 = theta(3);

  logc1_prior      = unifpdf(logc1, -6,2); 
  logc2_prior      = unifpdf(logc2, -6,2); 
  logc3_prior      = unifpdf(logc3, -6,2); 


  out = logc1_prior*logc2_prior*logc3_prior ;

elseif draw == 1
    
   logc1 = unifrnd(-6,2);
   logc2 = unifrnd(-6,2);
   logc3 = unifrnd(-6,2);
   
   out = [logc1,logc2,logc3];
   
end