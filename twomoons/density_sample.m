function dens=density_sample(AA,theta,pd,alpha,proposal_mean,proposal_cov,nu,type_copula,type_margin)
%%% Type_copula =0: Gaussian sampling
%%% Type copula = 1 Gaussian copula
%%% Type copula = 2 t copula
%%% Type copula = 3 Frank copula

%%% type_margin = 1: triangular
%%% type_margin = 2  tlocationscale
%%% type_margin = 3 logistic
%%% type_margin = 4 Gumbel
%%% type_margin = 5 uniform 
%%% type_margin = 6 gaussian

% 
%% Gaussian cop+Triangular/tlocationscale is slower than the unif!
%% Gaussian cop+Logistic is faster than the two above, slower than the others though
%% t cop+Triangular/tlocationscale is slower than the unif!
%% t cop+Logistic is faster than the two above, slower than the others though

switch type_copula 
    case 0
    dens=mvnpdf(theta,proposal_mean,proposal_cov);
    return
    case 1
        rho=corrcov(proposal_cov);
    denscop=copulapdf('Gaussian',AA,rho);
    case 2 
        rho=corrcov(proposal_cov);
          denscop=copulapdf('t',AA,rho,nu);
    case 3
                rho=corrcov(proposal_cov);
                denscop=copulapdf('frank',AA,alpha);

end
dim=length(proposal_mean);

switch type_margin
    case 5
                     densmarg= prod(unifpdf(theta,proposal_mean-sqrt(3*diag(proposal_cov)'),proposal_mean+sqrt(3*diag(proposal_cov)')));
    case 6
                     densmarg=prod(normpdf(theta,proposal_mean,sqrt(diag(proposal_cov)')));
    case 1     
        densmarg=1;
        for i=1: dim
         densmarg=densmarg*pdf(pd(i),theta(i));%*pdf(pd2,theta(2)); 
        end
         
    case 2    
        densmarg=1;
        for i=1: dim
         densmarg=densmarg*pdf(pd(i),theta(i));%*pdf(pd2,theta(2)); 
        end
  case 3 
      densmarg=1;
        for i=1: dim
         densmarg=densmarg*pdf(pd(i),theta(i));%*pdf(pd2,theta(2)); 
        end
    case 4     
        densmarg=1;
        for i=1: dim
         densmarg=densmarg*pdf(pd(i),theta(i));%*pdf(pd2,theta(2)); 
    end
end     
   
dens=denscop*densmarg;

              
    
            