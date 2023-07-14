function [AA,theta,pd,alphaF] = sample_theta(proposal_mean,proposal_cov,nu,type_copula,type_margin,nu_marg)
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
alphaF=[];
AA=[];
dim=size(proposal_cov);
dim=dim(1);
pd=struct([]);
switch type_copula
    case 0
        theta=mvnrnd(proposal_mean,proposal_cov);
        return
    case 1
    theta=zeros(1,dim);
    rho=corrcov(proposal_cov);
    AA = copularnd('Gaussian',rho,1);
    case 2
        theta=zeros(1,dim);
         rho=corrcov(proposal_cov);
        AA = copularnd('t',rho,nu,1);
    case 3
            rho=proposal_cov(1,2)/sqrt(proposal_cov(1,1)*proposal_cov(2,2));
             alphaF=copulaparam('frank',rho,'type','kendall');
              AA= copularnd('frank',alphaF,1);
end           
switch type_margin
    case 1
         for i=1:dim
         a = makedist('Triangular','a',proposal_mean(i)-sqrt(6*proposal_cov(i,i)),'b',proposal_mean(i),'c',proposal_mean(i)+sqrt(6*proposal_cov(i,i)));
         pd=[pd,a];
        theta(i)=icdf(a,AA(i));
         end
    case 2
           for i=1:dim
           a = makedist('tLocationScale','mu',proposal_mean(i),'sigma',sqrt((nu_marg-2)*proposal_cov(i,i)/nu_marg),'nu',nu_marg);
           theta(i)=icdf(a,AA(i));
           pd=[pd,a];
           end
    case 3
            for i=1:dim
            a = makedist('Logistic','mu',proposal_mean(i),'sigma',sqrt(3*proposal_cov(i,i))/pi);
            theta(i)=icdf(a,AA(i));
            pd=[pd,a];
            end
    case 4
        doubleeulergamma=double(eulergamma);
        for i=1:dim
             bet=sqrt(proposal_cov(i,i)*6)/pi;
             % noticen the MATLAB parameterization of the Gumble distributions is slightly defferent from the canonical one found eg in R or Wikipedia (below we use bet instead of -bet)
            a = makedist('ExtremeValue','mu', proposal_mean(i)+bet* doubleeulergamma,'sigma',bet);
            theta(i)=icdf(a,AA(i));
            pd=[pd,a];
         end
    case 5
        for i=1:dim 
                  theta(i)=unifinv(AA(i),proposal_mean(i)-sqrt(3*proposal_cov(i,i)),proposal_mean(i)+sqrt(3*proposal_cov(i,i)));
              end
    case 6
                   for i=1:dim
                       theta(i)=norminv(AA(i),proposal_mean(i),sqrt(proposal_cov(i,i)));
                   end
 end                   
