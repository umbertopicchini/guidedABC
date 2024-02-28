rng(100)  % was 100 for reproducibility

log_x10 = log(50);   % true parameter value used for data simulation
log_x20 = log(100);   % true parameter value used for data simulation
log_c1  = log(1);   % true parameter value used for data simulation
log_c2  = log(0.005);% true parameter value used for data simulation
log_c3  = log(0.6);   % true parameter value used for data simulation

sampletime = [0:1:31];  % the vector of sampling times
problem = 'lv';

bigtheta_true = [log_x10,log_x20,log_c1,log_c2,log_c3];   % store here all parameters needed for SDE simulation

%DEFINING PROPENSITY FUNCTIONS AS A FUNCTION HANDLE IN VECTORIZED MANNER
const_rates = [exp(log_c1) exp(log_c2) exp(log_c3)]; % these are the (c1,c2,c3) constant rates
prop = @(x,const_rates)([const_rates(1).*x(:,1),...
                 const_rates(2).*x(:,1).*x(:,2),...
                 const_rates(3).*x(:,2)]);
stoichiometry = [1 -1 0; 0 1 -1];
xinit = [exp(log_x10);exp(log_x20)];

covariates.sampletime = sampletime;
covariates.stoichiometry = stoichiometry;


%:::: generate data ::::::::::::::::::::::::::::::::::::::::::::::::::::
% or you can load it. It is attached as lv_data.dat
%[t,x] = directMethod(stoichiometry', prop, [sampletime(1),sampletime(end)], xinit', const_rates);
%% subselect the output
%yobs=interp1(t,x,sampletime);
%save('lv_data.dat','yobs','-ascii')
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

yobs = load('lv_data.dat');

%::: ABC PILOT to collect simulated summaries and compute a matrix that weights them

% THE PART BELOW CAN BE SKIPPED BY JUST LOADING THE ATTACHED
% SUMMARIES_PILOT MATRIX, THAT IS YOU CAN JUN RUN
% load('summaries_pilot.dat')

% size_pilot = 5000;
% summaries_pilot = zeros(size_pilot,9);
% abc_distances = zeros(size_pilot,1);
% summobs = lv_abc_summaries(yobs);
% for ii=1:size_pilot
%    ii
%    try
%      logc1_pilot   = unifrnd(-6,2); 
%      logc2_pilot   = unifrnd(-6,2); 
%      logc3_pilot   = unifrnd(-6,2);
%      const_rates_pilot = [exp(logc1_pilot) exp(logc2_pilot) exp(logc3_pilot)]; % these are the (c1,c2,c3) constant rates
%      prop_pilot = @(x,const_rates_pilot)([const_rates_pilot(1).*x(:,1),...
%                 const_rates_pilot(2).*x(:,1).*x(:,2),...
%                 const_rates_pilot(3).*x(:,2)]);
%      [t,x] = directMethod(stoichiometry', prop_pilot, [sampletime(1),sampletime(end)], xinit', const_rates_pilot);
%   catch
%       try
%          logc1_pilot   = unifrnd(-6,2); 
%          logc2_pilot   = unifrnd(-6,2); 
%          logc3_pilot   = unifrnd(-6,2);
%          const_rates_pilot = [exp(logc1_pilot) exp(logc2_pilot) exp(logc3_pilot)]; % these are the (c1,c2,c3) constant rates
%          prop_pilot = @(x,const_rates_pilot)([const_rates_pilot(1).*x(:,1),...
%          const_rates_pilot(2).*x(:,1).*x(:,2),...
%          const_rates_pilot(3).*x(:,2)]);
%          [t,x] = directMethod(stoichiometry', prop_pilot, [sampletime(1),sampletime(end)], xinit', const_rates_pilot);
%       catch
%          logc1_pilot   = unifrnd(-6,2); 
%          logc2_pilot   = unifrnd(-6,2); 
%          logc3_pilot   = unifrnd(-6,2);
%          const_rates_pilot = [exp(logc1_pilot) exp(logc2_pilot) exp(logc3_pilot)]; % these are the (c1,c2,c3) constant rates
%          prop_pilot = @(x,const_rates_pilot)([const_rates_pilot(1).*x(:,1),...
%          const_rates_pilot(2).*x(:,1).*x(:,2),...
%          const_rates_pilot(3).*x(:,2)]);
%          [t,x] = directMethod(stoichiometry', prop_pilot, [sampletime(1),sampletime(end)], xinit', const_rates_pilot);
%       end
%   end
%   % subselect the output
%   try
%       xhat=interp1(t,x,sampletime);
%   catch
%       xhat = NaN*ones(length(sampletime),2);
%   end
%   sim_summaries = lv_abc_summaries(xhat);
%   summaries_pilot(ii,:) = sim_summaries;
% end
% 
% save('summaries_pilot.dat','summaries_pilot','-ascii')
% 
% return
%:::::::::: END OF PILOT :::::::::::::::::::::::::::::::::::::::::::::

load('summaries_pilot.dat')

% remove NaNs
idnan = any(isnan(summaries_pilot),2); % find rows with NaNs
summaries_pilot_fixed = summaries_pilot(~idnan,:);  % remove those
% remove nasty outliers before computing a measure of variation for the
% summaries
%summaries_pilot_fixed = rmoutliers(summaries_pilot_fixed,'percentiles',[1.25,98.75]);

% obtain weighting matrix for ABC summaries
summ_weights = diag(mad(summaries_pilot_fixed,0).^2);

% Possible sampling: blocked, blockedopt, hybrid, fullcond, full_cond_opt, olcm, standard.
sampling = "blockedopt";
% as an example here we have set the guided "blockedopt" method. 
% We also specify it as a non-copula method, see the explanation further
% below
type_copula = 0;
type_marginals=5;  % any other number is fine, as it is not used (since type_copula = 0)
nu= 5; % any other number is fine, as it is not used (since type_copula = 0)
nu_marg=5;  % any other number is fine, as it is not used (since type_copula = 0)

% If you are interested in running any non-copula guided method (i.e. blocked, blockedopt, hybrid, fullcond, fullcondopt) or non-guided method (olcm, standard),
% then set type_copula=0, and in such case type_marginals will not be used and can be fixed to any number, e.g.
% type_marginals =6.
% The entries (nu, type_copula, type_marginals,nu_marginals) are relevant
% only when choosing copula-versions for "blocked", "blockedopt" and "hybrid" as samplers. If
% you are interested in the copula-based versions of them, then choose:
% type_copula = 1 for Gaussian copula, type_copula = 2 for t Copula, and 
% type_marginals = 1 for triangular , 2 for location-scale Student't,
% 3 for logistic, 4 for Gumbel, 5 for uniform, 6 for Gaussian.
% If you are interested in using the location-scale Student't as marginal distributions, you
% need to specify the degree of freedom nu_marg. If you are not interested in using such distribution, set nu_marg to any number,
% as its input is needed but it will not be used.
% If you are interested in using a t copula, you need to specify the degree of freedom nu. If you are not interested in using such copula, 
% set nu to any number, as its input is needed, despite not using it.

folder='blockedopt'
% "folder" is the name of the folder where you want your results to be stored in. It MUST be created before running this script or else it will result in an error.
% The whole path to access the folder can be given as input, allowing you to run the script from any other folder.
if ~exist(folder, 'dir') || ~strcmp(folder,sampling)
    error('A folder to store your results must be created first and must be named exactly as the sampler you wish to use.')
end


%::::::: ABC INFERENCE :::::::::::::::::::::::::::::::::::::::::::::::

rng(200)

                  % log_x10  log_x20      log_c1   log_c2   log_c3  
bigtheta_start = [  log(50)  log(100)    log(1)  log(0.005)  log(0.6)   ];
parmask        = [        0        0         1         1          1     ];

parbase = bigtheta_start;



% ABC settings
numparticles = 2000;  % notice the first SMC-ABC iteration will take several minutes to complete with these many particles
ABCthreshold = 30;  % starting threshold
alpha = 30;  % this is a percentile (\psi in the paper), determining the ABC threshold at each iteration. 
ABCthreshold_final = 3; % stop as soon as the threshold becomes smaller than this value 
                          % WARNING: you may want to enlarge ABCthreshold_final. Most samplers will struggle to reach this value.
numattempt = 1;           % number of independent repetitions of the experiment, on the same data 

tic
for attempt=1:numattempt
   ABCdraws = abc_sequential_varying_thresholds_lotkavolterra(problem,yobs,covariates,parmask,parbase,ABCthreshold,ABCthreshold_final,summ_weights,numparticles,alpha,sampling,attempt,nu,type_copula,type_marginals,nu_marg,folder);
end

eval = toc

