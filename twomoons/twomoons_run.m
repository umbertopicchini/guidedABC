% In the following, you can find instructions on how to run "abc_sequential_fixed_thresholds", allowing to sample from any of the
% samplers in Picchini and Tamborrino (2002) for the case of a vector of prefixed thresholds.

rng(100)

covariates = [];  % no covariates needed for this case-study

% Possible sampling: blocked, blockedopt, hybrid, fullcond, full_cond_opt, olcm, standard.
sampling = "blocked";
% as an example here we have set the guided "blocked" method. 
% We also specify it as a non-copula method, see the explanation further
% below
type_copula = 0;
type_marginals=6;  % any other number is fine, as it is not used (since type_copula = 0)
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

folder='blocked'
% "folder" is the name of the folder where you want your results to be stored in. It MUST be created before running this script or else it will result in an error.
% The whole path to access the folder can be given as input, allowing you to run the script from any other folder.
if ~exist(folder, 'dir') || ~strcmp(folder,sampling)
    error('A folder to store your results must be created first and must be named exactly as the sampler you wish to use.')
end

% set observed data (*MUST BE A VERTICAL VECTOR*)
yobs = [0;0];

%::: ABC PILOT to collect simulated summaries and compute a matrix that weights them
size_pilot = 5000;
summaries_pilot = zeros(size_pilot,2);
abc_distances = zeros(size_pilot,1);
summobs = twomoons_abc_summaries(yobs);
for ii=1:size_pilot
    prior_draw = twomoons_prior([],1);
    t1 = prior_draw(1);
    t2 = prior_draw(2);
    x = twomoons_model([t1,t2],covariates,1);
    sim_summaries = twomoons_abc_summaries(x);
    summaries_pilot(ii,:) = sim_summaries;
end

save('summaries_pilot','summaries_pilot')
save('summaries_pilot.dat','summaries_pilot','-ascii')

% obtain weighting matrix for ABC summaries
summ_weights = diag(mad(summaries_pilot,0).^2);

rng(200)
                  % t1   t2 
parbase        = [ 0      0];    % just some "baseline parameter values", e.g. prior means
parmask        = [ 1      1];  % 1= this is a parameter to infer; 0 = this is a known constant
% that is sum(parmask) is the number of parameters to infer

% ABC settings
numparticles = 1000;
ABCthreshold = [4 3 2 1 0.5 0.4 0.3 0.2 0.1 0.08 0.06 ];  % sequence of prefixed ABC thresholds, one for each iteration
problem = 'twomoons'; % a string specifying the name of a case study
numattempts = 1; % number of independent repetitions of the experiment, on the same data 

for attempt = 1:numattempts 
   ABCdraws = abc_sequential_fixed_thresholds(problem,yobs,covariates,parmask,parbase,ABCthreshold,summ_weights,numparticles,sampling,attempt,nu,type_copula,type_marginals,nu_marg,folder);
end


