% this is the RUN file for the hierarchical g-and-k model

rng(1234) 

numdatasets = 20; % WARNING <-- YOU MUST HARDCODE THIS VALUE IN GK_model and GK_PRIOR TOO

% ground-truth parameters
B_true = rand;
g_true = rand;
k_true = rand;
alpha_true = unifrnd(-10,10); 
As_true = alpha_true + randn(numdatasets,1);
covariates = [];  % no covariates needed for this model

% Possible sampling: blocked, blockedopt, hybrid, fullcond, full_cond_opt, olcm, standard.
sampling = "hybrid";
% as an example here we have set the guided "hybrid" method. 
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

folder='hybrid'
% "folder" is the name of the folder where you want your results to be stored in. It MUST be created before running this script or else it will result in an error.
% The whole path to access the folder can be given as input, allowing you to run the script from any other folder.
if ~exist(folder, 'dir') || ~strcmp(folder,sampling)
    error('A folder to store your results must be created first and must be named exactly as the sampler you wish to use.')
end

% store grund-truth parameters into bigtheta_true
bigtheta_true = [ B_true    g_true    k_true    alpha_true   As_true'    ];
% simulate data
yobs = gk_model(bigtheta_true);  %!!! MUST BE A VERTICAL VECTOR !!!!

                  % B             g    k        alpha         As 
parbase        = [ B_true    g_true    k_true    alpha_true   As_true'    ];    % just some "baseline parameter values", e.g. prior means
parmask        = [      0         0         0      1          ones(1,numdatasets)];  % 1= this is a parameter to infer; 0 = this is a known constant

bigtheta = parbase;

%::: ABC PILOT to collect simulated summaries and compute a matrix that weights them
% size_pilot = 5000;
% abc_distances = zeros(size_pilot,1);
% summobs = gk_abc_summaries(yobs);
% summaries_pilot = zeros(size_pilot,length(summobs));
% for ii=1:size_pilot
%     prior_draw = gk_prior([],1);
%     bigtheta = param_unmask(prior_draw,parmask,parbase); % combine fixed (known) parameters and free parameters
%     x = gk_model(bigtheta,covariates,1);
%     sim_summaries = gk_abc_summaries(x);
%     summaries_pilot(ii,:) = sim_summaries;
% end
% 
% 
% save('summaries_pilot','summaries_pilot')
% save('summaries_pilot.dat','summaries_pilot','-ascii')
% 
% % obtain weighting matrix for ABC summaries
% summ_weights = diag(mad(summaries_pilot,0).^2);
% 
% return

%summobs = gk_abc_summaries(yobs);
summ_weights = eye(9*numdatasets);  % dimension is (number of quantiles) * numdatasets

%::::::: ABC INFERENCE :::::::::::::::::::::::::::::::::::::::::::::::

rng(200)

% ABC settings
numparticles = 10000;  % notice the first SMC-ABC iteration will take several minutes to complete with these many particles
ABCthreshold = 50;  % starting threshold
alpha = 25;  % this is a percentile (\psi in the paper), determining the ABC threshold at each iteration. 
problem = 'gk'; % a string specifying the name of a case study
ABCthreshold_final = 0.62; % stop as soon as the threshold becomes smaller than this value 
                          % WARNING: you may want to enlarge ABCthreshold_final. Most samplers will struggle to reach this value.
numattempt = 1;           % number of independent repetitions of the experiment, on the same data 

tic
for attempt=1:numattempt
   ABCdraws = abc_sequential_varying_thresholds(problem,yobs,covariates,parmask,parbase,ABCthreshold,ABCthreshold_final,summ_weights,numparticles,alpha,sampling,attempt,nu,type_copula,type_marginals,nu_marg,folder);
end
eval = toc