% example from Price, L. F., Drovandi, C. C., Lee, A., & Nott, D. J. (2018). Bayesian synthetic likelihood. Journal of Computational and Graphical Statistics, 27(1), 1-11.
% code adapted from the folder "Scratch Assay" found at https://github.com/cdrovandi/Bayesian-Synthetic-Likelihood

rng(100)

% If you do not use Windows, you should first compile with the following
% mex -compatibleArrayDims simulate_mex.c mt19937ar.c

% a Windows64 compiled version is included. But if you work with a
% different platform you will have to run the above mex instruction.

% Now load ground truth parameters Pm = 0.36 and Pp=1e-3 and covariates.
% We use the same umber of summaries as in Price et al.
% the same initial number of cells Yinit, tau=1/24 and
% ground-truth parameters
load('all_locations_simulated.mat');
covariates.Yinit = Yinit;
covariates.Yobs = S;
covariates.obstimes = (1/12):(1/12):12;
data.Yinit = Yinit;
data.Y = S;

% Possible sampling: blocked, blockedopt, hybrid, fullcond, full_cond_opt, olcm, standard.
sampling = "olcm";
% as an example here we have set the guided "hybrid" method. 
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

folder='olcm'
% "folder" is the name of the folder where you want your results to be stored in. It MUST be created before running this script or else it will result in an error.
% The whole path to access the folder can be given as input, allowing you to run the script from any other folder.
if ~exist(folder, 'dir') || ~strcmp(folder,sampling)
    error('A folder to store your results must be created first and must be named exactly as the sampler you wish to use.')
end


% %::: ABC PILOT to collect simulated summaries and compute a matrix that weights them
% size_pilot = 10000;
% summaries_pilot = zeros(size_pilot,1+length(covariates.obstimes));
% abc_distances = zeros(size_pilot,1);
% summobs = cell_abc_summaries(data);
% for ii=1:size_pilot
%     prior_draw = cell_prior([],1);
%     Pm = prior_draw(1);
%     Pp = prior_draw(2);
%     data_sim = cell_model([Pm Pp],covariates,[]);
%     sim_summaries = cell_abc_summaries(data_sim);
%     summaries_pilot(ii,:) = sim_summaries;
% end
% save('summaries_pilot.dat','summaries_pilot','-ascii')
% 
% return

load('summaries_pilot.dat')
% obtain weighting matrix for ABC summaries
summ_weights = diag(mad(summaries_pilot,0).^2);

%::::::: ABC INFERENCE :::::::::::::::::::::::::::::::::::::::::::::::


problem = 'cell';
numparticles = 1000;
alpha = 1;         % a percentile value in (0,100) to compute percentiles of ABC distances
ABCthreshold = 	788;  
ABCthreshold_final = 1;  
numattempt = 1;

                  % Pm   Pp   
parbase        = [0.5   0.5 ];    % just some unused "baseline parameter values", e.g. prior means
parmask        = [ 1    1 ];  % 1= this is a parameter to infer; 0 = this is a known constant


tic
for attempt=1:numattempt
   ABCdraws = abc_sequential_varying_thresholds(problem,data,covariates,parmask,parbase,ABCthreshold,ABCthreshold_final,summ_weights,numparticles,alpha,sampling,attempt,nu,type_copula,type_marginals,nu_marg,folder);
end
eval = toc
%tic
%for attempt=1:numattempts
%  ABCdraws = abc_sequential_unifkernel_unnormalised(problem,data,covariates,parmask,parbase,ABCthreshold,ABCthreshold_final,summ_weights,numparticles,alpha,sampling,regressadj,attempt);
%end
%eval = toc



