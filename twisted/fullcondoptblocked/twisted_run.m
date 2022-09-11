% this run file is a bit different from the other ones in that it runs a
% specific sampler that is not considered among the "general" samplers
% file (abc_sequential_varying_thresholds). This one calls instead abc_sequential_fullcondopt_blockingtheta1theta2
% and is the one denoted 'fullcondoptblocked' in the paper. Unlike
% full_cond_opt (which proposes 1 parameter coordinate at the time),
% fullcondoptblocked proposes jointly the first two coordinates of theta (conditionally on the
% others and summaries) since these are highly correlated, while the 
% remaining three coordinates are sampled individually. 

rng(100)

covariates = [];  % no covariates needed for this model

% set data
yobs = zeros(5,1);  %!!! MUST BE A VERTICAL VECTOR !!!!
dimensionality = length(yobs);
yobs(1) = 10;

%::: ABC PILOT to collect simulated summaries and compute a matrix that
%weights them
size_pilot = 5000;
summaries_pilot = zeros(size_pilot,dimensionality);
abc_distances = zeros(size_pilot,1);
summobs = twisted_abc_summaries(yobs);
for ii=1:size_pilot
    prior_draw = twisted_prior([],1);
    x = twisted_model(prior_draw,covariates,1);
    sim_summaries = twisted_abc_summaries(x);
    summaries_pilot(ii,:) = sim_summaries;
end

save('summaries_pilot','summaries_pilot')
save('summaries_pilot.dat','summaries_pilot','-ascii')

% obtain weighting matrix for ABC summaries
summ_weights = diag(mad(summaries_pilot,0).^2);


rng(200)
                  % t1   t2 
parbase        = [ zeros(1,dimensionality)];    % just some "baseline parameter values", e.g. prior means
parmask        = [ ones(1,dimensionality)];  % 1= this is a parameter to infer; 0 = this is a known constant
% that is sum(parmask) is the number of parameters to infer

% ABC settings
numparticles = 1000;
ABCthreshold = 50;
alpha = 1;
sampling = 'full_cond_opt';
problem = 'twisted';
ABChreshold_final = 0.25; % stop as soon as the threshold becomes smaller than this value
                          % WARNING: you may want to enlarge ABCthreshold_final. Most samplers will struggle to reach this value.
numattempt = 1;

tic
for attempt=1:numattempt
    ABCdraws = abc_sequential_fullcondopt_blockingtheta1theta2(problem,yobs,covariates,parmask,parbase,ABCthreshold,ABChreshold_final,summ_weights,numparticles,alpha,sampling,attempt);
end
eval = toc


