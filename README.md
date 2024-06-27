# guidedABC
MATLAB code for Picchini &amp; Tamborrino (2022) "Guided sequential ABC schemes for intractable Bayesian models", forthcoming in Bayesian Analyis (https://doi.org/10.1214/24-BA1451), see also [arxiv:2206.12235](https://arxiv.org/abs/2206.12235)

Each subfolder implements a different case-study, see the paper for details. Go into the desired folder then locate the "run" file and launch it (but first create folders for the output files, see below). The run files are heavily commented to provide guidance. Notice that output files are automatically stored into an appropriate subfolder, which however *must first be created* by the user as a subfolder of the desired case-study, and should be named exactly as the desired proposal sampler, ie one among 'blocked', 'blockedopt', 'hybrid', 'fullcond', 'full_cond_opt', 'olcm', 'standard'. 

Example: if you wish to run the 'blocked' sampler for the 'twomoons' case study, first create a folder twomoons/blocked, then run twomoons_run.

Output files: the subfolder with inference results will contain many files upon completion of a run. Here is a description:

- posterior draws: a typical file name ABCdraws_stageX_attemptY_copulaZ_margA_nuB_numargC_partN
  This is a matrix of d x N values containing N parameters draws (particles) for each of the d parameters to infer. The X in 'stageX' is the iteration number the draws have been sampled at. The Y in 'attemptY' is the ID of the inference run (since it is possible to run several inference runs on the same dataset, one after the other). Z, A, B and C in 'copulaZ_margA_nuB_numargC' refer to the type of copula (Z=0 means no copula method is used. See the comments in the run file for info about Z, A, B and C). N in 'partN' denotes the desired number of accepted parameter draws, also named "particles".
  
The following files have a similar structure as the posterior draws files. We do not repeat the similar features and only focus on the differences:

- ess_stageX_attemptY_copulaZ_margA_nuB_numargC_partN: the value of the effective sample size (ESS) at the given iteration, at the given attempt etc.
- evaltime_stageX_attemptY_copulaZ_margA_nuB_numargC_partN: the number of seconds required to accept N draws at the given iteration, at the given attempt etc.
- numproposals_stageX_attemptY_copulaZ_margA_nuB_numargC_partN: the number of particles proposed (ie both the accepted and rejected), at the given iteration, attempt etc.
- weights_stageX_attemptY_copulaZ_margA_nuB_numargC_partN: the vector of N normalised importance weights associated to each of the N accepted particles at the given iteration, at the given attempt etc.
- ABCthreshold_stageX_attemptY_copulaZ_margA_nuB_numargC_partN: the value of the ABC threshold used at the given iteration, attempt etc.

