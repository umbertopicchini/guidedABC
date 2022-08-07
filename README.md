# guidedABC
Supporting code for Picchini &amp; Tamborrino (2022) "Guided sequential ABC schemes for intractable Bayesian models" [arxiv:2206.12235](https://arxiv.org/abs/2206.12235)

Each subfolder implements a different case-study. Go into the desired folder then locate the "run" file and launch it. The run files are heavily commented to provide guidance. Notice that output files are automatically stored into an appropriate subfolder, which however *must first be created* by the user as a subfolder of the desired case-study, and should be named exactly as the desired proposal sampler, ie one among 'blocked', 'blockedopt', 'hybrid', 'fullcond', 'full_cond_opt', 'olcm', 'standard'.  

Example: if you wish to run the 'blocked' sampler for the twomoons case study, first create a folder twomoons\blocked, then run twomoons_run.
