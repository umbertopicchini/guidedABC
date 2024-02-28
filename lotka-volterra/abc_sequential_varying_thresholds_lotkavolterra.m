function ABCdraws = abc_sequential_varying_thresholds(problem,data,covariates,parmask,parbase,ABCthreshold,ABCthreshold_final,summ_weights,numparticles,alpha,sampling,attempt,nu,type_copula,type_margin,nu_marg,folder)

% if any(isnan(data))
%     error('Supplied data seems to contain NaNs')
% end

summobs = feval([problem, '_abc_summaries'],data); % vector of summaries for the observed data


nfreepar = sum(parmask);
ABCdraws = zeros(nfreepar,numparticles);
nsummary = length(summobs); % the number of summary statistics
distance_all = [];
simsumm_all =  [];
simsumm_accepted = zeros(nsummary,numparticles);
distance_accepted = zeros(1,numparticles);  % this is only for the olcm proposal 

% initialization: here t=1
success = 0;
t = 1;
numproposals = 0;
tStart = tic;

while success < numparticles
    success
    numproposals = numproposals +1;
          theta = feval([problem, '_prior'],[],1);
          bigtheta = param_unmask(theta,parmask,parbase);
          simdata = feval([problem, '_model'],bigtheta,covariates,1);
          simsumm = feval([problem, '_abc_summaries'],simdata);
          xc = bsxfun(@minus,simsumm',summobs');
          distance = sqrt(sum((xc / summ_weights) .* xc, 2));
       if distance < ABCthreshold
           success = success+1;
           simsumm_accepted(:,success) =  simsumm;
           simsumm_all =  [simsumm_all,simsumm];
           ABCdraws(:,success) = theta;
           distance_all = [distance_all,distance'];
           switch sampling
               case 'olcm' % do this to compute the covariance for this specific sampler (only for accepted proposals)
                   distance_accepted(success) = distance;  
               case 'full_cond_opt'
                   distance_accepted(success) = distance;  
               case 'blockedopt'
                   distance_accepted(success) = distance;  
           end
        else
           distance_all = [distance_all,distance'];
           simsumm_all =  [simsumm_all,simsumm];
        %   simsumm_all(:,attempts+1:attempts+numsim) = feval([problem, '_abc_summaries'],simdata);
        end
end


% remove possible NaNs from simulated summaries
idnan = any(isnan(simsumm_all'),2); % find units with NaNs
summaries_fixed = simsumm_all(:,~idnan);  % remove those units
% remove nasty outliers 
%simsumm_all = rmoutliers(summaries_fixed','percentiles',[5,95]);
%simsumm_all = simsumm_all';
%size(simsumm_all)
simsumm_all = summaries_fixed;
% obtain weighting matrix for ABC summaries
summ_weights = diag(mad(simsumm_all',0).^2);
save(sprintf('%s/summ_weights.dat',folder),'summ_weights','-ascii')

eval_time = toc(tStart); 
save(sprintf('%s/evaltime_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%1.1f_part%d.dat',folder,t,attempt,type_copula,type_margin,nu,nu_marg,numparticles),'eval_time','-ascii')  
save(sprintf('%s/ABCthreshold_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%1.1f_part%d.dat',folder,t,attempt,type_copula,type_margin,nu,nu_marg,numparticles),'ABCthreshold','-ascii')
save(sprintf('%s/ABCdraws_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%1.1f_part%d.dat',folder,t,attempt,type_copula,type_margin,nu,nu_marg,numparticles),'ABCdraws','-ascii')
save(sprintf('%s/numproposals_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%1.1f_part%d.dat',folder,t,attempt,type_copula,type_margin,nu,nu_marg,numparticles),'numproposals','-ascii')
weights = ones(1,numparticles);


normweights = weights./sum(weights);
ess = 1/sum(normweights.^2)  % the Effective Sample Size
save(sprintf('%s/ess_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%1.1f_part%d.dat',folder,t,attempt,type_copula,type_margin,nu,nu_marg,numparticles),'ess','-ascii')
save(sprintf('%s/weights_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%1.1f_part%d.txt',folder,t,attempt,type_copula,type_margin,nu,nu_marg,numparticles),'normweights','-ascii')
all_accepted_thetasimsum = [ABCdraws;simsumm_accepted];
%all_cov_thetasimsum = cov(all_accepted_thetasimsum');

all_cov_thetasimsum = all_accepted_thetasimsum' - repmat(normweights * all_accepted_thetasimsum', numparticles, 1);   % subtract weighted mean                                                
all_cov_thetasimsum = all_cov_thetasimsum' * (all_cov_thetasimsum .* repmat(normweights', 1, nfreepar+nsummary));                                                  
all_cov_thetasimsum = all_cov_thetasimsum ./ (1-sum(normweights.^2)); % the weighted covariance matrix
all_cov_thetasimsum = 0.5 * (all_cov_thetasimsum + all_cov_thetasimsum');   % ensure symmetry

if t==1 && strcmp(sampling,'hybrid') % if we have an hybrid strategy then for t==2 (and only for t==2) propose from the blocked sampler
        sampling_temp = 'hybrid'; 
        sampling = 'blocked';
else
    sampling_temp = [];
end

switch sampling
    case 'blocked'
        mean_theta = ABCdraws * normweights';  % weighted mean across all accepted parameters
        mean_simsum = simsumm_accepted * normweights'; % weighted mean across all accepted summaries
        cov_theta = all_cov_thetasimsum(1:nfreepar,1:nfreepar);
        cov_simsum = all_cov_thetasimsum(nfreepar+1:end,nfreepar+1:end);
        cov_thetasimsum = all_cov_thetasimsum(1:nfreepar,nfreepar+1:end);
        cov_simsumtheta = cov_thetasimsum';
        proposal_mean = mean_theta + cov_thetasimsum * (cov_simsum \(summobs-mean_simsum));
        proposal_mean =  proposal_mean';  % must be a row vector when passed to mvnrnd()
        proposal_cov = cov_theta - cov_thetasimsum * (cov_simsum \ cov_simsumtheta);
        proposal_cov = (proposal_cov + proposal_cov.') / 2; % standard trick to ensure simmetry
        [~, notposdef] = cholcov(proposal_cov,0);
        if isnan(notposdef) || notposdef>0
            [L, DMC, P] = modchol_ldlt(proposal_cov); 
            proposal_cov = P'*L*DMC*L'*P;
        end 
    case 'blockedopt'
        mean_theta = ABCdraws * normweights';  % weighted mean across all accepted parameters
        mean_simsum = simsumm_accepted * normweights'; % weighted mean across all accepted summaries
        cov_simsum = all_cov_thetasimsum(nfreepar+1:end,nfreepar+1:end);
        cov_thetasimsum = all_cov_thetasimsum(1:nfreepar,nfreepar+1:end);
        cov_simsum = (cov_simsum+cov_simsum')/2;
        [~, notposdef] = cholcov(cov_simsum,0);
        if isnan(notposdef) || notposdef>0
          [L, DMC, P] = modchol_ldlt(cov_simsum); 
           cov_simsum = P'*L*DMC*L'*P;
        end
        proposal_mean = mean_theta + cov_thetasimsum * (cov_simsum \(summobs-mean_simsum));
        cov_blockedopt = zeros(nfreepar,nfreepar);
        ABCthreshold_temp = prctile(distance_all,alpha);
        ABCthreshold_old = ABCthreshold;
        ABCthreshold_new = min(ABCthreshold_temp,ABCthreshold_old); % ensure non-increasing thresholds
        if ABCthreshold_new == ABCthreshold_old
           ABCthreshold_new = 0.95*ABCthreshold_new;
        end
        id_blockedopt = find(distance_accepted < ABCthreshold_new);
        N0 = length(id_blockedopt);  % number of particles that satify the newest threshold
        if N0==0  % in this case there are no particles satifying the criterion above.
                  % we must exit
            fprintf('\n There are no particles satisfying the blockedopt criterion. We must leave...\n')
            return
        end
        weights_blockedopt = weights(id_blockedopt);
        normweights_blockedopt = weights_blockedopt/sum(weights_blockedopt);
        for jj=1:N0
              cov_blockedopt = cov_blockedopt + normweights_blockedopt(jj)*(ABCdraws(:,id_blockedopt(jj))-proposal_mean)*(ABCdraws(:,id_blockedopt(jj))-proposal_mean)';
        end
        cov_blockedopt = (cov_blockedopt+cov_blockedopt')/2;
        [~, notposdef] = cholcov(cov_blockedopt,0);
        if isnan(notposdef) || notposdef>0
          [L, DMC, P] = modchol_ldlt(cov_blockedopt); 
           cov_blockedopt = P'*L*DMC*L'*P;
        end 
        proposal_mean =  proposal_mean';  % must be a row vector when passed to mvnrnd()            
    case 'full_cond'
       mean_theta = ABCdraws * normweights';  % weighted mean across all accepted parameters
       mean_simsum = simsumm_accepted * normweights'; % weighted mean across all accepted summaries
       proposal_mean = zeros(nfreepar,numparticles);
       proposal_var  = zeros(nfreepar,numparticles);
       for ii = 1:numparticles
          for par=1:nfreepar
             % notice the ~=par operator means "select all entries in the array except those corresponding to index par"
             proposal_mean(par,ii) = mean_theta(par) + all_cov_thetasimsum(par,1:end~=par) * (all_cov_thetasimsum(1:end~=par,1:end~=par) \ ([ABCdraws(1:end~=par,ii);summobs] - [mean_theta(1:end~=par);mean_simsum]));
             proposal_var(par,ii) = all_cov_thetasimsum(par,par) - all_cov_thetasimsum(par,1:end~=par) * (all_cov_thetasimsum(1:end~=par,1:end~=par) \ all_cov_thetasimsum(1:end~=par,par));
          end
       end
    case 'standard'
        % compute the weighted covariance matrix C
        C = ABCdraws' - repmat(normweights * ABCdraws', numparticles, 1);   % subtract weighted mean                                                  % Remove mean (which is, also, weighted)
        C = C' * (C .* repmat(normweights', 1, nfreepar));                                                  % Weighted Covariance Matrix
        C = C ./ (1-sum(normweights.^2)); % the weighted covariance matrix
        C = 0.5 * (C + C');   % ensure symmetry
        Sigma = 2*C;
    case 'olcm'   % considers the local optimal covariance denoted OLCM in Filippi, S., Barnes, C. P., Cornebise, J., & Stumpf, M. P. (2013). On optimality of kernels for approximate Bayesian computation using sequential Monte Carlo. Statistical applications in genetics and molecular biology, 12(1), 87-107.
        cov_olcm = zeros(nfreepar,nfreepar);
        ABCthreshold_temp = prctile(distance_all,alpha);
        ABCthreshold_old = ABCthreshold;
        ABCthreshold_new = min(ABCthreshold_temp,ABCthreshold_old); % ensure non-increasing thresholds
        if ABCthreshold_new == ABCthreshold_old
           ABCthreshold_new = 0.95*ABCthreshold_new;
        end
        id_olcm = find(distance_accepted < ABCthreshold_new);
        N0 = length(id_olcm);  % number of particles that satify the newest threshold
        if N0==0  % in this case there are no particles satifying the criterion above.
                  % we must exit
            fprintf('\n There are no particles satisfying the OLCM criterion. We must leave...\n')
            return
        end
        weights_olcm = weights(id_olcm);
        normweights_olcm = weights_olcm/sum(weights_olcm);
        cov_olcm_all = [];
        for ii=1:numparticles
           for jj=1:N0
              cov_olcm = cov_olcm + normweights_olcm(jj)*(ABCdraws(:,id_olcm(jj))-ABCdraws(:,ii))*(ABCdraws(:,id_olcm(jj))-ABCdraws(:,ii))';
           end
          cov_olcm = (cov_olcm+cov_olcm')/2;
          cov_olcm_all = [cov_olcm_all,cov_olcm];
        end
    case 'full_cond_opt'
        mean_theta = ABCdraws * normweights';  % weighted mean across all accepted parameters
        mean_simsum = simsumm_accepted * normweights'; % weighted mean across all accepted summaries
        proposal_mean = zeros(nfreepar,numparticles);
        var_opt_all = zeros(nfreepar,numparticles);
        ABCthreshold_temp = prctile(distance_all,alpha);
        ABCthreshold_old = ABCthreshold;
        ABCthreshold_new = min(ABCthreshold_temp,ABCthreshold_old) % ensure non-increasing thresholds
        if ABCthreshold_new == ABCthreshold_old
           ABCthreshold_new = 0.95*ABCthreshold_new;
        end
        id_opt = find(distance_accepted < ABCthreshold_new);
        N0 = length(id_opt)  % number of particles that satify the newest threshold
        if N0==0  % in this case there are no particles satifying the criterion above.
                  % we must exit
            fprintf('\n There are no particles satisfying the OLCM criterion. We must leave...\n')
            return
        end
        weights_opt = weights(id_opt);
        normweights_opt = weights_opt/sum(weights_opt);
        all_cov_thetasimsum = (all_cov_thetasimsum+all_cov_thetasimsum')/2;
        [~, notposdef] = cholcov(all_cov_thetasimsum,0);
        if isnan(notposdef) || notposdef>0
          [L, DMC, P] = modchol_ldlt(all_cov_thetasimsum); 
           all_cov_thetasimsum = P'*L*DMC*L'*P;
        end
        for ii = 1:numparticles
           for par=1:nfreepar
              % notice the ~=par operator means "select all entries in the array except those corresponding to index par"
              proposal_mean(par,ii) = mean_theta(par) + all_cov_thetasimsum(par,1:end~=par) * (all_cov_thetasimsum(1:end~=par,1:end~=par) \ ([ABCdraws(1:end~=par,ii);summobs] - [mean_theta(1:end~=par);mean_simsum]));
              for jj=1:N0
                  var_opt_all(par,ii) = var_opt_all(par,ii) + normweights_opt(jj)*(ABCdraws(par,id_opt(jj))-proposal_mean(par,ii))^2;
              end
           end
        end
    otherwise
        error('SAMPLING must be one of the following options: "full_cond", "full_cond_opt" , "blocked", "blockedopt", "olcm", "standard".')
end

%ABCthreshold = min(prctile(distance_all,alpha),ABCthreshold);
%Sigma = 2*weightedcov(ABCdraws', normweights);  ´



while ABCthreshold > ABCthreshold_final  % keep doing as long as acceptance rate requirement is satsified
    t = t+1;
    
    if t==2 && strcmp(sampling_temp,'hybrid') % if we have an hybrid strategy then for t==2 (and only for t==2) propose from the blocked sampler
        sampling = 'blocked';
    elseif t>2 && strcmp(sampling_temp,'hybrid')
        sampling = 'blockedopt';
    end
    
    ABCthreshold_temp = prctile(distance_all,alpha);
    ABCthreshold_old = ABCthreshold;
    ABCthreshold = min(ABCthreshold_temp,ABCthreshold_old) % ensure non-increasing thresholds
    if ABCthreshold == ABCthreshold_old
        ABCthreshold = 0.95*ABCthreshold;
        fprintf('Forced decrease of ABC threshold to: %d',ABCthreshold)
    end


%    ABCthreshold_old = ABCthreshold;
%    ABCthreshold_temp = prctile(distance_all,alpha);
%    ABCthreshold = min(ABCthreshold_temp,ABCthreshold_old) % ensure non-incresing thresholds
    
   weights_old = weights;
    
    switch sampling
        case 'full_cond'
            proposal_mean_old = proposal_mean;
            proposal_var_old = proposal_var;
        case 'standard'
            ABCdraws_old = ABCdraws;
            Sigma_old = Sigma;
        case 'olcm'
            ABCdraws_old = ABCdraws;
            % now find the particles that would have been accepted (if any) under the NEW threshold
            id_olcm = find(distance_accepted < ABCthreshold);
            N0 = length(id_olcm)  % number of particles that satify the newest threshold
            if N0==0  % in this case there are no particles satifying the criterion above.
                      % we must exit
                fprintf('\n There are no particles satisfying the OLCM criterion. We must leave...\n')
                return
            end
            distance_accepted = zeros(1,numparticles);
            weights_olcm = weights_old(id_olcm);
            normweights_olcm = weights_olcm/sum(weights_olcm);   
        case 'full_cond_opt'
            ABCdraws_old = ABCdraws;
            proposal_mean_old = proposal_mean;
            % now find the particles that would have been accepted (if any) under the NEW threshold
            id_opt = find(distance_accepted < ABCthreshold);
            N0 = length(id_opt)  % number of particles that satify the newest threshold
            if N0==0  % in this case there are no particles satifying the criterion above.
                      % we must exit
                fprintf('\n N0=0. We must leave...\n')
                return
            end
            distance_accepted = zeros(1,numparticles);
            weights_opt = weights_old(id_opt);
            normweights_opt = weights_opt/sum(weights_opt);
    end
    
    distance_all = [];
    simsumm_all = [];
    simsumm_accepted = zeros(nsummary,numparticles);
    
    success = 0;
    numproposals = 0;
   
    tStart = tic;
    while success < numparticles
           success
           numproposals = numproposals +1;
          
           switch sampling
               case 'full_cond'
                   % sample a particle from the weighted set of particles
                   % from the previous generation
                   index = stratresample(normweights,1);
                   for par=1:nfreepar
                      theta(par) = proposal_mean_old(par,index) + sqrt(proposal_var_old(par,index)) * randn;
                   end
               case 'blocked'
                     %theta = mvnrnd(proposal_mean,proposal_cov);
                    [AA,theta,pd,alphaF]=sample_theta(proposal_mean,proposal_cov,nu,type_copula,type_margin,nu_marg);
               case 'blockedopt'
%                   theta = mvnrnd(proposal_mean,cov_blockedopt);
                    [AA,theta,pd,alphaF]=sample_theta(proposal_mean,cov_blockedopt,nu,type_copula,type_margin,nu_marg);
               case 'standard'
                   index = stratresample(normweights,1);
                   theta = mvnrnd(ABCdraws_old(:,index),Sigma_old);
               case 'olcm'
                   index = stratresample(normweights,1);
                   % id_olcm contains all the indeces of the particles from the previous generation that satify the
                   % current threshold. There are N0 of those.
                   cov_olcm = zeros(nfreepar,nfreepar);
                   for jj=1:N0
                     cov_olcm = cov_olcm + normweights_olcm(jj)*(ABCdraws_old(:,id_olcm(jj))-ABCdraws_old(:,index))*(ABCdraws_old(:,id_olcm(jj))-ABCdraws_old(:,index))';
                   end
                   cov_olcm = (cov_olcm+cov_olcm')/2;
                    [~, notposdef] = cholcov(cov_olcm,0);
                    if isnan(notposdef) || notposdef>0
                        [L, DMC, P] = modchol_ldlt(cov_olcm); 
                      %  perturbation = P'*L*DMC*L'*P - proposal_cov;
                         cov_olcm = P'*L*DMC*L'*P;
                      %     cov_olcm = nearestSPD(cov_olcm);
                    end
                   % the above covariance is not "global" but is instead specific for the sampled particle
                   theta = mvnrnd(ABCdraws_old(:,index),cov_olcm);
               case 'full_cond_opt'
                   index = stratresample(normweights,1);
                   % id_opt contains all the indeces of the particles from the previous generation that satify the
                   % current threshold. There are N0 of those.
                   var_opt = zeros(nfreepar,1);
                   for par = 1:nfreepar
                      for jj=1:N0
                         var_opt(par) = var_opt(par) + normweights_opt(jj)*(ABCdraws_old(par,id_opt(jj))-proposal_mean_old(par,index))^2;
                      end
                   end
                   % the above variance is not "global" but is instead
                   % specific for the sampled particle since it depends on proposal_mean_old(par,index)
                   for par=1:nfreepar
                      theta(par) = proposal_mean_old(par,index) + sqrt(var_opt(par)) * randn;
                   end
           end           

           prior = feval([problem, '_prior'],theta,0);
           bigtheta = param_unmask(theta,parmask,parbase);
           simdata = feval([problem, '_model'],bigtheta,covariates,1);
           simsumm = feval([problem, '_abc_summaries'],simdata);
           xc = bsxfun(@minus,simsumm',summobs');
           distance = sqrt(sum((xc / summ_weights) .* xc, 2));
%            loglike = -nsummary*log(ABCthreshold) +logsumexp(-distance.^2/(2*ABCthreshold^2));
%            logKernel_0 = -nsummary*log(ABCthreshold);
%            if log(rand) < loglike - logKernel_0
           if distance < ABCthreshold
              success = success+1;
              ABCdraws(:,success) = theta;
              simsumm_accepted(:,success) =  simsumm;
              distance_all = [distance_all,distance'];
              simsumm_all =  [simsumm_all,simsumm];
              switch sampling
               case 'full_cond'
                     % compute the weights denominator
                     denominator = 0;
                     for ii=1:numparticles
                            denominator_prod =1;
                            for par=1:nfreepar
                              denominator_prod = denominator_prod *normpdf(theta(par),proposal_mean_old(par,ii),sqrt(proposal_var_old(par,ii)));
                            end
                            denominator = denominator + weights_old(ii)*denominator_prod;
                     end
                     weights(success) = prior / denominator ;
               case 'full_cond_opt'
                     distance_accepted(success) = distance;
                     % compute the weights denominator
                     denominator = 0;
                     for ii=1:numparticles
                            denominator_prod =1;
                            for par=1:nfreepar
                              denominator_prod = denominator_prod *normpdf(theta(par),proposal_mean_old(par,ii),sqrt(var_opt_all(par,ii)));
                            end
                            denominator = denominator + weights_old(ii)*denominator_prod;
                     end
                     weights(success) = prior / denominator ;
               case 'blocked'
                      dens=density_sample(AA,theta,pd,alphaF,proposal_mean,proposal_cov,nu,type_copula,type_margin);
                     weights(success) = prior /dens;%prior/mvnpdf(theta,proposal_mean,proposal_cov); % weights of accepted particles     
               case 'blockedopt'
                     distance_accepted(success) = distance; 
                     dens=density_sample(AA,theta,pd,alphaF,proposal_mean,cov_blockedopt,nu,type_copula,type_margin);
                     weights(success) = prior /dens;%prior/mvnpdf(theta,proposal_mean,cov_blockedopt); % weights of accepted particles     
               case 'standard'
                     denominator =0;
                     for ii=1:numparticles
                         denominator = denominator + weights_old(ii)*mvnpdf(theta',ABCdraws_old(:,ii),Sigma_old);
                     end
                     weights(success) = prior / denominator ;   
               case 'olcm'
                     distance_accepted(success) = distance; 
                     denominator =0;
                   %  cov_olcm = (cov_olcm+cov_olcm')/2
                     for ii=1:numparticles
                         denominator = denominator + weights_old(ii)*mvnpdf(theta',ABCdraws_old(:,ii),cov_olcm_all(:,(ii-1)*nfreepar+1:ii*nfreepar));
                     end
                     weights(success) = prior / denominator ;
              end
           else  % proposal is rejected, but we still collect simulated summaries and distances for normalization purposes and threshold update
              distance_all = [distance_all,distance'];
              simsumm_all =  [simsumm_all,simsumm];
           end
    end


    normweights = weights./sum(weights);
    ess = 1/sum(normweights.^2)  % the Effective Sample Size
    eval_time = toc(tStart); 

    save(sprintf('%s/evaltime_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%1.1f_part%d.dat',folder,t,attempt,type_copula,type_margin,nu,nu_marg,numparticles),'eval_time','-ascii')  
    save(sprintf('%s/ABCthreshold_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%1.1f_part%d.dat',folder,t,attempt,type_copula,type_margin,nu,nu_marg,numparticles),'ABCthreshold','-ascii')
    save(sprintf('%s/ABCdraws_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%1.1f_part%d.dat',folder,t,attempt,type_copula,type_margin,nu,nu_marg,numparticles),'ABCdraws','-ascii')
    save(sprintf('%s/numproposals_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%1.1f_part%d.dat',folder,t,attempt,type_copula,type_margin,nu,nu_marg,numparticles),'numproposals','-ascii')
    save(sprintf('%s/ess_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%1.1f_part%d.dat',folder,t,attempt,type_copula,type_margin,nu,nu_marg,numparticles),'ess','-ascii')
    save(sprintf('%s/weights_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%1.1f_part%d.txt',folder,t,attempt,type_copula,type_margin,nu,nu_marg,numparticles),'normweights','-ascii')
    %  summ_weights = diag(mad(simsumm_all',0).^2);
  
    all_accepted_thetasimsum = [ABCdraws;simsumm_accepted];
    all_cov_thetasimsum = all_accepted_thetasimsum' - repmat(normweights * all_accepted_thetasimsum', numparticles, 1);   % subtract weighted mean                                     
    all_cov_thetasimsum = all_cov_thetasimsum' * (all_cov_thetasimsum .* repmat(normweights', 1, nfreepar + nsummary));                                                  
    all_cov_thetasimsum = all_cov_thetasimsum ./ (1-sum(normweights.^2)); % the weighted covariance matrix
    all_cov_thetasimsum = 0.5 * (all_cov_thetasimsum + all_cov_thetasimsum');   % ensure symmetry

    
    switch sampling
    case 'blocked'
       if t==2 && strcmp(sampling_temp,'hybrid')  % here we are in the hybrid case so we construct the blockedopt sampler for next iteration t>2
           %::: PREPARING FOR THE HYBRID STRATEGY ::::::::::::::::::::::
           mean_theta = ABCdraws * normweights';  % weighted mean across all accepted parameters
           mean_simsum = simsumm_accepted * normweights'; % weighted mean across all accepted summaries
           cov_simsum = all_cov_thetasimsum(nfreepar+1:end,nfreepar+1:end);
           cov_thetasimsum = all_cov_thetasimsum(1:nfreepar,nfreepar+1:end);
           cov_simsum = (cov_simsum+cov_simsum')/2;
           [~, notposdef] = cholcov(cov_simsum,0);
           if isnan(notposdef) || notposdef>0
             [L, DMC, P] = modchol_ldlt(cov_simsum); 
             cov_simsum = P'*L*DMC*L'*P;
           end
           proposal_mean = mean_theta + cov_thetasimsum * (cov_simsum \(summobs-mean_simsum));
           cov_blockedopt = zeros(nfreepar,nfreepar);
           ABCthreshold_temp = prctile(distance_all,alpha);
           ABCthreshold_old = ABCthreshold;
           ABCthreshold_new = min(ABCthreshold_temp,ABCthreshold_old); % ensure non-increasing thresholds
           if ABCthreshold_new == ABCthreshold_old
              ABCthreshold_new = 0.95*ABCthreshold_new;
           end
           id_blockedopt = find(distance_accepted < ABCthreshold_new);
           N0 = length(id_blockedopt);  % number of particles that satify the newest threshold
           if N0==0  % in this case there are no particles satifying the criterion above.
                  % we must exit
              fprintf('\n There are no particles satisfying the blockedopt criterion. We must leave...\n')
              return
           end
           weights_blockedopt = weights(id_blockedopt);
           normweights_blockedopt = weights_blockedopt/sum(weights_blockedopt);
           for jj=1:N0
               cov_blockedopt = cov_blockedopt + normweights_blockedopt(jj)*(ABCdraws(:,id_blockedopt(jj))-proposal_mean)*(ABCdraws(:,id_blockedopt(jj))-proposal_mean)';
           end
           cov_blockedopt = (cov_blockedopt+cov_blockedopt')/2;
           [~, notposdef] = cholcov(cov_blockedopt,0);
           if isnan(notposdef) || notposdef>0
             [L, DMC, P] = modchol_ldlt(cov_blockedopt); 
             cov_blockedopt = P'*L*DMC*L'*P;
           end 
           proposal_mean =  proposal_mean';  % must be a row vector when passed to mvnrnd()
       else
            %:::: this is standard BLOCKED sampler :::::::::::::::::::::::
            mean_theta = ABCdraws * normweights';  % weighted mean across all accepted parameters
            mean_simsum = simsumm_accepted * normweights'; % weighted mean across all accepted summaries
            cov_theta = all_cov_thetasimsum(1:nfreepar,1:nfreepar);
            cov_simsum = all_cov_thetasimsum(nfreepar+1:end,nfreepar+1:end);
            cov_thetasimsum = all_cov_thetasimsum(1:nfreepar,nfreepar+1:end);
            cov_simsumtheta = cov_thetasimsum';
            proposal_mean = mean_theta + cov_thetasimsum * (cov_simsum \(summobs-mean_simsum));
            proposal_mean =  proposal_mean';  % must be a row vector when passed to mvnrnd()
            proposal_cov = cov_theta - cov_thetasimsum * (cov_simsum \ cov_simsumtheta);
            proposal_cov = (proposal_cov + proposal_cov.') / 2; % standard trick to ensure simmetry
            [~, notposdef] = cholcov(proposal_cov,0);
            if isnan(notposdef) || notposdef>0
               [L, DMC, P] = modchol_ldlt(proposal_cov); 
               proposal_cov = P'*L*DMC*L'*P;
            end  
       end
    case 'blockedopt' 
        mean_theta = ABCdraws * normweights';  % weighted mean across all accepted parameters
        mean_simsum = simsumm_accepted * normweights'; % weighted mean across all accepted summaries
        cov_simsum = all_cov_thetasimsum(nfreepar+1:end,nfreepar+1:end);
        cov_thetasimsum = all_cov_thetasimsum(1:nfreepar,nfreepar+1:end);
        cov_simsum = (cov_simsum+cov_simsum')/2;
        [~, notposdef] = cholcov(cov_simsum,0);
        if isnan(notposdef) || notposdef>0
          [L, DMC, P] = modchol_ldlt(cov_simsum); 
           cov_simsum = P'*L*DMC*L'*P;
        end
        proposal_mean = mean_theta + cov_thetasimsum * (cov_simsum \(summobs-mean_simsum));
        cov_blockedopt = zeros(nfreepar,nfreepar);
        ABCthreshold_temp = prctile(distance_all,alpha);
        ABCthreshold_old = ABCthreshold;
        ABCthreshold_new = min(ABCthreshold_temp,ABCthreshold_old); % ensure non-increasing thresholds
        if ABCthreshold_new == ABCthreshold_old
           ABCthreshold_new = 0.95*ABCthreshold_new;
        end
        id_blockedopt = find(distance_accepted < ABCthreshold_new);
        N0 = length(id_blockedopt);  % number of particles that satify the newest threshold
        if N0==0  % in this case there are no particles satifying the criterion above.
                  % we must exit
            fprintf('\n There are no particles satisfying the blockedopt criterion. We must leave...\n')
            return
        end
        weights_blockedopt = weights(id_blockedopt);
        normweights_blockedopt = weights_blockedopt/sum(weights_blockedopt);
        for jj=1:N0
              cov_blockedopt = cov_blockedopt + normweights_blockedopt(jj)*(ABCdraws(:,id_blockedopt(jj))-proposal_mean)*(ABCdraws(:,id_blockedopt(jj))-proposal_mean)';
        end
        cov_blockedopt = (cov_blockedopt+cov_blockedopt')/2;
        [~, notposdef] = cholcov(cov_blockedopt,0);
        if isnan(notposdef) || notposdef>0
          [L, DMC, P] = modchol_ldlt(cov_blockedopt); 
       %  perturbation = P'*L*DMC*L'*P - proposal_cov;
           cov_blockedopt = P'*L*DMC*L'*P;
       %   proposal_cov = nearestSPD(proposal_cov);
        end 
        proposal_mean =  proposal_mean';  % must be a row vector when passed to mvnrnd()
    case 'full_cond'
        mean_theta = ABCdraws * normweights';  % weighted mean across all accepted parameters
        mean_simsum = simsumm_accepted * normweights'; % weighted mean across all accepted summaries
        proposal_mean = zeros(nfreepar,numparticles);
        proposal_var = zeros(nfreepar,numparticles);
        all_cov_thetasimsum = (all_cov_thetasimsum+all_cov_thetasimsum')/2;
        [~, notposdef] = cholcov(all_cov_thetasimsum,0);
        if isnan(notposdef) || notposdef>0
          [L, DMC, P] = modchol_ldlt(all_cov_thetasimsum); 
           all_cov_thetasimsum = P'*L*DMC*L'*P;
        end
        for ii = 1:numparticles
           for par=1:nfreepar
              % notice the ~=par operator means "select all entries in the array except those corresponding to index par"
              proposal_mean(par,ii) = mean_theta(par) + all_cov_thetasimsum(par,1:end~=par) * (all_cov_thetasimsum(1:end~=par,1:end~=par) \ ([ABCdraws(1:end~=par,ii);summobs] - [mean_theta(1:end~=par);mean_simsum]));
              proposal_var(par,ii) = all_cov_thetasimsum(par,par) - all_cov_thetasimsum(par,1:end~=par) * (all_cov_thetasimsum(1:end~=par,1:end~=par) \ all_cov_thetasimsum(1:end~=par,par));
           end
        end
    case 'full_cond_opt'
        mean_theta = ABCdraws * normweights';  % weighted mean across all accepted parameters
        mean_simsum = simsumm_accepted * normweights'; % weighted mean across all accepted summaries
        var_opt_all = zeros(nfreepar,numparticles);
        proposal_mean = zeros(nfreepar,numparticles);
        ABCthreshold_temp = prctile(distance_all,alpha);
        ABCthreshold_old = ABCthreshold;
        ABCthreshold_new = min(ABCthreshold_temp,ABCthreshold_old); % ensure non-increasing thresholds
        if ABCthreshold_new == ABCthreshold_old
           ABCthreshold_new = 0.95*ABCthreshold_new;
        end
        id_opt = find(distance_accepted < ABCthreshold_new);
        N0 = length(id_opt);  % number of particles that satify the newest threshold
        if N0==0  % in this case there are no particles satifying the criterion above.
                  % we must exit
            fprintf('\n There are no particles satisfying the OLCM criterion. We must leave...\n')
            return
        end
        weights_opt = weights(id_opt);
        normweights_opt = weights_opt/sum(weights_opt);
        all_cov_thetasimsum = (all_cov_thetasimsum+all_cov_thetasimsum')/2;
        [~, notposdef] = cholcov(all_cov_thetasimsum,0);
        if isnan(notposdef) || notposdef>0
          [L, DMC, P] = modchol_ldlt(all_cov_thetasimsum); 
           all_cov_thetasimsum = P'*L*DMC*L'*P;
        end
        for ii = 1:numparticles
           for par=1:nfreepar
              % notice the ~=par operator means "select all entries in the array except those corresponding to index par"
              proposal_mean(par,ii) = mean_theta(par) + all_cov_thetasimsum(par,1:end~=par) * (all_cov_thetasimsum(1:end~=par,1:end~=par) \ ([ABCdraws(1:end~=par,ii);summobs] - [mean_theta(1:end~=par);mean_simsum]));
              for jj=1:N0
                  var_opt_all(par,ii) = var_opt_all(par,ii) + normweights_opt(jj)*(ABCdraws(par,id_opt(jj))-proposal_mean(par,ii))^2;
              end
           end
        end
    case 'standard'
        % compute the weighted covariance matrix C
        C = ABCdraws' - repmat(normweights * ABCdraws', numparticles, 1);   % subtract weighted mean                                                  % Remove mean (which is, also, weighted)
        C = C' * (C .* repmat(normweights', 1, nfreepar));                                                  % Weighted Covariance Matrix
        C = C ./ (1-sum(normweights.^2)); % the weighted covariance matrix
        C = 0.5 * (C + C');   % ensure symmetry
        Sigma = 2*C;    
     case 'olcm'
        ABCthreshold_temp = prctile(distance_all,alpha);
        ABCthreshold_old = ABCthreshold;
        ABCthreshold_new = min(ABCthreshold_temp,ABCthreshold_old); % ensure non-increasing thresholds
        if ABCthreshold_new == ABCthreshold_old
           ABCthreshold_new = 0.95*ABCthreshold_new;
        end
        id_olcm = find(distance_accepted < ABCthreshold_new);
        N0 = length(id_olcm);  % number of particles that satify the newest threshold
        if N0==0  % in this case there are no particles satifying the criterion above.
                  % we must exit
            fprintf('\n There are no particles satisfying the OLCM criterion. We must leave...\n')
            return
        end
        weights_olcm = weights(id_olcm);
        normweights_olcm = weights_olcm/sum(weights_olcm);
        cov_olcm_all = [];
        cov_olcm = zeros(nfreepar,nfreepar);
        for ii=1:numparticles
           for jj=1:N0
              cov_olcm = cov_olcm + normweights_olcm(jj)*(ABCdraws(:,id_olcm(jj))-ABCdraws(:,ii))*(ABCdraws(:,id_olcm(jj))-ABCdraws(:,ii))';
           end
          cov_olcm = (cov_olcm+cov_olcm')/2;
          cov_olcm_all = [cov_olcm_all,cov_olcm];
        end
    end
end
  
end

