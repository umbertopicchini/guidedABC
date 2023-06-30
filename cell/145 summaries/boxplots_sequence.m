iterations =5;
paramplot = 1;  % the index of the parameter you want boxplots for
attempt = 1;
numparticles = 1000;
theta = [];
g = [];
for ii = 1:iterations
    filename=sprintf('ABCdraws_stage%d_attempt%d_copula0_marg5_nu5_numarg5.0_part1000.dat',ii,attempt);
    theta_temp = load(filename);
    filename=sprintf('weights_stage%d_attempt%d_copula0_marg5_nu5_numarg5.0_part1000.txt',ii,attempt);
    weights = load(filename);
    index = randsample(numparticles,numparticles,true,weights); % resample according to weights
    theta_temp = theta_temp(:,index); 
    theta_temp = theta_temp';
    theta_par = theta_temp(:,paramplot);
    theta = [theta;theta_par];
    g = [g; ii*ones(size(theta_par))];
end
%figure
boxplot(theta,g)
xticks([1:5])
xticklabels({'1','2','3','4','5'})
%xticklabels({'1','3','5','7','9','11'})
