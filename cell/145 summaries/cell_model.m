function data = cell_model(bigtheta,covariates,~)

Yinit = covariates.Yinit;
Yobs = covariates.Yobs;
obstimes = covariates.obstimes;

tau = 1/24; %Duration of each time step
num_obs = size(Yobs,3);

sampling_rate = max(obstimes)/num_obs; %Time between consecutive observations
sim_iters = sampling_rate/tau;

[rows, cols] = find(Yinit); %Finds non-zero indices (positions of cells)
rows = int32(rows-1);
cols = int32(cols-1);
X = int32(Yinit); %signed 32 bit integers

S = simulate_cell(bigtheta,X,rows,cols,num_obs,sim_iters);

data.Y = S;
data.Yinit = Yinit;

end

