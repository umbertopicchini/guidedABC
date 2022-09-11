function x = gk_model(bigtheta,~,numsim)

if nargin == 1
    numsim=1;
end

if nargin < 2
    numsim=1;
end


B = bigtheta(1);
g = bigtheta(2);
k = bigtheta(3);
alpha = bigtheta(4);
As = bigtheta(5:end);


% we are going to simulate numdatasets datasets each having length nobs. Each 
% dataset is generated with a different A coefficient (As(ii))
numdatasets = 20;
nobs = 1000;  % number of observations for each data set
x = zeros(nobs,numdatasets);

for ii=1:numdatasets
   z = randn(nobs,1);
   % stack the datasets on top of eachother
   x(:,ii) = As(ii) + B * (1 + 0.8 * (1-exp(-g*z))./(1+exp(-g*z))) .* (1 + z.^2).^k .* z;
end

%x = x(:);

end

