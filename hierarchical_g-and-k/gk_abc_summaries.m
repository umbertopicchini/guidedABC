function sim_summ = gk_abc_summaries(x)

 sim_summ = zeros(9,size(x,2));
% 
% for ii=0:8
%     sim_summ(ii+1) = quantile(x,ii/8) ;
% end

for ii=1:size(x,2)
   sim_summ(:,ii) = quantile(x(:,ii),[0:8]/8);
end

sim_summ = sim_summ(:);

end

