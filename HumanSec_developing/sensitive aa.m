load('/Users/rasools/Box Sync/phd/ppf/TPM_v03/secretion_network_modeling/HMR2s/EPO_Snsitivity_analysis/Sensitivity results.mat');

for i=1:size(length, 1)
	aacomp = aacount(char(seq(i)));
	comp(i,:) = struct2table(aacomp);
end
comp_ar = table2array(comp);

for i=1:size(length, 1)
	comp_nor(i,:) = comp_ar(i,:)/length(i);
end

comp_nor(find(isnan(si)),:) = [];
si(find(isnan(si)),:) = [];

for i=1:20
	cor(1,i) = corr(comp_nor(:,i),si);
end