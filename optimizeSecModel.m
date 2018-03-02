function [new_secModel, sol]=optimizeSecModel(secModel, protein)
  %
  %
  %
  %
  %
  %
  %
  %
  %
  %
  %
  %
  %
  %
  % Rasool Saghaleyni, 2018-02-16
  rxnName = strcat(protein ,'_Final_demand');
  new_secModel = secModel;
  new_secModel = setParam(new_secModel, 'obj', rxnName, [1]); % set intrested protein production reaction as objective function
  sol = solveLP(new_secModel,1);
end
