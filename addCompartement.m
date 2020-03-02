function newModel = addCompartement(model,comps,compNames,compOutside)
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
%
%
%
% Rasool Saghaleyni, 2020-02-22
  newModel=model;
  newModel.comps = vertcat(model.comps,comps);
  newModel.compNames = vertcat(model.compNames,compNames);
end
