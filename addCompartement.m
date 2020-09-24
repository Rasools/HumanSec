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
% Rasool Saghaleyni,                              2020-02-22 
% don't add repetitive compartments               2020-09-23
  newComps = setdiff(comps, model.comps);
  [~, b] = ismember(newComps, comps);
  newCompNames = compNames(b);
  newModel = model;
  newModel.comps = vertcat(model.comps,newComps);
  newModel.compNames = vertcat(model.compNames,newCompNames);
end
