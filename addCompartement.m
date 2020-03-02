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
%
  newModel=model;
  newModel.comps = vertcat(model.comps,comps);
  newModel.compNames = vertcat(model.compNames,compNames);
end
