function newModel = copyRxn(model,refModel,rxnList)
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
  rxns_ind = find(ismember(refModel.rxns,rxnList));
  rxnsToAdd.rxns = refModel.rxns(rxns_ind);
  rxnsToAdd.equations = constructEquations(refModel,rxnList);
  rxnsToAdd.rxnNames = refModel.rxnNames(rxns_ind);
  rxnsToAdd.lb = refModel.lb(rxns_ind);
  rxnsToAdd.ub = refModel.ub(rxns_ind);
  rxnsToAdd.c = refModel.c(rxns_ind);
  rxnsToAdd.eccodes = refModel.eccodes(rxns_ind);
  rxnsToAdd.subSystems = refModel.subSystems(rxns_ind);
  rxnsToAdd.grRules = refModel.grRules(rxns_ind);
  %rxnsToAdd.rxnMiriams = refModel.rxnMiriams(rxns_ind);
  rxnsToAdd.rxnComps = refModel.rxnComps(rxns_ind);
  %rxnsToAdd.rxnNotes = refModel.rxnNotes(rxns_ind);
  %rxnsToAdd.rxnReferences = refModel.rxnReferences(rxns_ind);
  %rxnsToAdd.rxnConfidenceScores = refModel.rxnConfidenceScores(rxns_ind);
  newModel = addRxns(model,rxnsToAdd,3);
end
