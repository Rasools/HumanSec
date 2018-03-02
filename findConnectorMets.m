function connectors=findConnectorMets(model, addedRxns, protein)
  % for finding connector metabolites between protein secretory pathway and other pathways. Maybe useful for gap filling purposes.
  %
  %
  % model       A genome scale metaboic model like hmr2. This model should not have secretory pathway.
  % addedRxns   a data structure produced by humanSec function.
  %
  % connectors    an string array of connector metabolites
  %
  %
  %
  %
  %
  %
  %
  % Rasool Saghaleyni, 2018-02-16
  outModel   = addCompartmentalizedMetNames(model);
  metabolites = parseRxnEqu(addedRxns.rxnFormula);
  connectors = intersect(outModel.metNamesC,metabolites);
  a = contains(connectors,protein);
  b = find(a==0);
  connectors = connectors(b);
end
