function [newListRxns, newListNames, newListGPRs, newListComps]=addCanonicalRxns(listOfRxns,listOfRxnNames,listOfGPRs,listOfComps,tempRxns)
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
  newListRxns = listOfRxns;
  newListNames = listOfRxnNames;
  newListGPRs = listOfGPRs;
  newListComps = listOfComps;
  %Add  canonical reactions
  [newListRxns, newListNames] = addPathway('Canonical',listOfRxns,listOfRxnNames,tempRxns);
  %Add canonical GPRs and compartements
  r = {};
  n = {};
  [r,n] = addPathway('Canonical',r,n,tempRxns);
  newGPRs = getGPRsFromRxnNames(n,tempRxns);
  newListGPRs = vertcat(newListGPRs, newGPRs);

  r = {};
  n = {};
  [r,n] = addPathway('Canonical',r,n,tempRxns);
  newComps = getCompsFromRxnNames(n,tempRxns);
  newListComps = vertcat(newListComps, newComps);
end
