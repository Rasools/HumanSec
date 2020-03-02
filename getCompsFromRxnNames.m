function comp_list=getCompsFromRxnNames(listOfRxnsNames, tempRxns)
  %Function for creating a list of components given a list of Rxn names
  %
  %
  %
  %
  %
  %
  %
  %
  % Rasool Saghaleyni, 2018-02-16
  comp_list = {};
  %print listOfRxnsNames
  for i = 1:numel(listOfRxnsNames)
    if ismember(listOfRxnsNames(i),tempRxns.Abbreviation) == 1
      comp_list(i,1) = tempRxns.Comp(find(ismember(tempRxns.Abbreviation,listOfRxnsNames(i))));
    else
      comp_list(i,1) = {''};
    end
  end
end
