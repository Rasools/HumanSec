function GPR_list=getGPRsFromRxnNames(listOfRxnsNames, tempRxns)
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
  GPR_list = {};
  %print listOfRxnsNames
  for i=1:numel(listOfRxnsNames)
    if ismember(listOfRxnsNames(i),tempRxns.Abbreviation) == 1
      GPR_list(i,1) = tempRxns.GPR(find(ismember(tempRxns.Abbreviation,listOfRxnsNames(i))));
    else
      GPR_list(i,1)={''};
    end
  end
end
