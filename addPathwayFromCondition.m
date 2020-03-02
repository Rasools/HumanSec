function [newList,newList2]=addPathwayFromCondition(conditionName,listOfRxns,listOfRxnsNames,tempRxns)
%Function for adding the reactions of a given PathwayName given a condition
  %
  %
  %
  %
  %
  %
  %
  %
  % Rasool Saghaleyni, 2020-02-22
  newList = listOfRxns;
  newList2 = listOfRxnsNames;
  for i=1:length(tempRxns.Conditions)
    if strcmp(tempRxns.Conditions(i),conditionName) == 1
      newList=vertcat(newList,tempRxns.Formula(i));
      newList2=vertcat(newList2,tempRxns.Abbreviation(i));
    end
  end
end
