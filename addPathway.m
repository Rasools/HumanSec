function [newList, newList2]=addPathway(pathwayName,listOfRxns,listOfRxnsNames,tempRxns)
  %Function for adding the reactions of a given PathwayName
  %
  %   pathwayName         a cell string containing name of pathway
  %   listOfRxns          a cell array which contains list of formulas for current reactions
  %   listOfRxnsNames     a cell array which contains list of names for current reactions
  %   tempRxns
  %   tempRxns            the tempRxns structure can have the following fields:
  %
  %       Abbreviation
  %       Components
  %       Comp
  %       Conditions
  %       Formula
  %       GPR
  %       Pathway
  %
  %   newList
  %   newList2
  %
  %
  %
  % Rasool Saghaleyni, 2018-02-16
  newList = listOfRxns;
  newList2 = listOfRxnsNames;
  for i=1:numel(tempRxns.Pathway)
    if strcmp(tempRxns.Pathway(i),pathwayName) == 1
      newList = vertcat(newList,tempRxns.Formula(i));
      newList2 = vertcat(newList2,tempRxns.Abbreviation(i));
    end
  end
end
