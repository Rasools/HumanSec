function tempRxns=importTempRxnsFromExcel(fileName)
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
  % Rasool Saghaleyni, 2018-02-16
  [num,txt,~] = xlsread(fileName);
  tempRxns.Abbreviation = txt(2:end,1);
  tempRxns.Formula = txt(2:end,2);
  tempRxns.Comp = txt(2:end,3);
  tempRxns.Components = txt(2:end,4);
  tempRxns.Conditions = txt(2:end,5);
  tempRxns.GPR = txt(2:end,6);
  tempRxns.Pathway = txt(2:end,7);

  clear num txt
end
