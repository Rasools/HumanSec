function [fixedModel, fixingAddedRxns] = secModelGapFilling2(secModel,refModel,protein,psim,tempRxns,addedRxns)
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
  objRxn = strcat(protein,'_Final_demand');
  outMet = parseRxnEqu(addedRxns.rxnFormula(find(ismember(addedRxns.rxnsNames,objRxn))));
  [secRefModel,~] = humanSec(refModel, protein, psim, tempRxns);
  [~, sol] = optimizeSecModel(secRefModel, protein);
  [in_ex, in_ex_ind] = getExchangeRxns(secRefModel,'in');
  [out_ex, out_ex_ind] = getExchangeRxns(secRefModel,'out');
  in_ex_havFlx = in_ex(find(sol.x(in_ex_ind)));
  out_ex_havFlx = out_ex(find(sol.x(out_ex_ind)));

  in_mets = parseRxnEqu(constructEquations(secRefModel,in_ex_havFlx));
  out_mets = parseRxnEqu(constructEquations(secRefModel,out_ex_havFlx));

%   taskStructure.inputs = in_mets;
%   taskStructure.LBin(1:numel(taskStructure.inputs),1) = 0;
%   taskStructure.UBin(1:numel(taskStructure.inputs),1) = 1000;

  taskStructure.inputs = {'O2[s]';'glucose[s]';'NH3[s]';'H2O[s]';'arginine[s]';'histidine[s]';'lysine[s]';'methionine[s]';'phenylalanine[s]';'tryptophan[s]';'tyrosine[s]';'alanine[s]';'glycine[s]';'serine[s]';'threonine[s]';'aspartate[s]';'glutamate[s]';'asparagine[s]';'glutamine[s]';'isoleucine[s]';'leucine[s]';'proline[s]';'valine[s]';'cysteine[s]'};
  taskStructure.LBin(1:numel(taskStructure.inputs),1) = 0;
  taskStructure.UBin(1:numel(taskStructure.inputs),1) = 1000;

%   taskStructure.outputs = out_mets;
%   taskStructure.LBout(1:numel(taskStructure.outputs)-1,1) = 0;
%   taskStructure.LBout(numel(taskStructure.outputs),1) = 0.0001;
%   taskStructure.UBout(1:numel(taskStructure.outputs),1) = 1000;

  taskStructure.outputs = vertcat({'CO2[s]';'H2O[s]';'H2S[s]';'urate[s]'},outMet);
  taskStructure.LBout = [0;0;0;0;0.001];
  taskStructure.UBout(1:numel(taskStructure.outputs),1) = 1000;

  taskStructure.id = {'sec'};
  taskStructure.description = strcat({'production of '},protein);
  taskStructure.shouldFail = 0;
  taskStructure.printFluxes = 0;
  taskStructure.comments = '';
  taskStructure.equations = {};
  taskStructure.LBequ = [];
  taskStructure.UBequ = [];
  taskStructure.changed = {};
  taskStructure.LBrxn = [];
  taskStructure.UBrxn = [];

  [fixedModel, fixingAddedRxns] = fitTasks(secModel,refModel,[],[],[],taskStructure,[]);
  fixingAddedRxns = refModel.rxns(find(fixingAddedRxns));

  neededRxns = {'HMR_9034';'HMR_9038';'HMR_9039';'HMR_9040';'HMR_9041';'HMR_9042';'HMR_9043';'HMR_9044';'HMR_9045';'HMR_9046';'HMR_9047';'HMR_9048';'HMR_9058';'HMR_9061';'HMR_9062';'HMR_9063';'HMR_9064';'HMR_9065';'HMR_9066';'HMR_9067';'HMR_9068';'HMR_9069';'HMR_9070';'HMR_9071';'HMR_9073';'HMR_9075';'HMR_9103'};
  gapRxns = neededRxns(find(~ismember(neededRxns,fixedModel.rxns)));
  fixedModel = copyRxn(fixedModel,refModel,gapRxns);
  fixingAddedRxns = vertcat(fixingAddedRxns,gapRxns);
end
