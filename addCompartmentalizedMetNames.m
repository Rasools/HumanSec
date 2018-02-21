function outModel   = addCompartmentalizedMetNames(inModel)
  % addCompartmentalizedMetNames
  % This function adds to the model structure a field called metNamesC
  % identical to metNames but with the compartment appended at the end.
  %
  % INPUT
  % inModel:  a genome-scale metabolic model structure. Note: if a field
  %           called metNamesC is already in the model, it will be
  %           overwritten.
  %
  % OUTPUT
  % outModel:  a genome-scale metabolic model structure with a metNamesC
  %            field.
  %
  % Usage:  outModel = addCompartmentalizedMetNames(inModel)
  % 2013-07-17 Francesco Gatto
  outModel = inModel;
  if isfield(outModel,'metNamesC')
    outModel = rmfield(outModel,'metNamesC');
  end
  mets = outModel.metNames;
  nMets = length(mets);
  comps = outModel.comps(outModel.metComps);
  metNamesC = cellfun(@(a,b,c,d) [a,b,c,d],mets,repmat({'['},nMets,1),comps,repmat({']'},nMets,1),'uni',false);
  outModel.metNamesC = metNamesC;
end
