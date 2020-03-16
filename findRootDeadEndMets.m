function [justConsumedMets, justProducedMets, id_cons, id_prod]=findRootDeadEndMets(model)
  %   This script will find metabolites which are just produced or just 
  %   consumed in the reactions. These metabolites cause
  %   deadEnd reactions.
  %
  %   model                 a model structure
  %   justConsumedMets      list of metabolites which are just consumed by
  %   the model
  %   justProducedMets      list of metabolites which are just produced by
  %   the model
  %
  %   Rasool Saghaleyni, 2020-03-07
  
  % Adding all boundry metabolites to the model
  model = addBoundaryMets(model);
  % remove all both directional reactions from model
  uniModel=removeReactions(model,model.rxns(model.rev==1));
  % generate a model just containing both directional reactions
  forwModel=removeReactions(model,model.rxns(model.rev==0));
  % fixing directionality of reaction to uni-directional
  forwModel.rev(:,1) = 0;
  % generate a model just containing both directional reactions
  revModel=removeReactions(model,model.rxns(model.rev==0));
  % change directionality of reaction to reverse
  revModel.S=-revModel.S;
  % fixing directionality of reaction to uni-directional
  revModel.rev(:,1) = 0;
  % adding removed reactions to main model
  newModel=mergeModels({uniModel,forwModel,revModel});
  % finding metabolites which are just produced or consumed in all
  % reactions
  for i=1:size(newModel.mets)
    if size(find(newModel.S(i,:)<0),2)==0
        all_mets(i,1)={'prod'};
    elseif size(find(newModel.S(i,:)>0),2)==0
        all_mets(i,1)={'cons'};
    end
  end
  % if you add this metabolites to the model using addExchangeRxns function, there will
  % not be any deadEnd reaction in the model.
  id_cons = find(strcmp(all_mets,'cons'));
  id_prod = find(strcmp(all_mets,'prod'));
  justConsumedMets=newModel.mets(strcmp(all_mets,'cons'));
  justProducedMets=newModel.mets(strcmp(all_mets,'prod'));
end
