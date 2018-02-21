function [new_model, newMets] = addSecretionReactions(model, rxnFile, rxnNames, rxnGPRs, rxnComps)
%
% addSecretionReactions Expands a metabolic model to include secretory
% reactions for a given protein
%
%INPUTS
% model = COBRA model structure
% rxnFile = String of .txt file containing reactions list
% rxnNames = String of .txt file containing reactions names
% rxnGPRs = String of .txt file containing reactions GPRs
% geneList = String of .txt file containing genes to be added to model

%OUTPUT
% new_model = Expanded model with Metabolic and Secretory capabilities

%EXAMPLES
% for CHO:
% choEPO = addSecretionReactions(model, 'CHO_P01588_reactions.txt', 'CHO_P01588_names.txt', 'CHO_P01588_GPRs.txt', 'geneIDs_CHO.txt')

% for HUMAN:
% humanEPO = addSecretionReactions(recon1, 'HUMAN_P01588_reactions.txt', 'HUMAN_P01588_names.txt', 'HUMAN_P01588_GPRs.txt', 'geneIDs_HUMAN.txt')

% NOTE: You need to run 'initCobraToolbox' first!

% Jahir M. Gutierrez 08/24/2016

warning off all

new_model = model;
name = {};
reaction = {};
gpr = {};
%-----------------------
% Import reaction names into a cell array
fileID = fopen(rxnNames);
i=1;
while feof(fileID)==0
    name{i,1} = fgetl(fileID);
    i = i+1;
end
fclose(fileID);
name=strrep(name,'[','');
name=strrep(name,']','');
%-----------------------
% Import reaction into a cell array
fileID = fopen(rxnFile);
i=1;
while feof(fileID)==0
    reaction{i,1} = fgetl(fileID);
    i=i+1;
end
fclose(fileID);
reaction=strrep(reaction,'[e]','[s]');
%---------------------
% Import reaction GPRs into a cell array
fileID = fopen(rxnGPRs);
i=1;
while feof(fileID)==0
    gpr{i,1} = fgetl(fileID);
    i=i+1;
end
fclose(fileID);
%_--------------------
fileID = fopen(rxnComps);
i=1;
while feof(fileID)==0
    comp{i,1} = fgetl(fileID);
    i=i+1;
end
fclose(fileID);
%--------------------- add new genes
gprList=gpr;
gprList=strrep(gprList,'(','');
gprList=strrep(gprList,')','');
gprList=strrep(gprList,';','');
gprList=strrep(gprList,' or ',' ');
gprList=strrep(gprList,' and ',' ');

for i=1:length(gprList)
    genes(1,:)=strsplit(char(gprList(i)));
    if i==1
        geneList=genes(:);
    else
        geneList=vertcat(geneList,genes(:));
    end
    clear genes
end


geneList=unique(geneList);

genesToAdd.genes = setdiff(geneList,new_model.genes);
genesToAdd.genes = genesToAdd.genes(~cellfun('isempty', genesToAdd.genes));
new_model=addGenes(new_model,genesToAdd);
%-------------------correcting reactions filrmula to can be added to model
reaction=strrep(reaction,'  ',' ');
reaction=strrep(reaction,'=> +','=>');
%-----------------------new metaboilites for function output
rxnList=reaction;
rxnList=strrep(rxnList,' + ',' ');
rxnList=strrep(rxnList,' => ',' ');
rxnList=strrep(rxnList,' <=> ',' ');
rxnList=strrep(rxnList,'=>',' ');

for i=1:length(rxnList)
    mets(1,:)=strsplit(char(rxnList(i)));
    if i==1
        metList=mets(:);
    else
        metList=vertcat(metList,mets(:));
    end
    clear mets
end

metList=unique(metList);
fullmetList = metList(~cellfun('isempty', metList));
indC = strfind(fullmetList, '[');
fullmetList = fullmetList(find(not(cellfun('isempty', indC))));
outModel=addCompartmentalizedMetNames(model);
newMets=setdiff(fullmetList,outModel.metNamesC);

%-----------------------new reactions

rxnsToAdd.rxns = name(ismember(name,model.rxns)==0);
rxnsToAdd.equations = reaction(ismember(name,model.rxns)==0);
rxnsToAdd.rxnNames = name(ismember(name,model.rxns)==0);
rxnsToAdd.grRules = gpr(ismember(name,model.rxns)==0);

comp = comp(ismember(name,model.rxns)==0);
%comp = cellfun(@(x)str2double(x), comp);

% Add protein secretion reactions to HMR2-drived model and modify grule Matrix
new_model = addRxns(new_model,rxnsToAdd,3,[],true);
% Add compartments for new reactions
new_model.rxnComps(length(model.rxns)+1:length(new_model.rxns)) = comp(1:end);

warning on all
