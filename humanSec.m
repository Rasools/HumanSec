function [secModel, addedRxns]=humanSec(model, protein, psim, tempRxns)
  % humanSec
  %   This script contains functions to create a list of protein-specific reactions
  %   for the protein secretion pathway in HMR2s.
  %
  %   model           a model structure, could be HMR2 or a HMR2 drived model
  %   protein         the Uniprot ID of your interested protein. if protein ID is not in PSIM_HUMAN you should create
  %                   a specific PSIM for your interested protein
  %   psim            the PSIM structure can have the following fields:
  %
  %       entry
  %       name
  %       length
  %       mass
  %       sp
  %       dsb
  %       gpi
  %       ng
  %       og
  %       tmd
  %       location
  %       sequence
  %
  %   tempRxns        the tempRxns structure can have the following fields:
  %
  %       Abbreviation
  %       Components
  %       Comp
  %       Conditions
  %       Formula
  %       GPR
  %       Pathway
  %
  %   secModel        updated model structure
  %   addedRxns       a structure array contains reactions which added to model to generate secModel
  %
  %       rxnsNames
  %       rxnFormula
  %       rxnGPRs
  %       rxnComps
  %       newMets
  %       newGens
  %
  %
  %
  %   Usage: [secModel, addedRxns]=humanSec(model, proToAdd)
  %
  %   Rasool Saghaleyni, 2018-02-16
  %                      2019-02-07   
  %                      2020-03-01
  
  p_row = find(ismember(psim.entry,protein)); %protein row in psim
  sequence = char(psim.sequence(p_row));
  length = psim.length(p_row);
  mass = psim.mass(p_row);
  Kv = 0.7;
  V = mass * 1.21 / 1000.0; % Protein Volume in nm^3
  clathrin_coeff = floor(round(29880.01 * Kv / V)); % Number of proteins per clathrin vesicle
  copi_coeff = floor(round(143793.19 * Kv / V));
  copii_coeff = floor(round(268082.35 * Kv / V));
  %Prepare vectors that will store reactions and components
  rxns={};
  rxnNames={};
  %Add translation reaction
  translationFormula = translationRxn(psim, protein);
  rxns(1) = cellstr(translationFormula);
  rxnNames(1) = {'TRANSLATION_protein'};
  secModel = model;
  secModel.description = char('human protein-specific secretory model');
  secModel.id = char(strcat('HumanSec_',protein));
  % Add secretory compartements
  comps = {'sv';'cv';'rm';'gm';'pm'};
  compNames = {'secretory vesicle';'clathrin vesicle';'ER membrane';'Golgi membrane';'Plasma membrane'};
  compOutside = {'c';'c';'c';'c';'s'};
  secModel = addCompartement(secModel,comps,compNames,compOutside);
  %------------------------------------------------------------------------------
  %If protein doesnt have signal peptide then ignore
  if psim.sp(p_row) == 0
    GPRs = getGPRsFromRxnNames(rxnNames, tempRxns);
    Comps = getCompsFromRxnNames(rxnNames, tempRxns);
    %Change the reaction names and their formulas to include the UniProt Identifier
    for i = 1:numel(rxns)
      rxns(i) = insert_prot_name_in_rxnFormula((rxns(i)),protein);
    end
    for i = 1:numel(rxns)
      rxnNames(i) = insert_prot_name_in_rxnName(rxnNames(i),protein);
    end
    rxns = vertcat(rxns,strcat(protein, {'[c] => '}));
    rxnNames = vertcat(rxnNames,strcat(protein ,{'_Final_demand'}));
    GPRs = vertcat(GPRs,{''});
    Comps = vertcat(Comps,9);
  else %Translocate protein if it has signal peptide
    if length <= 160
      [rxns,rxnNames] = addPathway('Post-translational Translocation',rxns,rxnNames,tempRxns);
      if strcmp(psim.location(p_row), '[s]') || strcmp(psim.location(p_row), '')
        [rxns,rxnNames] = addPathway('Post-translational Translocation (Secretory protein)',rxns,rxnNames,tempRxns);
      end
      if strcmp(psim.location(p_row), '[pm]') || psim.og(p_row) ~= 0
        [rxns,rxnNames] = addPathway('Post-translational Translocation (Tail anchored membrane protein)',rxns,rxnNames,tempRxns);
      end
    else
      [rxns,rxnNames] = addPathway('Translocation',rxns,rxnNames,tempRxns);
    end
    number_BiP = length/40; %Number of BiPs depends on protein length http://www.cshperspectives.com/content/5/5/a013201.full
    for i = 1:numel(rxns)
      rxns(i) = strrep(rxns(i), '!',num2str(number_BiP));
    end
    connector = 'XXX[r]';
    %------------------------------------------------------------------------------
    %Add Disulphide Bond Reactions
    if psim.dsb(p_row) ~= 0
      number_DSB = psim.dsb(p_row);
      DSBrxns = {};
      DSBrxnNames = {};
      DSBrxns = vertcat(DSBrxns,strcat(connector ,{' => XXX_preDSB[r]'}));
      DSBrxnNames = vertcat(DSBrxnNames,{'Start_DSB'});
      [DSBrxns,DSBrxnNames]=addPathwayFromCondition('DSB>0',DSBrxns,DSBrxnNames,tempRxns);
      for i = 1:numel(DSBrxns)
        if contains(DSBrxns(i),'?') == 1
          if number_DSB == 1
            DSBrxns(i) = strrep(DSBrxns(i),'? ','');
          else
            DSBrxns(i) = strrep(DSBrxns(i),'?',num2str(number_DSB));
          end
        end
      end
      connector = 'XXX_DSB[r]';
      rxns = vertcat(rxns,DSBrxns);
      rxnNames = vertcat(rxnNames,DSBrxnNames);
    end
    %------------------------------------------------------------------------------
    %Add GPI reactions
    if psim.gpi(p_row) == 1
      rxns = vertcat(rxns,strcat(connector,{' => XXX_preGPI[r]'}));
      rxnNames = vertcat(rxnNames,{'Start_GPI'});
      [rxns,rxnNames] = addPathwayFromCondition('GPI=1',rxns,rxnNames,tempRxns);
      connector = 'XXX-dgpi_hs[r]';
    end
    %------------------------------------------------------------------------------
    %Add N-glycosylation reactions
    if psim.ng(p_row) ~= 0 %protein has N-glycans!
      rxns = vertcat(rxns,strcat(connector, {' => XXX_preNG[r]'}));
      rxnNames = vertcat(rxnNames,{'Start_NG'});
      number_Nglycans = psim.ng(p_row); %Get number of N-Glycans
      NGlyrxns = {};
      NGlyrxnNames = {};
      [NGlyrxns,NGlyrxnNames] = addPathwayFromCondition('NG>0',NGlyrxns,NGlyrxnNames,tempRxns);
      [NGlyrxns,NGlyrxnNames] = addPathway('Golgi processing N',NGlyrxns,NGlyrxnNames,tempRxns);
      for i = 1:numel(NGlyrxns) %Change the '?' for the number of N-glycans
        if strcmp(NGlyrxns(i),{'XXX-M5-unfold-UBIQP[c] + ? H2O[c] + RAD23A[c] => XXX-unfold-UBIQP-RAD23A[c] + ? N-acetylglucosamine[c] + ? mannose[c]'}) == 1
          h2o = num2str(7*number_Nglycans);
          acgam = num2str(2*number_Nglycans);
          man = num2str(5*number_Nglycans);
          NGlyrxns(i) = strcat({'XXX-M5-unfold-UBIQP[c] + '}, h2o,{' H2O[c] + RAD23A[c] => XXX-unfold-UBIQP-RAD23A[c] + '},acgam,{' N-acetylglucosamine[c] + '},man,{' mannose[c]'});
        end
        if contains(NGlyrxns(i),'?') == 1
          if number_Nglycans == 1
            NGlyrxns(i) = strrep(NGlyrxns(i),'? ','');
          else
            NGlyrxns(i) = strrep(NGlyrxns(i),'?',num2str(number_Nglycans));
          end
        end
      end
      rxns = vertcat(rxns,NGlyrxns);
      rxnNames = vertcat(rxnNames,NGlyrxnNames);
      connector = 'XXX-M8B[r]';
    end
    %------------------------------------------------------------------------------
    %Add  COPII reactions
    if strcmp(connector,'XXX-M8B[r]') == 1
      copii_rxns = {};
      copii_names = {};
      [copii_rxns, copii_names] = addPathway('COPII_NG',copii_rxns,copii_names,tempRxns);
      for i = 1:numel(copii_rxns)
        copii_rxns(i) = strrep(copii_rxns(i),'!',num2str(copii_coeff));
      end
      rxns = vertcat(rxns,copii_rxns);
      rxnNames = vertcat(rxnNames,copii_names);
      connector = 'XXX-M3-GN2[g]';
      %-------------- ADD IF PROTEIN = EPO (HUMAN)---------------------------------
      % IMPORTANT! => Reactions in this pathway assume EPO has NG = 3 and Og = 1 (see PSIM)
      if strcmp(protein,'P01588') == 1
        rxns = vertcat(rxns,{'P01588[c] => EPO_Human[c]'});
        rxnNames = vertcat(rxnNames,{'make_EPO'});
        [rxns,rxnNames] = addPathway('Golgi processing (EPO specific)',rxns,rxnNames,tempRxns);
        connector = 'XXX-M3-GN4-GL4-NA4-F[g]';
        protein = 'EPO_Human';
      end
    elseif strcmp(connector,'XXX-dgpi_hs[r]') == 1
      copii_rxns = {};
      copii_names = {};
      [copii_rxns, copii_names] = addPathway('COPII_GPI',copii_rxns,copii_names,tempRxns);
      for i = 1:numel(copii_rxns)
        copii_rxns(i) = strrep(copii_rxns(i),'!',num2str(copii_coeff));
      end
      rxns = vertcat(rxns,copii_rxns);
      rxnNames = vertcat(rxnNames,copii_names);
      connector = 'XXX-dgpi_hs[g]';
    elseif strcmp(connector,'XXX_DSB[r]') == 1
      [rxns,rxnNames] = addPathway('COPII_DSB',rxns,rxnNames,tempRxns);
      copii_rxns = {};
      copii_names = {};
      [copii_rxns, copii_names] =  addPathway('COPII-canonical',copii_rxns,copii_names,tempRxns);
      for i = 1:numel(copii_rxns)
        copii_rxns(i) = strrep(copii_rxns(i),'!',num2str(copii_coeff));
      end
      rxns = vertcat(rxns,copii_rxns);
      rxnNames = vertcat(rxnNames,copii_names);
      connector = 'XXX[g]';
    elseif strcmp(connector,'XXX[r]') == 1
      [rxns,rxnNames] = addPathway('COPII-normal',rxns,rxnNames,tempRxns);
      copii_rxns = {};
      copii_names = {};
      [copii_rxns, copii_names] =  addPathway('COPII-canonical',copii_rxns,copii_names,tempRxns);
      for i = 1:numel(copii_rxns)
        copii_rxns(i) = strrep(copii_rxns(i),'!',num2str(copii_coeff));
      end
      rxns = vertcat(rxns,copii_rxns);
      rxnNames = vertcat(rxnNames,copii_names);
      connector = 'XXX[g]';
    end
    %------------------------------------------------------------------------------
    %Add O-glycosylation reactions
    if psim.og(p_row) ~= 0 %Has O-glycans?
      rxns = vertcat(rxns,strcat(connector,{' => XXX_preOG[g]'}));
      rxnNames = vertcat(rxnNames,{'Start_OG'});
      number_Oglycans = psim.og(p_row); %Get number of O-Glycans
      OGlyrxns = {};
      OGlyrxnNames = {};
      [OGlyrxns,OGlyrxnNames] = addPathwayFromCondition('OG>0',OGlyrxns,OGlyrxnNames,tempRxns);
      for i = 1:numel(OGlyrxns) %Change the '?' for the number of O-glycans
        if contains(OGlyrxns(i),'?') == 1
          if number_Oglycans == 1
            OGlyrxns(i) = strrep(OGlyrxns(i),'? ', '');
          else
            OGlyrxns(i) = strrep(OGlyrxns(i),'?',num2str(number_Oglycans));
          end
        end
      end
      rxns = vertcat(rxns,OGlyrxns);
      rxnNames = vertcat(rxnNames,OGlyrxnNames);
      connector = 'XXX-Core2[g]';
    end
    %------------------------------------------------------------------------------
    %determining final location of the protein
    if strcmp(psim.location(p_row),'[g]') == 1 || strcmp(psim.location(p_row),'[gm]') == 1
      %Stay in Golgi if protein is localized there
      location = psim.location(p_row);
      if ismember('SP_degradation',rxnNames) == 1
        SPaas = aacount(sequence(1:22)); %Amino acids in signal peptide assuming length is 22 on average
        rxns(find(ismember(rxnNames,'SP_degradation'))) = substitute_AAs_count(rxns(find(ismember(rxnNames,'SP_degradation'))),SPaas);
      end
      if ismember('Ubiquitination_degradation',rxnNames) == 1
        new_aas = aacount(sequence(1:22));
        rxns(find(ismember(rxnNames,'Ubiquitination_degradation'))) = strrep(rxns(find(ismember(rxnNames,'Ubiquitination_degradation'))),'!', '?');
        rxns(find(ismember(rxnNames,'Ubiquitination_degradation'))) = substitute_AAs_count(rxns(find(ismember(rxnNames,'Ubiquitination_degradation'))),new_aas);
        rxns(find(ismember(rxnNames,'Ubiquitination_degradation'))) = strrep(rxns(find(ismember(rxnNames,'Ubiquitination_degradation'))),'?',num2str(length-22));
      end
      %Add GPRs
      GPRs = getGPRsFromRxnNames(rxnNames, tempRxns);
      Comps = getCompsFromRxnNames(rxnNames, tempRxns);
      %Change the reaction names and their formulas to include the UniProt Identifier
      for i = 1:numel(rxns)
        rxns(i) = insert_prot_name_in_rxnFormula((rxns(i)),protein);
      end
      for i = 1:numel(rxnNames)
        rxnNames(i) = insert_prot_name_in_rxnName(rxnNames(i),protein);
      end
      rxns = vertcat(rxns,strcat(protein,location, {' => '}));
      rxnNames = vertcat(rxnNames,strcat(protein, {'_Final_demand'}));
      GPRs = vertcat(GPRs,{});
      Comps = vertcat(Comps,0);
      [rxns,rxnNames,GPRs, Comps] = addCanonicalRxns(rxns,rxnNames,GPRs, Comps, tempRxns);
      %------------------------------------------------------------------------------
      %Add COPI
    elseif strcmp(psim.location(p_row),'[r]') == 1 || strcmp(psim.location(p_row),'[rm]') == 1
      rxns = vertcat(rxns,strcat(connector, {' => XXX_preCOPI[g]'}));
      rxnNames = vertcat(rxnNames,{'Start_COPI'});
      copi_rxns = {};
      copi_names = {};
      [copi_rxns, copi_names] = addPathway('COPI',copi_rxns,copi_names,tempRxns);
      for i = 1:numel(copi_rxns)
        copi_rxns(i) = strrep(copi_rxns(i),'!',num2str(copi_coeff));
      end
      rxns = vertcat(rxns,copi_rxns);
      rxnNames = vertcat(rxnNames,copi_names);
      connector = 'XXX_mature[r]';
      location = psim.location(p_row);
      if strcmp(psim.location(p_row),'[r]') == 1
        rxns = vertcat(rxns,strcat(connector, {' => '}));
        rxnNames = vertcat(rxnNames,strcat(protein,{'_Final_demand'}));
      elseif strcmp(psim.location(p_row),'[rm]') == 1
        rxns = vertcat(rxns,strcat(connector,{' => XXX_mature'},location));
        rxnNames = vertcat(rxnNames,strcat({'Final_location_'},location));
        rxns = vertcat(rxns,strcat({'XXX_mature'},location,{' => '}));
        rxnNames = vertcat(rxnNames,strcat(protein,{'_Final_demand'}));
      end
      %------------------------------------------------------------------------------
      %Add Clathrin vesicles
    elseif  strcmp(psim.location(p_row),'[p]') == 1 || strcmp(psim.location(p_row),'[l]') == 1
      rxns = vertcat(rxns,strcat(connector, {' => XXX-preClathrin[g]'}));
      rxnNames = vertcat(rxnNames,{'Start_Clathrin_vesicle'});
      clath_rxns = {};
      clath_names = {};
      [clath_rxns, clath_names] = addPathway('Clathrin vesicles',clath_rxns,clath_names,tempRxns);
      for i = 1:numel(clath_rxns)
        clath_rxns(i) = strrep(clath_rxns(i),'!',num2str(clathrin_coeff));
      end
      rxns = vertcat(rxns,clath_rxns);
      rxnNames = vertcat(rxnNames,clath_names);
      connector = 'XXX_mature[cv]';
      location = psim.location(p_row);
      rxns = vertcat(rxns,strcat(connector,{' => XXX_mature'},location));
      rxnNames = vertcat(rxnNames,strcat({'Final_location_'},location));
      rxns = vertcat(rxns,strcat({'XXX_mature'},location,{' => '}));
      rxnNames = vertcat(rxnNames,strcat(protein,{'_Final_demand'}));
     %------------------------------------------------------------------------------
      %Send to corresponding location
    else
      location = psim.location(p_row);
      rxns = vertcat(rxns,strcat(connector, {' => XXX-preSV[g]'}));
      rxnNames = vertcat(rxnNames,{'Start_Secretion'});
      sv_rxns = {};
      sv_names = {};
      [sv_rxns, sv_names] = addPathway('SV',sv_rxns,sv_names,tempRxns);
      for i = 1:numel(sv_rxns)
        sv_rxns(i) = strrep(sv_rxns(i),'!',num2str(clathrin_coeff));
      end
      rxns = vertcat(rxns,sv_rxns);
      rxnNames = vertcat(rxnNames,sv_names);
      if strcmp(psim.location(p_row),'') == 1
        location = '[s]';
      end
      rxns = vertcat(rxns,strcat({'XXX_mature[sv]'},{' => XXX_mature'},location));
      rxnNames = vertcat(rxnNames,strcat({'Final_location_'},location));
      rxns = vertcat(rxns,strcat({'XXX_mature'},location,{' => '}));
      rxnNames = vertcat(rxnNames,strcat(protein,{'_Final_demand'}));
    end
    %------------------------------------------------------------------------------
    %Add coeffcients to SP_degradation and Ubiquitin_degradation reactions (if applicable)
    if ismember({'SP_degradation'}, rxnNames) == 1
      SPaas = aacount(sequence(1:22)); %Amino acids in signal peptide assuming length is 22 on average
      rxns(find(ismember(rxnNames,'SP_degradation'))) = substitute_AAs_count(rxns(find(ismember(rxnNames,'SP_degradation'))),SPaas);
    end
    if ismember('Ubiquitination_degradation',rxnNames) == 1
      new_aas = aacount(sequence(1:22));
      rxns(find(ismember(rxnNames,'Ubiquitination_degradation'))) = strrep(rxns(find(ismember(rxnNames,'Ubiquitination_degradation'))),'!', '?');
      rxns(find(ismember(rxnNames,'Ubiquitination_degradation'))) = substitute_AAs_count(rxns(find(ismember(rxnNames,'Ubiquitination_degradation'))),new_aas);
      rxns(find(ismember(rxnNames,'Ubiquitination_degradation'))) = strrep(rxns(find(ismember(rxnNames,'Ubiquitination_degradation'))),'?',num2str(length-22));
    end
    %------------------------------------------------------------------------------
    GPRs = getGPRsFromRxnNames(rxnNames, tempRxns);
    Comps = getCompsFromRxnNames(rxnNames, tempRxns);
    %------------------------------------------------------------------------------
    %Change the reaction names and their formulas to include the UniProt Identifier
    for i = 1:numel(rxns)
      rxns(i) = insert_prot_name_in_rxnFormula((rxns(i)),protein);
    end
    for i = 1:numel(rxnNames)
      rxnNames(i) = insert_prot_name_in_rxnName(rxnNames(i),protein);
    end
    [rxns,rxnNames,GPRs, Comps] = addCanonicalRxns(rxns,rxnNames,GPRs, Comps,tempRxns);
  end
  %------------------------------------------------------------------------------
  rxnNames = strrep(rxnNames,'[','');
  rxnNames = strrep(rxnNames,']','');
  %------------------------------------------------------------------------------
  %add generated reactions to the model and generate secModel
  %------------------------------------------------------------------------------
  %first add list of genes which are not present in model and are in GPRs
  gprList = GPRs;
  gprList = strrep(gprList,'(','');
  gprList = strrep(gprList,')','');
  gprList = strrep(gprList,';','');
  gprList = strrep(gprList,' or ',' ');
  gprList = strrep(gprList,' and ',' ');
  for i = 1:numel(gprList)
    genes(1,:) = strsplit(char(gprList(i)));
    if i == 1
      geneList = genes(:);
    else
      geneList = vertcat(geneList,genes(:));
    end
    clear genes
  end
  geneList = unique(geneList);
  genesToAdd.genes = setdiff(geneList,secModel.genes);
  genesToAdd.genes = genesToAdd.genes(~cellfun('isempty', genesToAdd.genes));
  secModel = addGenesRaven(secModel,genesToAdd);
  %------------------------------------------------------------------------------
  %find and add new metabolites to secModel
  metList = parseRxnEqu(rxns);
  outModel = addCompartmentalizedMetNames(secModel);
  newMets = setdiff(metList,outModel.metNamesC); %list of new metabolites with their compartment.
  if isempty(newMets) == 0
    [newNames, newComps] = splitComp(newMets);
    metsToAdd.compartments = newComps;
    metsToAdd.metNames = newNames;
    metsToAdd.mets = strcat(metsToAdd.metNames,metsToAdd.compartments);
    secModel = addMets(secModel,metsToAdd);
  end
  %------------------------------------------------------------------------------
  %finding reactions which are not in the model
  rxnsToAdd.rxns = rxnNames(ismember(rxnNames,secModel.rxns) == 0);
  rxnsToAdd.equations = rxns(ismember(rxnNames,secModel.rxns) == 0);
  rxnsToAdd.rxnNames = rxnNames(ismember(rxnNames,secModel.rxns) == 0);
  rxnsToAdd.grRules = GPRs(ismember(rxnNames,secModel.rxns) == 0);
  rxnsToAdd.subSystems(1:numel(rxnsToAdd.rxnNames)) = {'protein secretion'};

  com = Comps(ismember(rxnNames,secModel.rxns) == 0);
  %Add protein secretion reactions to human1-drived model and modify grule Matrix
  secModel = addRxns(secModel,rxnsToAdd,3,[],false);
  %Add compartments for new reactions
  %secModel.rxnComps(numel(model.rxns)+1:numel(secModel.rxns)) = com(1:end);
  %------------------------------------------------------------------------------
  addedRxns.rxnsNames = rxnsToAdd.rxns;
  addedRxns.rxnFormula = rxnsToAdd.equations;
  addedRxns.rxnGPRs = rxnsToAdd.grRules;
  addedRxns.rxnComps = com;
  addedRxns.newMets = newMets;
  addedRxns.newGens = genesToAdd.genes;
end
