function newAbbreviation=insert_prot_name_in_rxnName(rxnAbbrev,protein)
  % Function that replaces the XXX template name for the Reaction abbreviation
  %
  %
  %
  %
  % Rasool Saghaleyni, 2018-02-16
  if contains(rxnAbbrev,protein) == 0
    newAbbreviation = strcat(protein,'_',rxnAbbrev);
  else
    newAbbreviation = rxnAbbrev;
  end
end
