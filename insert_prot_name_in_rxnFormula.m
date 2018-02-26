function newFormula=insert_prot_name_in_rxnFormula(formula,protein)
  %Function that replaces the XXX template name for the UniProt ID
  %
  %
  %
  %
  % Rasool Saghaleyni, 2018-02-16
  newFormula = strrep(formula,'XXX',protein);
end
