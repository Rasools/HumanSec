function translationFormula=translationRxn(psim, protein)
  % Function that generates the translation reaction of a protein given its UniProt ID
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
  %Obtain protein sequence and then length
  p_row = find(ismember(psim.entry,protein)); %protein row in psim
  sequence = char(psim.sequence(p_row));
  AAcounts = aacount(sequence);
  templateFormula = string('? H2O[c] + ? ATP[c] + ? GTP[c] + ? glycine[c] + ? alanine[c] + ? valine[c] + ? leucine[c] + ? isoleucine[c] + ? methionine[c] + ? tryptophan[c] + ? phenylalanine[c] + ? proline[c] + ? serine[c] + ? threonine[c] + ? cysteine[c] + ? tyrosine[c] + ? asparagine[c] + ? glutamine[c] + ? glutamate[c] + ? aspartate[c] + ? lysine[c] + ? arginine[c] + ? histidine[c] => ? H+[c] + ? AMP[c] + ADP[c] + ? Pi[c] + ? GDP[c] + ? PPi[c] + XXX[c]');
  N = numel(sequence); %Number to replace ATP, GTP, PPi, Pi, H2O, AMP, and GTP
  %then replace stochiometric coefficients based on calculated N index
  translationFormula = substitute_AAs_count(templateFormula, AAcounts);
  translationFormula = strrep(translationFormula, '? H2O[c]', strcat(num2str(2*N-1),' H2O[c]'));
  translationFormula = strrep(translationFormula, '? H+[c]', strcat(num2str(2*N-1),' H+[c]'));
  translationFormula = strrep(translationFormula, '? PPi[c]', strcat(num2str(N),' PPi[c]'));
  translationFormula = strrep(translationFormula, '? Pi[c]', strcat(num2str(2*N-1),' Pi[c]'));
  translationFormula = strrep(translationFormula, '? GDP[c]', strcat(num2str(2*N-2),' GDP[c]'));
  translationFormula = strrep(translationFormula, '? AMP[c]', strcat(num2str(N),' AMP[c]'));
  translationFormula = strrep(translationFormula, '? GTP[c]', strcat(num2str(2*N-2),' GTP[c]'));
  translationFormula = strrep(translationFormula, '? ATP[c]', strcat(num2str(N+1),' ATP[c]'));
  translationFormula = strrep(translationFormula, 'XXX[c]', strcat(protein,'[c]'));
end
