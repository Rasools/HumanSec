function newFormula=substitute_AAs_count(formula,AAcounts)
  %Define a function that substitutes question marks for AA counts in a formula (string)
  %Useful for translation and protein degradation pathway
  %
  %
  %
  %
  %
  %
  %
  %
  % Rasool Saghaleyni, 2018-02-16
  %replace "?" for aa based on aa count data
  newFormula = strrep(formula, '? glycine[c]', strcat(num2str(AAcounts.G),' glycine[c]'));
  newFormula = strrep(newFormula, '? alanine[c]', strcat(num2str(AAcounts.A),' alanine[c]'));
  newFormula = strrep(newFormula, '? valine[c]', strcat(num2str(AAcounts.V),' valine[c]'));
  newFormula = strrep(newFormula, '? leucine[c]', strcat(num2str(AAcounts.L),' leucine[c]'));
  newFormula = strrep(newFormula, '? isoleucine[c]', strcat(num2str(AAcounts.I),' isoleucine[c]'));
  newFormula = strrep(newFormula, '? methionine[c]', strcat(num2str(AAcounts.M),' methionine[c]'));
  newFormula = strrep(newFormula, '? tryptophan[c]', strcat(num2str(AAcounts.W),' tryptophan[c]'));
  newFormula = strrep(newFormula, '? phenylalanine[c]', strcat(num2str(AAcounts.F),' phenylalanine[c]'));
  newFormula = strrep(newFormula, '? proline[c]', strcat(num2str(AAcounts.P),' proline[c]'));
  newFormula = strrep(newFormula, '? serine[c]', strcat(num2str(AAcounts.S),' serine[c]'));
  newFormula = strrep(newFormula, '? threonine[c]', strcat(num2str(AAcounts.T),' threonine[c]'));
  newFormula = strrep(newFormula, '? cysteine[c]', strcat(num2str(AAcounts.C),' cysteine[c]'));
  newFormula = strrep(newFormula, '? tyrosine[c]', strcat(num2str(AAcounts.Y),' tyrosine[c]'));
  newFormula = strrep(newFormula, '? asparagine[c]', strcat(num2str(AAcounts.N),' asparagine[c]'));
  newFormula = strrep(newFormula, '? glutamine[c]', strcat(num2str(AAcounts.Q),' glutamine[c]'));
  newFormula = strrep(newFormula, '? glutamate[c]', strcat(num2str(AAcounts.E),' glutamate[c]'));
  newFormula = strrep(newFormula, '? aspartate[c]', strcat(num2str(AAcounts.D),' aspartate[c]'));
  newFormula = strrep(newFormula, '? lysine[c]', strcat(num2str(AAcounts.K),' lysine[c]'));
  newFormula = strrep(newFormula, '? arginine[c]', strcat(num2str(AAcounts.R),' arginine[c]'));
  newFormula = strrep(newFormula, '? histidine[c]', strcat(num2str(AAcounts.H),' histidine[c]'));
end
