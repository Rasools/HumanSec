function psim=importPsimFromExcel(fileName)
  %
  %   order of columns in PSIM matrix should be
  %    'Entry'	'Protein names'	'Length'	'Mass'	'SP'	'DSB'	'GPI'	'NG'	'OG'	'TMD'	'Location'	'Sequence'
  %   names could be different but order should be same.
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
  psim.entry = txt(2:end,1);
  psim.names = txt(2:end,2);
  psim.length = num(1:end,1);
  psim.mass = num(1:end,2);
  psim.sp = num(1:end,3);
  psim.dsb = num(1:end,4);
  psim.gpi = num(1:end,5);
  psim.ng = num(1:end,6);
  psim.og = num(1:end,7);
  psim.tmd = num(1:end,8);
  psim.location = txt(2:end,11);
  psim.sequence = txt(2:end,12);
  clear num txt
end
