function [mets, comps]=splitComp(metComps)
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
  % Rasool Saghaleyni, 2020-03-03
  for i=1:numel(metComps)
    a1 = char(metComps(i));
    if isempty(strfind(a1,'pm')) == 1 && isempty(strfind(a1,'cv')) == 1 && isempty(strfind(a1,'sv')) == 1 && isempty(strfind(a1,'gm')) && isempty(strfind(a1,'rm'))
      a2 = a1(1:end-3);
      a3 = a1(end-2:end);
    else
      a2 = a1(1:end-4);
      a3 = a1(end-3:end);
    end
    mets(i,1) = cellstr(a2);
    comps(i,1) = cellstr(a3);
  end
  comps = strrep(comps,'[','');
  comps = strrep(comps,']','');
end
