function _str = HT_StrTrimLeft(_str, _chars)
  if nargin < 2,
    _chars = " \r\n\t";
  endif

  if iscell(_str)
    assert(all(cellfun(@(v) ischar(v), _str)));
    for i=1:numel(_str)
      _str{i} = BH1_StrTrimLeft(_str{i}, _chars);   % Recursive call
    endfor
  else
    assert(ischar(_str));
    if isempty(_str), return; endif

    ind = 0;
    while ~isempty(strchr(_str(ind+1), _chars))
      ind++;
      if ind >= numel(_str)
        break;
      endif
    endwhile

    if ind > 0, _str(1:ind) = []; endif
  endif
endfunction

