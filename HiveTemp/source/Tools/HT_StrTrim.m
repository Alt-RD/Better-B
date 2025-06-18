function _str = HT_StrTrim(_str, _chars)
  if nargin < 2,
    _chars = " \r\n\t";
  endif

  if iscell(_str)
    assert(all(cellfun(@(v) ischar(v), _str)));
    for i=1:numel(_str)
      _str{i} = BH1_strTrim(_str{i}, _chars);   % Recursive call
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

    if isempty(_str), return; endif

    ind = numel(_str)+1;
    while ~isempty(strchr(_str(ind-1), _chars))
      ind--;
      if ind <= 1
        break;
      endif
    endwhile

    if ind <= numel(_str), _str(ind:end) = []; endif;
  endif

##  assert(nargin == 2, 'Missing parameter');
##
##  if iscell(S)
##    for i=1:numel(S)
##      S{i} = HT_StrTrim(S{i}, charList);
##    endfor
##  else
##    while ~isempty(S) && any(S(1) == charList)
##      S(1) = [];
##    endwhile
##
##    while ~isempty(S) && any(S(end) == charList)
##      S(end) = [];
##    endwhile
##  endif
endfunction
