function [_ind, _str] = HT_CellStrFilter(_str, _wildcard)
  assert(nargin == 2, 'Invalid parameter');

  _wildcard = strrep(_wildcard, ')', '\)');
  _wildcard = strrep(_wildcard, '(', '\(');
  _wildcard = strcat('^', regexptranslate("wildcard", _wildcard), '$');

  s = regexp(_str, _wildcard);
  s = cellfun(@(v) ~isempty(v), s);
  _ind = find(s);

  if isargout(2)
    _str = _str(s);
  endif
endfunction
