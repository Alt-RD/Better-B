% Résolution numérique d'un circuit électrique 1D.
% Chaque noeud est défini par par une capacité, une conductance avec le noeud suivant,
% un apport de flux de chaleur, une résistance thermique à une source externe constante
% Retourne une matrice T dim Nnoeud x Ntimes
function T = HT_NodalNetwork1D(t, Cvec, Gvec, Phivec, Hvec, Tvec, Tinit, Rvec, Trad)
  n = max([numel(Cvec), numel(Gvec), numel(Phivec), numel(Hvec), size(Tvec, 1)]);

  % If the only vector is Gvec, it's a special case since numel(Gvec) = n-1
  if n == numel(Gvec)
    n = n+1;
  else
    assert((numel(Gvec) == n-1) || (numel(Gvec)==1), sprintf("Size of vector Gvec (=%d) should be n-1. (n=%d)", numel(Gvec), n));
  endif

  TvecIsFunc = strcmp(class(Tvec), 'function_handle');
  PhivecIsFunc = strcmp(class(Phivec), 'function_handle');
  TradvecIsFunc = strcmp(class(Trad), 'function_handle');
  PhivecIsNonlinear = false;
  if PhivecIsFunc
    lPhivecInfos = functions(Phivec);
    [~, ~, ~, ~, t, ~, ~] = regexp(lPhivecInfos.function, '@\(([,\w]+)\)');
    lArgs = cellfun(@(v) strtrim(v), strsplit(t{1}{1}, ','), 'UniformOutput', false);
    if numel(lArgs) > 1
      assert(strcmp(lArgs{2}, 'T'), 'Invalid argument of phi function');
      PhivecIsNonlinear = true;
    endif
  endif

  if numel(Cvec) == 1, Cvec = repmat(Cvec, n,1); end;
  if numel(Gvec) == 1, Gvec = repmat(Gvec, n-1,1); end;
  if numel(Hvec) == 1, Hvec = repmat(Hvec, n,1); end;
  if numel(Tinit) == 1, Tinit = repmat(Tinit, n,1); end;
  if numel(Rvec) == 1, Rvec = repmat(Rvec, n,1); end;
  if numel(Trad) == 1, Trad = repmat(Trad, n,1); end;

  if ~PhivecIsFunc && numel(Phivec) == 1, Phivec = repmat(Phivec, n,1); end;
  if ~TvecIsFunc && numel(Tvec) == 1, Tvec = repmat(Tvec, n,1); end;

  if size(Cvec,2) > 1, Cvec = Cvec'; end;
  if size(Gvec,2) > 1, Gvec = Gvec'; end;
  if size(Phivec,2) > 1, Phivec = Phivec'; end;
  if size(Hvec,2) > 1, Hvec = Hvec'; end;
  if (size(Tvec,2) > 1) && (size(Tvec,1) == 1), Tvec = Tvec'; end;
  if size(Tinit,2) > 1, Tinit = Tinit'; end;
  if size(Rvec,2) > 1, Rvec = Rvec'; end;
  if size(Trad,2) > 1, Trad = Trad'; end;

  Amat = sparse(n,n);
  Amat(1,1) = -Gvec(1); ... - Hvec(1);
  Amat(1,2) = Gvec(1);
  Amat(n,n) = -Gvec(n-1); ... - Hvec(n);
  Amat(n,n-1) = Gvec(n-1);

  for i=2:(n-1)
    Amat(i,i) = -Gvec(i) - Gvec(i-1); ... - Hvec(i);
    Amat(i,i+1) = Gvec(i);
    Amat(i,i-1) = Gvec(i-1);
  endfor

  T = NaN(n,numel(t));
  T(:,1) = Tinit;

  Cmat = spdiags(Cvec, 0, n,n);
  Hmat = spdiags(Hvec, 0, n,n);
  Rmat = spdiags(Rvec, 0, n,n);

  if TvecIsFunc
    TvecOld = Tvec(1);
  elseif size(Tvec,2) > 1
    TvecOld = Tvec(:,1);
  else
    TvecOld = Tvec;
  endif

  if PhivecIsNonlinear
    Phivecloc = Phivec(i, TvecOld);
  elseif PhivecIsFunc
    PhivecOld = Phivec(i);
  elseif size(Phivec,2) > 1
    PhivecOld = Phivec(:,i);
  else
    PhivecOld = Phivec;
  endif

  if TradvecIsFunc
    TradOld = Trad(i);
  elseif size(Trad,2) > 1
    TradOld = Trad(:,i);
  else
    TradOld = Trad;
  endif

  OptOptions = optimset('AutoScaling', 'off', ...
                        'MaxIter', 500, ...
                        'MaxFunEvals', 800*n, ...
                        'Display', 'iter', ...
                        'TolX', 1E-6);

  for i=2:numel(t)
    dt = t(i) - t(i-1);   % Pas de temps courant

    if TvecIsFunc
      Tvecloc = Tvec(i);
    elseif size(Tvec,2) > 1
      Tvecloc = Tvec(:,i);
    else
      Tvecloc = Tvec;
    endif

    if PhivecIsNonlinear
      Phivecloc = Phivec(i, TvecOld);
    elseif PhivecIsFunc
      Phivecloc = Phivec(i);
    elseif size(Phivec,2) > 1
      Phivecloc = Phivec(:,i);
    else
      Phivecloc = Phivec;
    endif

    if TradvecIsFunc
      Tradloc = Trad(i);
    elseif size(Trad,2) > 1
      Tradloc = Trad(:,i);
    else
      Tradloc = Trad;
    endif

    if isempty(Trad)
      T(:,i) = linsolve(2*Cmat-Amat*dt+Hmat*dt, ...
                      (2*Cmat+Amat*dt-Hmat*dt)*T(:,i-1) + ...
                      Hmat * dt * (TvecOld + Tvecloc) + dt * (PhivecOld + Phivecloc));
    else
      R = Rmat * dt * 5.67E-8*273.15^4;

      A = 2*Cmat-Amat*dt+Hmat*dt;
      B = (2*Cmat+Amat*dt-Hmat*dt)*T(:,i-1) - R*(T(:,i-1)/273.15).^4 + Hmat * dt * (TvecOld + Tvecloc) + dt * (PhivecOld + Phivecloc);
      fobj = @(T) sqrt(sumsq(A*T+R*(T/273.15).^4 - B - R*((TradOld/273.15).^4 + (Tradloc/273.15).^4)));
      [T(:,i) fval info] = fminunc(fobj, T(:,i-1), OptOptions);
      if ~any(info == [1 2 3])
        disp(sprintf('Convergence failed: %d at iteration %d', info, i));
      endif
    endif

    TvecOld = Tvecloc;
    PhivecOld = Phivecloc;

    if mod(i-1,ceil((numel(t)-1)/10)) == 0
      disp(sprintf('Processing %d/%d', i-1, numel(t)-1));
    endif
  endfor

end
