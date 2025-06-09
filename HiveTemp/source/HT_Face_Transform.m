%  This file is part of project HiveTemp.
%  This work was supported by the Better-B project, which has received funding
%  from the European Union, the Swiss State Secretariat for Education, Research
%  and Innovation (SERI) and UK Research and Innovation (UKRI) under the UK
%  government's Horizon Europe funding guarantee (grant number 10068544). Views
%  and opinions expressed are however those of the author(s) only and do not
%  necessarily reflect those of the European Union, European Research Executive
%  Agency (REA), SERI or UKRI. Neither the European Union nor the granting
%  authorities can be held responsible for them.
%
%  Copyright (c) 2022: Montpellier University
%  Copyright (c) 2023-2025: CoActions-AltRD-Emmanuel Ruffio
%  Author: emmanuel.ruffio@alt-rd.com
%
%  HiveTemp is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  HiveTemp is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with HiveTemp.  If not, see <https://www.gnu.org/licenses/>
% ========================================================================
function F = HT_Face_Transform(F, varargin)
  [lKeys, lValues, ~, lOptions] = HT_ReadVarargin(varargin);
  assert(HT_CheckType(F, "face"));

  for i=1:numel(lKeys)
    if strcmp(lKeys{i}, 'position')
      if isnumeric(lValues{i})
        assert(numel(lValues{i}) == 3);
        F = Int_SetAll(F, "globalPosition", lValues{i});
      else
        assert(isstruct(lValues{i}));
        lPositionInfos = lValues{i};
        lPositionInfos = HT_CheckField(lPositionInfos, "position", [], {'exist', @(v) (numel(v) == 3) && iscolumn(v)});
        lPositionInfos = HT_CheckField(lPositionInfos, "axis", []);
        lPositionInfos = HT_CheckField(lPositionInfos, "center", [0;0], {@(v) iscolumn(v)});

        for k=1:numel(F)
          lObj = Int_GetObject(F(i));
          lAxis = lPositionInfos.axis;

          if isempty(lAxis)
            lAxis = lObj.axis;
          else
            M = abs(lObj.axis' * lAxis);
            assert(sum(sum(abs(M-eye(2)))) < 1E-12 || sum(sum(abs(M-fliplr(eye(2))))) < 1E-12, 'Invalid axis');
          endif

          lCurrentPosition = lObj.globalPosition + lAxis * lPositionInfos.center;
          lOffset = lPositionInfos.position - lCurrentPosition;
          lOffset(isna(lOffset)) = 0;

          if iscell(F)
            F{i}.globalPosition += lOffset;
          else
            F(i).globalPosition += lOffset;
          endif

        endfor
      endif
    elseif strcmp(lKeys{i}, 'move')
      for i=1:numel(F)
        if iscell(F)
          F{i}.globalPosition += lValues{i};
        else
          F(i).globalPosition += lValues{i};
        endif
      endfor
    elseif strcmp(lKeys{i}, 'rotate')
      assert(isstruct(lValues{i}),'Invalid rotation information');
      for k=1:numel(F)
        lObj = Int_GetObject(F(i));

        lRotationInfos = HT_CheckField(lValues{i}, "center", lObj.globalPosition);
        lRotationInfos = HT_CheckField(lRotationInfos, "angle", 0);
        lRotationInfos = HT_CheckField(lRotationInfos, "axis", lObj.norm);
        if ischar(lRotationInfos.axis)
          if strcmp(lRotationInfos.axis, 'X')
            lRotationInfos.axis = [1; 0; 0];
          elseif strcmp(lRotationInfos.axis, 'Y')
            lRotationInfos.axis = [0; 1; 0];
          elseif strcmp(lRotationInfos.axis, 'Z')
            lRotationInfos.axis = [0; 0; 1];
          elseif strcmp(lRotationInfos.axis, 'U')
            lRotationInfos.axis = lObj.axis(:,1);
          elseif strcmp(lRotationInfos.axis, 'V')
            lRotationInfos.axis = lObj.axis(:,2);
          elseif strcmp(lRotationInfos.axis, 'W')
            lRotationInfos.axis = lObj;
          else
            error("Invalid axis");
          endif
        endif

        lVec = lObj.globalPosition - lRotationInfos.center;
        M = rotv(lRotationInfos.axis(:)', lRotationInfos.angle);

        if iscell(F)
          F{i}.axis = M * F{i}.axis;
          F{i}.norm = M * F{i}.norm;
          F{i}.globalPosition = lRotationInfos.center+M*lVec;
        else
          F(i).axis = M * F(i).axis;
          F(i).norm = M * F(i).norm;
          F(i).globalPosition = lRotationInfos.center+M*lVec;
        endif

      endfor
    elseif strcmp(lKeys{i}, 'scale')
      assert(isstruct(lValues{i}),'Invalid rotation information');
      for k=1:numel(F)
        lObj = Int_GetObject(F(i));

        lScaleInfos = lValues{i};
        lScaleInfos = HT_CheckField(lScaleInfos, "center", [0;0]);
        lScaleInfos = HT_CheckField(lScaleInfos, "factor", [1;1]);
        lScaleInfos = HT_CheckField(lScaleInfos, "direction", []);

        if (~isempty(lScaleInfos.direction))
          assert(numel(lScaleInfos.direction) == 3);
          assert(isempty(lScaleInfos.factor));
          lScaleInfos.factor = lObj.axis' * lScaleInfos.direction;
        endif

        if iscell(F)
          F{i}.size = lObj.size .* lScaleInfos.factor(:);
          F{i}.globalPosition -= lObj.axis * (lScaleInfos.center .* lObj.size .* (lObj.factor-1));
        else
          F(i).size = F(i).size .* lScaleInfos.factor(:);
          F(i).globalPosition -= lObj.axis * (lScaleInfos.center .* lObj.size .* (lObj.factor-1));
        endif

      endfor
    elseif strcmp(lKeys{i}, 'resize')
      if isnumeric(lValues{i})
        assert(numel(lValues{i}) == 2,'Invalid resize information');
        lSizeInfos = struct('size', lValues{i}, 'direction', [], 'center', [0.5;0.5]);
      else
        assert(isstruct(lValues{i}));
        lSizeInfos = lValues{i};
        lSizeInfos = HT_CheckField(lSizeInfos, "size", [0;0]);
        lSizeInfos = HT_CheckField(lSizeInfos, "direction", []);
        lSizeInfos = HT_CheckField(lSizeInfos, "center", [0.5;0.5]);
      endif

      assert(isscalar(lSizeInfos.size) || isempty(lSizeInfos.direction));

      for k=1:numel(F)
        lObj = Int_GetObject(F(i));
        lNewSize = lObj.size(:);

        if isempty(lSizeInfos.direction)
          lNewSize(~isna(lSizeInfos.size)) = lSizeInfos.size(~isna(lSizeInfos.size));
        elseif abs(dot(lSizeInfos.direction, lObj.axis(:,1)) - 1) < 1E-12
          lNewSize(1) = lSizeInfos.size;
        elseif abs(dot(lSizeInfos.direction, lObj.axis(:,2)) - 1) < 1E-12
          lNewSize(2) = lSizeInfos.size;
        else
          error('Invalid direction vector for resize operation');
        endif

        lNewPosition = lObj.globalPosition - lObj.axis * (lSizeInfos.center(:) .* (lNewSize - lObj.size));

        if iscell(F)
          F{i}.size = lNewSize;
          F{i}.globalPosition = lNewPosition;
        else
          F(i).size = lNewSize;
          F(i).globalPosition = lNewPosition;
        endif

      endfor
    else
      error(sprintf('Invalid property name <%s>', lKeys{i}));
    endif
  endfor
endfunction

function F = Int_SetAll(F, name, value)
  for i=1:numel(F)
    if iscell(F)
      F{i} = setfield(F{i}, name, value);
    else
      F(i) = setfield(F(i), name, value);
    endif
  endfor
endfunction

function F = Int_GetObject(F)
  if iscell(F)
    F = F{1};
  endif
endfunction
