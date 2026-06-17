clear variables;

% Specify the main directory of HiveTemp library (relative or absolute)
HT_VAR_LIB_PATH = getenv('HTDIR');
addpath(make_absolute_filename(HT_VAR_LIB_PATH));

% Init the library
HT_Init();

lFilePath = 'C:\Users\JB_Ch\Documents\Collectif\ADA\AURA\2025\Mesures\Synthese\SN10\SN10.csv'


  lParameters = struct( 'cacheFile', './../SN10.mat', ...
                        'maxHeaderLine', 100, ...
                        'comment', '#', ...
                        'delimiter', ',', ...
                        'options', '', ...
                        'columns', [], ...
                        'forceReload', false, ...
                        'strict', true, ...   % Abort if error occur
                        'fillColumn', true, ... % If strict is false, fill the missing column with <fillValue>
                        'fillValue', NA, ...
                        'bufferLine', 1000, ...
                        'timeColumn', [], ... Specify the time column (necessary of subsampling if the format is not epoch)
                        'timeFormat', 'yyyy/mm/ddTHH:MM:SS.FFF', ...
                        'subSampling', 0);

[D T] = HT_ReadCsvFile(lFilePath, lParameters);

##tind = find(strcmpi(T, 'Epoch'));
tind = find(strcmpi(T, 'TimeStamp'));
t = D{tind};

##[tCell tind tepoch] = HT_SplitTimeVector(t, 'epochType', '1970');
##[tinfos tepoch] = HT_SplitTimeVector(t, 'position', 'month', 'boundaryPolicy', 'overlap', 'emptyBin', 'reject');

lColumnList = {'NI23.Top.T(A47;2.3d)', 'NI23.Bot.T(A65;2.3f)'};
[S binInfos] = HT_CsvBinAverage(D, T, 'bin', 'month', ...
                       'unit', 'day', ...
                       'maxHoleTime', 10, ...
                       'timeColumn', 'Epoch', ...
                       'timeType', '1970', ...
                       'columns', lColumnList,
                       'reject', [],...
                       'exportRawData', true);

%HT_CsvWrite('TestCsvOperation_OutFile.csv', D, T);
lColors = zeros(columns(S), 3);
lLegendList = cell(rows(S), 1);

figure(1);
clf;
for i=1:rows(S) % Iterate over bins
  subplot(1,rows(S),i);

  for j=1:columns(S) % Iterate over sensors
    h = plot(S(i,j).t, S(i,j).data, 'displayname', S(i,j).name, ...
                                    'linewidth', 2); hold on;

    if i==1
      lColors(j,:) = get(h, 'color');
    endif

    set(h, 'color', lColors(j,:));

    lLegendList{i} = [lLegendList{i}; h];

    plot(S(i,j).t, S(i,j).rawData, '-', 'color', 1-(1-lColors(j,:))*0.4);
  endfor


##  hl = legend();
##  set(hl, 'location', 'south');
##  p = get(hl, 'dataaspectratio')
##  p(1) *= 0.8;
##  set(hl, 'dataaspectratio', p);
endfor

for i=1:rows(S) % Iterate over sensors
  subplot(1,rows(S),i);
  grid on;
  title(binInfos(i).name);
  xlabel('Temps (hour)');
  xtick = get(gca, 'xtick');
  xticklabel = arrayfun(@(v) sprintf('%.1f', v), xtick, 'UniformOutput', false);
  set(gca, 'xticklabel', xticklabel);

  hl = legend(lLegendList{i});
##  set(hl, 'location', 'south');
##  r(1) *= 0.7;
  set(hl, 'dataaspectratio', [10 25 1], 'location', 'south');
endfor


