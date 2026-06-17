clear variables;

% Specify the main directory of HiveTemp library (relative or absolute)
HT_VAR_LIB_PATH = getenv('HTDIR');
addpath(make_absolute_filename(HT_VAR_LIB_PATH));

% Init the library
HT_Init();

lOperations = struct( 'name', 'Substract ambiant air temperature', ...
                      'outputTitle', @(k, infos) sprintf('Diff_%s', infos.sourceName1{k}), ... 'Diff_$(sourceName1)', ...
                      'description', '', ...
                      'function', @(i, D1, D2) D1-D2,... %[D1(:,1)-D1(:,2), D1(:,1)+D1(:,2)], ...
                      'options1', {'reject'}, ...
                      'options2', '', ...
                      'source1', {{'Timestamp', 'Index', 'Epoch'}}, ...
                      'source2', 'AirExt.T(T112@3.4f)');


[D T H] = HT_CsvOperation({'TestCsvOperation.csv'}, {}, lOperations, 'timeColumn1', 'Timestamp', 'timeColumn2', 'Timestamp');

%HT_CsvWrite('TestCsvOperation_OutFile.csv', D, T);




