clear variables;

% Specify the main directory of HiveTemp library (relative or absolute)
HT_VAR_LIB_PATH = getenv('HTDIR');
addpath(make_absolute_filename(HT_VAR_LIB_PATH));

% Init the library
HT_Init();

x = (0:2:10)';
y = 0.5*(0:2:10)';
y(3) = NA;
##y = ones(size(x));

##t = [1.2 3 4]';
t = [0.5 1 3 4]';

x = [0.5 1 1.5  2.5 3 4 ]';
y = [1   1 NA   2   2 3 ]';

##z = [0; (y(2:end) + y(1:end-1))/2 .* diff(x) ]
##z(isna(z)) = 0;
##cumsum(z)

##u = x;
##u(isna(y)) = NA;
##u = [0; diff(u)]
##u(isna(u)) = 0;
##cumsum(u)

##x = [0 1 1  2 2 3 ]';
##y = [1 1 -2 2 2 2 ]';

##x = [0 1 2 3 ]';
##y = [1 1 2 2 ]';
##
##cumsum((y(2:end) + y(1:end-1))/2 .* diff(x))


figure(1);
clf;
plot(x,y);
ylim([-2 3]);
grid on;

[u w]= HT_Resample(x,y,t, 'average', true, 'omitNaN', false)
##u = u ./ diff(t);

figure(1);
clf;
plot(x,y);
grid on;
