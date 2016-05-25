% Fit parameters
a       = 2;
b       = 3;

% "data"
N = 100; % Number of points
x   = linspace(-1,1,N)';
y_exact = a * x + b;
y_pert  = y_exact + 1 / 10 * randn( N, 1 );
%sig = 1/1000 * ones(N,1); 
sig = ones(N,1); 
%sig = 1/20 * sqrt( (1:N) )';

% errorbar( x, y_pert, sig,'o' );

[fitlin, goodfit]  = lsfLin( x, y_pert, sig );
[fitpoly, ~]  = lsfPoly( x, y_pert, sig, 1 );
fitNoSig = lsfLin( x, y_pert );
p        = polyfit( x, y_pert,1 );

poly.a0 = p(2);
poly.a1 = p(1);

fo = fitoptions('Weights',sig);
ft = fittype('poly1');

fitmat = fit( x, y_pert, ft, fo);
fitn = polyfitn(x, y_pert, 'constant, x' );





