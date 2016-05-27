% Fit parameters
a       = 2;
b       = 3;

% "data"
N = 100; % Number of points
x   = linspace(-1,1,N)';
y_exact = a * x + b;
y_pert  = y_exact + 1 / 10 * randn( N, 1 );
%sig = 1/1000 * ones(N,1); 
%sig = ones(N,1); 
sig = 1/20 * sqrt( (1:N) )';

% errorbar( x, y_pert, sig,'o' );


% My fit routines
[fitlin, goodfit]  = lsfLin( x, y_pert, sig );
[fitpoly, ~]  = lsfPoly( x, y_pert, sig, 1 );
fitNoSig = lsfLin( x, y_pert );

% Polyfit
p        = polyfit( x, y_pert,1 );
poly.a0 = p(2);
poly.a1 = p(1);

% fit.m with weights
fo = fitoptions('Weights', 1 ./ sig.^2 );
ft = fittype('poly1');
[fitmat, gof, fitinfo]     = fit( x, y_pert, ft, fo);
% Uncertainites in fit.m 
ci = confint( fitmat, 0.67 );
matfit.Coeff = [fitmat.p2 fitmat.p1]; 
matfit.CoeffErr = [ ( ci(2,2) - ci(1,2) ) / 2  ( ci(2,1) - ci(1,1) ) / 2];

% fit.m with no weights
[fitmatNoSig, gof, fitinfoNs] = fit( x, y_pert,ft );
% Uncertainites
ci = confint( fitmatNoSig, 0.67 );
matfitNoSig.Coeff = [fitmatNoSig.p2 fitmatNoSig.p1]; 
matfitNoSig.CoeffErr = [ ( ci(2,2) - ci(1,2) ) / 2  ( ci(2,1) - ci(1,1) ) / 2];

% polyfitn
fitn = polyfitn(x, y_pert, 'constant, x' );

% Display Everything
% No error bar stuff
fprintf('fitNoSig: no sigma, my linear \n');
disp(fitNoSig);

fprintf('matfitNoSig: no sigma, uses fit\n');
disp(matfitNoSig);

fprintf('fitn: no sigma, uses polyfitn\n');
disp(fitn);

% Error bar stuff
%fprintf('fitlin: uses sigma, my linear \n'); 
%disp(fitlin);
%fprintf('fitpoly: uses sigma, my curvature matrix fit \n');
%disp(fitpoly);
%fprintf('poly: no sigma, uses polyfit\n');
%disp(poly);
%fprintf('fitmat: uses sigma, uses fit\n');
%disp(matfit);

% Playground

% build X
m = 2;
X = zeros(N,2);

yhat = fitNoSig.Coeff(1) + fitNoSig.Coeff(2) * x;
s    = norm(y_pert - yhat );

for i = 1:m
  X(:,i) = x .^ (i-1);
end
w = 1;
w = 1 ./ sig .^ 2;

alpha = X' * X;
[Q,R] = qr(X);

Rinv   =  eye(m) / R;
errMat = inv( alpha );

Var  = s^2 / (N-m) * sum(Rinv.^2,2);



