% [ fitobj, gof ] = lsfPoly( x, y, sig, deg )
%
% Decription: LSF for a polynomial of degree n, y = a_m x^m + ...a_1 x + a0
% Returns an object containing fit parameters, uncertainties, chi square. 
%
% fitobj: parameters and uncertainties if given sig
% gof: chi square and residual 
% 
% Reference: Data Reduction and Error Analysis, Bevington and Robinson
% Author: Michael Stefferson

function [fitobj, gof] = lsfPoly(x,y,sig,deg)

% Weights
if nargin == 3;
  [xr, xc] = size( x );
  w =  ones( xr, xc );
else
  w = 1 ./ sig .^ 2;
end

% Size of mat
N = length(x);

% Check that all vectors are the same length
[xr, xc] = size( x );
[yr, yc] = size( y );
[sigr, sigc] = size( sig );
if xr ~= yr || xc ~= yc || xr ~= sigr || xc ~= sigc
  error( 'vectors are not the same size. check row vs column vec');
end

% Set up the matrix beta = curvMat * coeff; 
% f_i = x^(i-1)
% beta_i = sum ( y .* f_i(x) / sig^2 )
beta    = zeros( deg + 1, 1);
coeff   = zeros( deg + 1, 1);
curvMat = zeros( deg + 1, deg + 1);
f       = zeros( N, deg + 1 );

% Build f, beta
for i = 1:deg+1
  f(:,i) = x .^ (i - 1);
  beta(i) =  sum( w .* y .* f(:,i) );
end

% Build curvature matrxi
for i = 1:deg+1
  for j = 1:deg+1
    curvMat(i,j) = sum( w .* f(:,i) .* f(:,j) );
  end
end
  
% Error matrix, inverse of a
% Inverting is slow, but we need diagonal of inverse for uncertainties
err = inv( curvMat );
coeff = err * beta;

sig_coeff = sqrt( diag(err) );

fitobj.coeff     = coeff';
fitobj.sig_coeff = sig_coeff';
fitobj.covarmat  = err;

gof = 0;
