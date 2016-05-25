% [fitobj, gof] = lsfLin(x,y,sig)
%
% Description: LSF for a straight line, y = a1 x + a.
% Returns an object containing fit parameters, uncertainties, chi square. 
%
% fitobj: parameters and uncertainties if given sig
% gof: chi square and residual 
% 
% Reference: Data Reduction and Error Analysis, Bevington and Robinson
% Author: Michael Stefferson

function [fitobj, gof] = lsfLin(x,y,sig)

% if no sig, we assume it's the same for all points
if nargin == 2
  % Check that all vectors are the same length
  [xr, xc] = size( x );
  [yr, yc] = size( y );
  if xr ~= yr || xc ~= yc 
    error( 'vectors are not the same size. check row vs column vec');
  end
  
  % Calculate sums
  x2_sig2_sum = sum(  x .^ 2 ); 
  y_sig2_sum  = sum( y );
  x_sig2_sum  = sum( x );
  xy_sig2_sum = sum(  x .* y );
  one_sig2_sum  = length(x);
else
  % Check that all vectors are the same length
  [xr, xc] = size( x );
  [yr, yc] = size( y );
  [sigr, sigc] = size( sig );
  if xr ~= yr || xc ~= yc || xr ~= sigr || xc ~= sigc
    error( 'vectors are not the same size. check row vs column vec');
  end

  % Weights
  w = 1 ./ sig .^ 2;
  
  % Calculate sums
  x2_sig2_sum = sum( w .* ( x .^ 2  )  ); 
  y_sig2_sum  = sum( w .* y );
  x_sig2_sum  = sum( w .* x );
  xy_sig2_sum = sum( w .* ( x .* y ) );
  one_sig2_sum  = sum( w ); 
end

% Calculate fit parameters. Method of max likelihood. Assuming we are sampling
% a Guassian. Find parameters that maximize the probabiliy of getting those
% parameters, ie., minimizes chi squared.
Delta = one_sig2_sum * x2_sig2_sum - (x_sig2_sum) .^ 2;
a0 = 1 / Delta * ( x2_sig2_sum * y_sig2_sum - x_sig2_sum * xy_sig2_sum );
a1 = 1 / Delta * ( one_sig2_sum * xy_sig2_sum - x_sig2_sum * y_sig2_sum );

% Uncertainties/Chi Squared. Only defined if we have sigma
if nargin == 3
  % Chi square
  gof.chiSq = sum( ( ( y - (a1 * x  + a0) ) ./ sig ) .^ 2 );
  gof.chiSq_red = gof.chiSq ./ ( length(x) - 2 ); % Reduced chiSq
  % Uncertainties. Calucated from error propagation
  %sig_a0 = sqrt( 1 ./ Delta * x2_sig2_sum );
  %sig_a1 = sqrt( 1 ./ Delta * one_sig2_sum ); 
  xprime = x - mean(x);
  sig_a1 = sqrt( 1 ./ sum( xprime .^ 2 ) ); 
  sig_a0 = sqrt( one_sig2_sum ); 

end

gof.res  = sum( ( y - (a1 * x  + a0) )  .^ 2 ); 

% Put it in a struct
fitobj.a0 = a0;
if nargin == 3; fitobj.sig_a0 = sig_a0; end
fitobj.a1 = a1;
if nargin == 3; fitobj.sig_a1 = sig_a1; end

