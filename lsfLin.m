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

% vector length
N = length( x );

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
  one_sig2_sum  = N;
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

  % Uncertainties. Calculated from error propagation
  sig2_a0_ep =  1 ./ Delta * x2_sig2_sum ;
  sig2_a1_ep =  1 ./ Delta * one_sig2_sum ; 

  % Uncertainties. Calculated using Gaussian properties
  aveXw  = x_sig2_sum / one_sig2_sum;
  aveX2w = x2_sig2_sum / one_sig2_sum;

  % Shifted x. Done to get uncorrelated errors
  xprime = x - aveXw;

  sig2_a1_g = 1 / sum( w .* xprime .^2 ) ;
  sig2_a0p_g = 1 / one_sig2_sum ;
  sig2_a0_g = sig2_a0p_g  +  aveXw .^ 2 * sig2_a1_g ;

end

% Shifted y intercept, should be equal to <y>
a0p = a0 + a1 .* x_sig2_sum / one_sig2_sum;

% Calculate errors assuming all errors are equal- scatter of points
sig2_guess = 1 / ( N - 2 ) .* sum( ( a0 + a1 .* x - y ) .^ 2 );
aveXw  = sum(x) / N;
xprime = x - aveXw ;

sig2_a1_sp = sig2_guess  ./ sum( xprime .^ 2 ) ;
sig2_a0p_sp = sig2_guess / N ;
sig2_a0_sp = sig2_a0p_sp + aveXw .^ 2 * sig2_a1_sp ;

% Put it in a struct

% a0
fitobj.a0 = a0;
if nargin == 3; 
  fitobj.sig_a0_ep = sqrt( sig2_a0_ep ); 
  fitobj.sig_a0_g = sqrt( sig2_a0_g ); 
end
  fitobj.sig_a0_sp = sqrt( sig2_a0_sp ); 

% a0' (prime)
fitobj.a0p = a0p;
if nargin == 3; 
  fitobj.sig_a0p_g = sqrt( sig2_a0p_g ); 
end
  fitobj.sig_a0p_sp = sqrt( sig2_a0p_sp ); 

% a1  
fitobj.a1 = a1;
if nargin == 3; 
  fitobj.sig_a1_ep = sqrt( sig2_a1_ep ); 
  fitobj.sig_a1_g = sqrt( sig2_a1_g ); 
end
  fitobj.sig_a1_sp = sqrt( sig2_a1_sp ); 

  fitobj.sig_Chi2  = sqrt(sig2_guess);
  


% Residual
gof.res  = sum( ( y - (a1 * x  + a0) )  .^ 2 ); 


