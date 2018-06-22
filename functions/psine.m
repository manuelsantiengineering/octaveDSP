function [pout] = psine(T)

%
% USAGE: pout = psine(T)
%
% INPUTS:
%              T - Width of pulse (samples)
%
% OUTPUTS:
%             pout - Pulse with T/2 non-zero samples.
%
% psine:  This routine generates a sinusoidal
%              pulse of T samples wide.
%

if      ~nargin; help psine; return; end
if      nargin < 1;     help psine; return; end

pout=sin(pi*[0:T-1]/T);


return

