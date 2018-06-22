function [pout] = pnrz(T)

%
% USAGE: pout = pnrz(T)
%
% INPUTS:
%              T - Width of pulse (samples)
%
% OUTPUTS:
%             pout - Pulse of with T [samples].
%
% pnrz:  This routine generates a non-return-to-zero
%              pulse of T samples wide.
%

if      ~nargin; help pnrz; return; end
if      nargin < 1;     help pnrz; return; end

pout=ones(1,T);


return

