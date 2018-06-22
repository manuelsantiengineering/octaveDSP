function [y] = prcos(rollfac,len,T)

%
% USAGE: pout = prcos(T)
%
% INPUTS:
%         rollfac - Rolloff Factor (0 to 1).
%             len - Onesided pulse length (len=2T+1)
%               T - Oversampling Rate
%
% OUTPUTS:
%               y - Cosine pulse.
%
% prcos:  This routine generates a cosine pulse.
%

if      ~nargin; help prcos; return; end
if      nargin < 3;     help prcos; return; end

y=rcosfir(rollfac,len,T,1,'normal');

return

