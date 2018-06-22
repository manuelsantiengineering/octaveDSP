% (file name: ustep.m)
% The unit step function
% Usage y = ustep(t);
% t must be real and can be a vector or a matrix
% Obtained from:
% Lathi,B.P. & Zhi Ding. (2009). Modern Digital and Analog Communication
% Systems. New York, NY. Oxford University Press.

function y=ustep(t)
    y = (t>=0);
end

