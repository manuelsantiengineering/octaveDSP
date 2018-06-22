% The rectangular function 
% Usage y = rect(t);
% t must be real and can be a vector or a matrix
% Obtained from:
% Lathi,B.P. & Zhi Ding. (2009). Modern Digital and Analog Communication
% Systems. New York, NY. Oxford University Press.

function y=rect(t)
    y = (sign(t+0.5) - sign(t-0.5) > 0);
end

