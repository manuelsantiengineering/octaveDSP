% The triangle function
% Usage y = triangl(t);
% t must be real and can be a vector or a matrix
% Obtained from:
% Lathi,B.P. & Zhi Ding. (2009). Modern Digital and Analog Communication
% Systems. New York, NY. Oxford University Press.

function y=triangl(t)
    y = (1-abs(t)).*(t>=-1).*(t<1);
end

