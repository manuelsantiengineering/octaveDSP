function [bb, tim] = rcosfir(r, N_T, rate, T, fil_type, col)
%RCOSFIR Design a raised cosine FIR filter.
%   B = RCOSFIR(R, N_T, RATE, T) designs and returns a raised cosine FIR filter.
%   A raised cosine filter is typically used to shape and oversample a symbol
%   stream before modulation/transmission as well as after reception and
%   demodulation.  It is used to reduce the bandwidth of the oversampled symbol
%   stream without introducing intersymbol interference.
%
%   The time response of the raised cosine filter is,
%
%     h(t) = SINC(t/T) COS(pi R t/T)/(1 - 4 R^2 t^2 /T^2)
%
%   The frequency domain has the spectrum
%
%            /  T                                 when 0 < |f| < (1-r)/2/T
%            |          pi T         1-R    T          1-R         1+R
%     H(f) = < (1 + cos(----) (|f| - ----) ---    when --- < |f| < ---
%            |            r           2T    2          2 T         2 T
%            \  0                                 when |f| > (1+r)/2/T
%
%
%   T is the input signal sampling period, in seconds.  RATE is the
%   oversampling rate for the filter (or the number of output samples per input
%   sample). The rolloff factor, R, determines the width of the transition
%   band.  R has no units.  The transition band is (1-R)/(2*T) < |f| <
%   (1+R)/(2*T).
%
%   N_T is a scalar or a vector of length 2.  If N_T is specified as a
%   scalar, then the filter length is 2 * N_T + 1 input samples.  If N_T is
%   a vector, it specifies the extent of the filter.  In this case, the filter
%   length is N_T(2) - N_T(1) + 1 input samples (or
%   (N_T(2) - N_T(1))* RATE + 1 output samples).
%
%   The default value for N_T is 3.  The default value of RATE is 5.
%   The default value of T is 1.
%
%   B = RCOSFIR(R, N_T, RATE, T, FILTER_TYPE) designs and returns a
%   square root raised cosine filter if FILTER_TYPE == 'sqrt'. The default
%   value of FILTER_TYPE, 'normal', returns a normal raised cosine filter.
%
%   RCOSFIR(R, N_T, RATE, T, FILTER_TYPE, COL) produces the time response
%   and frequency response with the curve color as specified in the string
%   variable COL. The string in COL can be any type as defined in
%   PLOT.  If COL is not present, the default color will be used in the plot
%
%   [B, Sample_Time] = RCOSFIR(...) returns the FIR filter and the output sample
%   time for the filter. Note that the filter sample time is T / RATE.
%
%   See also RCOSIIR, RCOSFLT, RCOSINE, FIRRCOS, RCOSDEMO, GAUSSFIR.

%   Copyright 1996-2005 The MathWorks, Inc.
%   $Revision: 1.15.4.1 $   $Date: 2004/12/10 19:22:59 $

%routine check
if nargin < 1
    error('Not enough input variables for RCOSFIR')
elseif nargin < 2
    N_T = [3 3]; rate = 5; T = 1; fil_type = 'normal';
elseif nargin < 3,
    rate = 5; T = 1; fil_type = 'normal';
elseif nargin < 4,
    T = 1; fil_type = 'normal';
elseif nargin < 5,
    fil_type = 'normal';
end;

if (r < 0) | (r > 1) | ~isreal(r)
    error('The Rolloff factor in RCOSFIR must be a positive integer in the range, [0, 1].')
end;

[N_T, rate, T, fil_type] = checkinp(N_T, rate, T, fil_type,...
    [3 3], 5,  1, 'normal');
if length(N_T) < 2
    N_T = [N_T N_T];
end;

if (rate <= 1) | (ceil(rate) ~= rate)
    error('RATE in RCOSFIR must be an integer greater than 1')
end

% calculation

N_T(1) = -abs(N_T(1));
time_T = [0 : 1/rate : max(N_T(2), abs(N_T(1)))];
cal_time = time_T * T;
time_T_r = r * time_T;
if ~isempty(findstr(fil_type,'root')) | ~isempty(findstr(fil_type,'sqrt'))
    % square root raised cosine
    b=firrcos(rate*(N_T(2)-N_T(1)),1/(2*T),r,rate/T,'r','sqrt',-N_T(1)*rate)*sqrt(rate);
else
    % regular raised cosine
    b=firrcos(rate*(N_T(2)-N_T(1)),1/(2*T),r,rate/T,'r',[],-N_T(1)*rate)*rate;
end

tim = cal_time(2) - cal_time(1);

% In the case needs a plot
if nargout < 1
    if nargin < 6
        col = '';
    end;

    % the time response part
    hand = subplot(211);
    % dont filter, plot using plot([0 : 1/rate : N_T(2) - N_T(1)],b) insteat
    out = filter(b, 1, [1, zeros(1, length(cal_time) - 1)]);
    plot(cal_time, out, col)
    % if not hold, change the axes
    hol = get(hand,'NextPlot');
    if (hol(1:2) ~= 'ad') | (max(get(hand,'Ylim')) < max(b))
        axis([min(cal_time), max(cal_time), min(out) * 1.1, max(out) * 1.1]);
        xlabel('time');
        title('Impulse Response of the Raised Cosine Filter (with time shift)')
    end;

    % the frequency response part
    hand = subplot(212);
    len = length(b);
    P = abs(fft(b)) * abs(N_T(2) - N_T(1)) / len * T;
    f = (0 : len / 2) / len * rate / T;
    ind = find(f < 1.5 / T);
    f = f(ind);
    P = P(ind);
    plot(f, P, col);
    hol = get(hand, 'NextPlot');
    if hol(1:2) ~= 'ad'
        xlabel('frequency');
        ylabel('Amplitude');
        title('Frequency Response of the Raised Cosine Filter')
    end;
else
    bb = b;
end;
%--end of rcosfir.m--

return

%%%%%%

function [y1,y2,y3,y4,y5,y6,y7,y8,y9,y10] = checkinp(x1,...
    x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20)
%CHECKINP Checks input and fills in the default.
%       [y1, .., yN] = checkinp(x1, ..., xN, z1, ..., zN) checks if xi is
%       empty. If xi is empty, the function fills yi with the given default
%       value zi. Otherwise, the function assigns yi=xi. The number of
%       variable N is limited to less than 10.

%       Wes Wang 5/26/94, 9/30/95 
%       Copyright 1996-2002 The MathWorks, Inc.
%       $Revision: 1.9 $

if (nargin ~= 2* nargout)
    error('The number of input and output variables in CHECKINP does not match.');
end

for i = 1: nargout
    eval(['test_flag = x', num2str(i),';']);
    if isempty(test_flag)
        eval(['y',num2str(i),'=x',num2str(nargout+i),';']);
    else
        eval(['y',num2str(i),'=x',num2str(i),';']);
    end
end

return

%%%%%%
