

% This example shows how to numerically evaluate the exponential
% Fourier series coefficients Dn directly.
% Obtained from Modern Digital and Analog Communication Systems
% Fourth Edition B.P. Lathi & Zhi Ding
% The user needs to define a symbolic function g(t).

echo off; clear; clf;

    j=sqrt(-1); %Defines j for complex algebra.
    b =2; a=-1; %Determines the signal period
    to1=1.e-5; % Set the integration error tolerance
    T=b-a; %Length of the period
    N=11; %Number of Fourier Series coefficients on each side of zero freq.
    
    Fi=[-N:N]*2*pi/T; %Sets the frequency range
    
%Now calculate D-0 and store it in D(N+1);
    Func = @(t) funct_tri(t/2);
    D(N+1)=1/T*quad(Func,a,b,to1);
    
    for i=1:N
        Func= @(t) exp(-j*2*pi*t*i/T).*funct_tri(t/2);
        D(i+N+1)=1/T*quad(Func,a,b,to1);
        Func= @(t) exp(j*2*pi*t*(N+1-i)/T).*funct_tri(t/2);
        D(i)=1/T*quad(Func,a,b,to1);
    end
    
    figure(1);
    subplot(211);s1=stem([-N:N], abs(D));
    subplot(212);s2=stem([-N:N], angle(D));

        
        
        