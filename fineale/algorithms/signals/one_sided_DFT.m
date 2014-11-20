function [f,Y] =one_sided_DFT(dt,x)
% Compute the one-sided Discrete Fourier transform using FFT.
% 
% function Y =one_sided_DFT(dt,x)
% 
% x= this is the signal to transform
% dt= the time step
    NFFT = 2^nextpow2(length(x)); % the number of points in the FFT
    fs =(1/dt);% sampling frequency
    f = fs/2*linspace(0,1,NFFT/2);% Well-sampled frequencies are only up to Nyquist frequency
    Y = fft(x,NFFT)/NFFT;% Matlab FFT does not divide the Fourier coefficients by NFFT
    Y =2*abs(Y(1:NFFT/2));% one-sided FFT,  Factor 2: compensate for taking just one half
end