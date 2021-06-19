function [delay,rxy] = tder(s1,s2,max_delay,segsamp,overlap,nfft)
%TDER	Time-Delay Estimation using ML windowed cross-correlation
%	[delay,rxy] = tder(s1,s2,max_delay,segsamp,overlap,nfft)
%	s1       - data at sensor 1
%	s2       - data at sensor 2
%	           s1 and s2 should have the same dimensions
%	max_delay  - maximum delay
%	segsamp - if not specified, set equal to the power of 2 just
%	          greater than 4*max_delay + 1;
%	        - if s1 is a matrix, segsamp is set to the number of rows
%	overlap - percentage overlap, allowed range [0,99]. [default = 50];
%	        - if s1 is a matrix, overlap is set to 0.
%	nfft - FFT length [default = power of two > segsamp]
%	delay      - estimated delay (positive means that y lags x)
%	rxy        - estimate of windowed cross-correlation function.

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.6 $
%  A. Swami   January 20, 1995

%     RESTRICTED RIGHTS LEGEND
% Use, duplication, or disclosure by the Government is subject to
% restrictions as set forth in subparagraph (c) (1) (ii) of the
% Rights in Technical Data and Computer Software clause of DFARS
% 252.227-7013.
% Manufacturer: United Signals & Systems, Inc., P.O. Box 2374,
% Culver City, California 90231.
%
%  This material may be reproduced by or for the U.S. Government pursuant
%  to the copyright license under the clause at DFARS 252.227-7013.

% --------- parameter checks ---------

[m1,n1] = size(s1);
[m2,n2] = size(s2);
if (m1 ~= m2 | n1 ~= n2)
   error('TDER: s1 and s2 should have identical dimensions');
end
if (exist('max_delay') ~= 1)
   error('TDER: max_delay must be specified ');
end

if (m1 == 1)   s1 = s1.' ; s2 = s2.' ;  end
[lx,nrecs] = size(s1);

lfft = 2^(nextpow2(4*max_delay+1));

if (exist('segsamp') ~= 1) segsamp = lfft; end
if (exist('overlap') ~= 1) overlap = 50;   end
overlap = max(0,min(overlap,99));

if (nrecs > 1)   segsamp = lx; overlap = 0; end
noverlap = fix(overlap/100 * segsamp);
wind     = hanning(segsamp);
if (exist('nfft') ~= 1) nfft = 0 ; end
if (nfft < segsamp)
   nfft = 2 ^ nextpow2(segsamp) ;
end

%------------ compute auto- and power-spectra ------
% note: spectrum computes lots of things that we do not really need.

P = spectrum(s1,s2,nfft,noverlap,wind);
pxy = P(:,3);                            % cross-spectrum
Cxy = P(:,5);                            % coherence
n = length(pxy);
Rw = Cxy ./ ( (1 - Cxy) .* abs(pxy) );   % ML-window
Rw = Rw .* pxy;                          % windowed cross-spectrum
Rw(1) = 0;                               % undo what spectrum does
% Rw = [Rw; 0; conj(Rw(n:-1:2))];
Rw = [Rw; conj(Rw(n-1:-1:2))];
rxy = real(fftshift(ifft(Rw)));          % windowed cross-correlation


[val,d] = max(abs(rxy));                 % locate peak
if (d > 1 & d < length(rxy))             % 3-point interpolation
  delay   = (2*d-1)/2 - (rxy(d)-rxy(d-1)) / (rxy(d+1)-2*rxy(d)+rxy(d-1));
else
  delay   = d;
end
delay   = delay - n ;                  % take care of offset

disp(['Estimated delay= ',num2str(delay)])

plot(-n+1:n-2,rxy, -n+1:n-2,rxy,'o'), grid on 
title(['TDER: Windowed Rxy: delay = ',num2str(delay)])
set(gcf,'Name','Hosa TDER')
return
