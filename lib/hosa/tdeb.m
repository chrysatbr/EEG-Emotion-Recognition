function [delay,ctau] = tdeb (x,y, max_delay,  nfft, wind, nsamp, overlap)
%TDEB	Time Delay Estimation using conventional bispectrum method.
%	[delay,ctau] = tdeb (x,y,  max_delay, nfft, wind, segsamp, overlap)
%	x         - signal at sensor 1
%	y         - signal at sensor 2
%	max_delay - maximum expected delay in samples
%	nfft - fft length [default = power of two >= segsamp]
%	wind - window specification for frequency-domain smoothing
%	       if 'wind' is a scalar, it specifies the length of the side
%	          of the square for the Rao-Gabr optimal window  [default=5]
%	       if 'wind' is a vector, a 2D window will be calculated via
%	          w2(i,j) = wind(i) * wind(j) * wind(i+j)
%	       if 'wind' is a matrix, it specifies the 2-D filter directly
%	segsamp - samples per segment [default: 4*max_delay+1]
%	        - if x is a matrix, segsamp is set to the number of rows
%	overlap - percentage overlap, allowed range [0,99]. [default = 50];
%	        - if x is a matrix, overlap is set to 0.
%	delay     - estimated delay (in samples)
%       ctau      - estimated third-order hologram.

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.9 $
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

% parameter checks ---------------------------------------------------------

[lx,nrecs] = size(x);
[ly,ny] = size(y);
if (ly ~= lx | ny ~= nrecs)
    error('matrices x and y should have the same dimensions')
end
if (lx == 1), lx = nrecs; nrecs = 1; x = x(:); y = y(:); end

if (exist('max_delay') ~= 1)
   error('max_delay must be specified');
end

if (exist('nsamp') ~= 1) nsamp =  2^nextpow2(4*max_delay+1); end
if (exist('overlap') ~= 1) overlap = 50; end
overlap = min(99,max(overlap,0));

if (exist('wind') ~= 1)       wind = 5; end

if (nrecs > 1)
   overlap = 0;      nsamp   = lx;
end
if (exist('nfft') ~= 1)  nfft  =   0; end
if (nfft < nsamp)
   nfft = 2 ^ nextpow2(nsamp) ;
end

plotflag = 0;
Bxyx = bispecdx (x, y, x, nfft, wind, nsamp, overlap, plotflag) ;
Bxxx = bispecdx (x, x, x, nfft, wind, nsamp, overlap, plotflag) ;

Bxyx = fftshift(Bxyx);     % odd nfft ?
Bxxx = fftshift(Bxxx);

 coher = sum(Bxyx ./ Bxxx);     % sum integrates over omega_2
 ctau  = fftshift(real(ifft(coher)));
 [val, delay] = max(abs(ctau));
 delay = delay - nfft/2 - 1;
 disp(['Delay estimated by TDEB is ',int2str(delay)])
 plot(-nfft/2:nfft/2-1,ctau, -nfft/2:nfft/2-1,ctau,'o'), grid on 
 title(['TDEB: Hologram delay = ',int2str(delay)])
 set(gcf,'Name','Hosa TDEB')
