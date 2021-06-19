function [Tspec,waxis] = trispect (ma, ar, nfft, f3)
%TRISPECT  Theoretical trispectrum of an ARMA model: 2-D slice
%	[Tspec,waxis] = trispect (ma, ar,  nfft, f3)
%	ma   - ma parameter vector
%	ar   - ar parameter vector [default: [1]]
%	nfft - FFT length          [default: 512]
%       f3   - fixed frequency of third argument; default value is 0.
%		nominal range is [-0.5,0.5]
%	Tspec - Trispectrum of ARMA model, an nfft by nfft array
%	        with axes pointing down and right
%		Contains the trispectrum S3(f1,f2,f3), where
%			f3 is fixed.
%	waxis - frequency axis associated with the rows/columns of Tspec

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

% --------------

if (exist('ma') ~= 1)
   error('insufficient number of parameters')
end
q = length(ma)-1;
if (exist('ar') ~= 1) ar = 1; end
p = length(ar)-1;
pq = max(p,q) + 1;
if (exist('nfft') ~= 1) nfft = max(2^nextpow2(pq), 512); end
nfft  = max(nfft, 2^nextpow2(pq) );

if (exist('f3') ~= 1) f3 = 0; end
f3 = rem(f3,1);  f3 = f3 + (f3 < -1/2) - (f3 > 1/2);
ind = rem (fix(f3*nfft), nfft);
ind = ind + nfft * (ind < 0);

% ---------------

Xf    = freqz(ma,ar,nfft,'whole') ;
Xfc   = conj(Xf);
colind = rem( [1:nfft] + ind-1,      nfft) + 1;
rowind = rem( [nfft,1:nfft-1]+ind-1, nfft) + 1;
Tspec = Xf(ind+1) * (Xf * Xf') .* hankel(Xfc(colind),Xfc(rowind));

Tspec = fftshift(Tspec);     % center the origin; normal axes

if (rem(nfft,2) == 0)
   waxis = [-nfft/2:nfft/2-1]'/nfft;
else
   waxis = [-(nfft-1)/2:(nfft-1)/2]'/nfft;
end

%contour(abs(Tspec),6,waxis,waxis),grid
contour(waxis,waxis,abs(Tspec),6), grid on 
title(['Sliced Trispectrum, f3 = ',num2str(f3),' Hz'])

% It is fun to create a movie as follows:
% ma = [1 -1];  ar = [1 -0.8 0.65];   nfft = 64;
%  n = 10;  M = moviein(2*n+1);
%  for k=-n:n
%     trispect(ma,ar,nfft,k/(2*n));
%     M(:,k+n+1) = getframe;
%  end
%  movie(M)
%  clear M
% compare the f3=0 slice with bispect output.
