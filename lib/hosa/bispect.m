function [Bspec,waxis] = bispect (ma, ar, nfft)
%BISPECT  Theoretical bispectrum of an ARMA model
%	[Bspec,waxis] = bispect (ma, ar,  nfft)
%	ma   - ma parameter vector
%	ar   - ar parameter vector [default: [1]]
%	nfft - FFT length          [default: 512]
%	Bspec - Bispectrum of ARMA model, an nfft by nfft array
%	        with axes pointing down and right
%	waxis - frequency axis associated with the rows/columns of Bspec

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc.
%       $Revision: 1.7 $
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

Xf    = freqz(ma,ar,nfft,'whole') ;
Xfc   = conj(Xf);
Bspec = (Xf * Xf') .* hankel(Xfc,Xfc([nfft,1:nfft-1]));

Bspec = fftshift(Bspec);     % center the origin; normal axes

if (rem(nfft,2) == 0)
   waxis = [-nfft/2:nfft/2-1]'/nfft;
else
   waxis = [-(nfft-1)/2:(nfft-1)/2]'/nfft;
end

%contour(abs(Bspec),6,waxis,waxis),grid
contour(waxis,waxis,abs(Bspec),6), grid on
xlabel('f1'), ylabel('f2')

x = [1  1  2/3  1  1/3   0  -1  -1  -2/3  -1  -1/3  0] ;
y = [1  0 -1/3 -1 -2/3  -1  -1   0   1/3   1   2/3  1] ;
hold on
for k=1:12
   plot([0 x(k)], [0 y(k)], '--')
end
hold off
set(gcf,'Name','Hosa BISPECT')
