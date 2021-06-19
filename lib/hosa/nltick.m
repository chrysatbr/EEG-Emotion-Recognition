function [ht,qt] = nltick(x,y,nfft,wind,segsamp,overlap)
%NLTICK Second-order Volterra System Identification, for Gaussian inputs
%	[h, q] = nltick(x,y,nfft,wind,segsamp,overlap)
%	x  - input to the Volterra system
%	y  - output of the Volterra system
%	x and y must have identical dimensions;  multiple realizations
%	     are required, with columns corresponding to realizations.
%	nfft - FFT length to use for computing power spectra/bispectra
%	wind - window specification for frequency-domain smoothing
%	       if 'wind' is a scalar, it specifies the length of the side
%	          of the square for the Rao-Gabr optimal window  [default=5]
%	       if 'wind' is a vector, a 2D window will be calculated via
%	          w2(i,j) = wind(i) * wind(j) * wind(i+j)
%	       if 'wind' is a matrix, it specifies the 2-D filter directly
%	segsamp - samples per segment [default: so as to have 8 records]
%	        - if x is a matrix, segsamp is set to the number of rows
%	overlap - percentage overlap, allowed range [0,99]. [default = 50];
%	        - if x is a matrix, overlap is set to 0.
%	h   - estimated IR of the linear part
%	q   - estimated IR of the quadratic part

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.8 $
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

% Parameter checking ------------------------------------------

if (exist('x') ~= 1  | exist('y') ~= 1)
   error('both x and y must be specified')
end
if (size(x) ~= size(y))
   error('x and y must have the same dimensions')
end

[lx,nrecs] = size(x);
if (lx == 1) lx = nrecs; nrecs = 1; x = x(:); y = y(:); end

if (exist('overlap') ~= 1) overlap = 50; end
overlap = min(99,max(overlap,0));
if (exist('wind') ~= 1)       wind = 5; end
if (exist('segsamp') ~= 1)  segsamp = 0; end
if (nrecs == 1 & segsamp <= 0)
    segsamp = fix(lx/ (8 - 7 * overlap/100));
end
if (nrecs > 1) overlap = 0;  segsamp = lx; end
if (exist('nfft') ~= 1)      nfft    = 0;  end
if (nfft < segsamp)   nfft = 2^nextpow2(segsamp); end

% Estimate auto and cross-spectra --------------------------
noverlap = fix(overlap/100 * segsamp);
w        = hamming(segsamp);
nadvance = segsamp - noverlap;

pxx = zeros(nfft,1);
pyy = zeros(nfft,1);
pxy = zeros(nfft,1);

ncols = fix(  (lx*nrecs - segsamp) / nadvance) + 1;
ind = [1:segsamp]';
for k=1:ncols
    Xf = fft( x(ind).* w , nfft);
    Yf = fft( y(ind).* w , nfft);
    pxx = pxx + abs(Xf).^2;
    pyy = pyy + abs(Yf).^2;
    pxy = pxy + Yf .* conj(Xf);
    ind = ind + nadvance ;
end
pxx = pxx / (segsamp * ncols * norm(w)^2);
pyy = pyy / (segsamp * ncols * norm(w)^2);
pxy = pxy / (segsamp * ncols * norm(w)^2);


% Estimate cross-bispectrum --------------------------------

plotflag = 0;
[qf,w] = bispecdx (x,x, y, nfft, wind, segsamp, overlap, plotflag);
qf     = fftshift(qf);                   % now in raw 2-D FFT form
bxxy   = qf;

ind    = [1, (2*nfft):nfft:(nfft^2)] - [0,0:nfft-2];
ind    = ind(:);
mu_y = mean(y(:));
qf(ind) = qf(ind) - pxx * mu_y/segsamp;        % w1+w2 = 0 line correction

% Estimate the linear and quadratic TF's  -------------------
hf = pxy  ./ pxx;
qf = qf ./ (2 * pxx * pxx');

% Frequencies are defined only over [-0.25, 0.25]
% The rest is garbage.

ind = (nfft/4+2):(nfft*3/4-1);
win = ones(nfft,nfft);
win(ind,ind) = zeros(nfft/2-2,nfft/2-2);

qf  = qf .* win;

% Estimate IR's
qt = real(ifft2(qf));
ht = real(ifft(hf));

qt = flipud(fliplr(qt));
qt = qt(1:nfft/4, 1:nfft/4);


% ------------ Display estimated TF's and IR's
w1 = [1:nfft/2]/nfft;
w2 = [-nfft/2:nfft/2-1]/nfft;

clf
subplot(221)
semilogy(w1,abs(hf(1:nfft/2))), title('linear TF'), xlabel('f'), grid on
subplot(222)
%contour(abs(fftshift(qf)), 6, w2, w2),  title('quadratic TF')
contour(w2,w2,abs(fftshift(qf)), 6),  title('quadratic TF'), grid on 
xlabel('f1'), ylabel('f2')
subplot(223)
plot(ht),  title('linear part: IR'), xlabel('t'), grid on 
subplot(224)
contour(qt), title('quadratic part: IR'), grid on 
xlabel('t1'), ylabel('t2')
set(gcf, 'Name','Hosa NLTICK')
