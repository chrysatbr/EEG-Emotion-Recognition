function [wx,waxis] = wig4c (x0,nfft,sigma,flag)
%WIG4C	 Computes the Sliced Reduced Interference Wigner Trispectrum
%	[wx,waxis] = wig4c (x, nfft, sigma,flag)
%	x     - time series, must be a vector
%	nfft  - FFT length to use; default is the power of 2 just larger
%	        than four times the length of x.
%	sigma - sigma for the Choi-Williams window  [default = 0.05]
%		large values of sigma are equivalent to no windowing;
%	flag  - By default, if signal 'x' is real valued, its analytic form
%	        is used to compute the WD; this helps supress cross terms
%	        around D.C.;  if flag is 0, the analytic form is not used.
%	wx    - the Sliced Reduced Interference Wigner Trispectrum
%	        rows correspond to time, columns to frequencies
%	        time increases with row number, frequencies with col number
%	waxis - the frequency axis associated with the WD
%       It is recommended that the analytic form of the signal be used.

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

% --------------------- parameter checks --------------------------
[m, n] = size(x0);
if (min(m,n) ~= 1)
   disp(['wig4c: input argument x is a ',int2str(m),' by ',int2str(n), ...
         ' array'])
   error('Input argument x must be a vector');
end

if (exist('sigma') ~= 1) sigma = .05; end
if (exist('flag') ~= 1) flag = 1; end
if (all(imag(x0)==0) & flag ~= 0) x0 = hilbert(x0); end

% ------------- find power of two for FFT --------------------------
% signal must be zero-padded to twice the length to avoid aliasing

lx0 = length(x0);
x0 = conv(x0,x0);                             % the basic relationship
lx = length(x0);

lfft = 2^nextpow2(2*lx);                      % minimum FFT length
if (exist('nfft') ~= 1) nfft = lfft; end
if (isempty(nfft)) nfft = lfft; end

if (nfft < 2*lx)
   disp(['WIG4C: FFT length must exceed four times the signal length'])
   disp(['     resetting FFT length to ',int2str(lfft)])
   nfft = lfft;
end

x  = zeros(nfft,1);  x(1:lx) = x0(:);       cx = conj(x);
wx = zeros(nfft,nfft);
L1 = lx - 1;

% --------- compute r(tau,t) = cx(t-tau/2) * x(t+tau/2) -------------

for n=0:lx
   indm = max(-n,-L1+n) : min(n,L1-n);
   indy = indm + (indm < 0) * nfft ;   % output indices y(m;n)
   y = zeros(nfft,1);
   y(indy + 1) = x(n+indm + 1) .* cx(n-indm + 1);
   wx(:,n+1)   = y;
end


% ---------- A(tau,theta) = IFT (t-->theta) r(tau,t) ---------------
% compute the ambiguity function

wx  = ifft(wx.').';

% ---------------- Apply the Choi Williams window  ---------------
if (finite(sigma))
     win =   [ (1:nfft)-nfft/2-1 ]' * [ (1:nfft)-nfft/2-1] ...
     		/ nfft;
     win = fftshift( exp (- win.^2 / sigma) );
else
     win = ones(nfft,nfft);
end
wx  = wx .* win;

% -------------- WD(t,f) = FT A(theta,tau) -------------------
wx  = fft2(wx);         % fft along both time and frequency


% --------------- Display the WD -----------------------------
% note the frequency scaling by 2M
% WD(f,t) = X(2f,t), where X(f,t) = FT ( r(tau,t) )

nfftby2 = nfft/2;
wx  = real(wx(:,1:2:lx).'); % throw away odd indices
wx = wx(:,[nfftby2+1:nfft,1:nfftby2]) ;
waxis = [-nfftby2:nfftby2-1] / (2*nfft);
taxis = 1:lx0;

%contour(abs(wx),8,waxis,taxis), grid,
contour(waxis,taxis, abs(wx),8), grid on 
ylabel('time in samples')
xlabel('frequency')
title(['CWT s=',num2str(sigma)])
set(gcf,'Name','Hosa WIG4C')
