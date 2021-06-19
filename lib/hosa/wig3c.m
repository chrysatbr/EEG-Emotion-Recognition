function [wx, waxis] = wig3c (x0,nfft,sigma,flag)
%WIG3C	 Computes the (f,f) slice of the Third-order Wigner Distribution,
%        with Choi-Williams smoothing
%	[wx,waxis] = wig3c (x, nfft, sigma,flag)
%	x     - time series, must be a vector
%	nfft  - FFT length to use; default is the power of 2 just larger
%	        than thrice the length of x.
%	sigma - sigma for the Choi-Williams window; default = 0.05
%	        large values of sigma are equivalent to no windowing;
%	        values greater than or equal to 100 are treated as infinity
%	flag  - By default, if signal 'x' is real valued, its analytic form
%	        is used to compute the WD; this helps supress cross terms
%	        around D.C.;  if flag is 0, the analytic form is not used.
%	wx    - The Thrid-Order Wigner distribution (TWD)
%	        The TWD is a function of time and two-frequencies;
%	        here we compute the `diagonal slice', f1=f2.
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
   disp(['wig3c: input argument x is a ',int2str(m),' by ',int2str(n), ...
         ' array'])
   error('Input argument x must be a vector');
end

if (exist('sigma') ~= 1) sigma = .05; end

if (exist('flag') ~= 1) flag = 1; end
if (all(imag(x0)==0) & flag ~= 0) x0 = hilbert(x0); end

% ------------- find power of two for FFT --------------------------
% signal must be zero-padded to twice the length to avoid aliasing

lx = length(x0);
lfft = 2 ^ nextpow2(3*lx);
if (exist('nfft') ~= 1) nfft = lfft; end
if (isempty(nfft)) nfft = lfft; end


if (nfft < 3*lx)
   disp(['wig3c: FFT length must exceed thrice the signal length'])
   disp(['     resetting FFT length to ',int2str(lfft)])
   nfft = lfft;
end

x = zeros(nfft,1);   x(1:lx) = x0(:);               cx = conj(x);
wx = zeros(nfft,nfft);
L1 = lx-1;

% --------- compute r(tau,t) = cx(t-tau/2) * cx(t+tau/2) -------------

for n = 0:L1
    for m = -L1+n:n
        ind3k = max(m-n,-L1+n+2*m) : min(L1+m-n,n+2*m);
        s = sum (x(n-m+ind3k +1) .* x(n+2*m-ind3k + 1) );
        m1 = m + (m<0) * nfft + 1;
        wx(m1, n+1) = s * cx (n-m+1);
    end
end

% ----------- AF (tau,w) = IFT (t-->w) r(tau,t) -------------------
wx = ifft(wx.').';

% ----------- Apply the Choi-Williams window function -------------

if (finite(sigma))      % the window function in the ambiguity domain
     win =   [ (1:nfft)-nfft/2-1]' * [ (1:nfft)-nfft/2-1] /nfft;
     win = fftshift( exp (- win.^2 / sigma) );
else
     win = ones(nfft,nfft);
end
wx  = wx .* win;        % apply the window

% ---------------- FT to obtain the WB ---------------------------
wx  = fft2(wx);         % fft along both time and frequency

% ----------- display the WD ---------------------------------------
%  note the frequency scaling by 2M
%  WD(f,t) =  X(2f,t),  where X(f,t) = FT r(tau,t) .

nfftby2 = nfft/2;
wx    = wx(:,1:lx).' ;
wx = wx(:,[nfftby2+1:nfft,1:nfftby2]) ;
waxis = [-nfftby2:nfftby2-1] / (2*nfft);
taxis = 1:lx;

%contour(abs(wx),8,waxis,taxis), grid,
contour(waxis,taxis, abs(wx),8), grid on 
ylabel('time in samples')
xlabel('frequency')
title(['CWB f1=f2  s=',num2str(sigma)])
set(gcf,'Name','Hosa WIG3C')
