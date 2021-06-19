function [ar_vec,Bspec] = qpctor(y,maxlag,ar_order,nfft,nsamp,overlap,flag)
%QPCTOR	 Estimation of quadratic-phase coupling using third-order cumulants.
%	[arvec,Bspec] = qpctor(y, maxlag,ar_order, nfft, nsamp,overlap,flag)
%	y        - data vector or matrix
%	maxlag   - maximum number of cumulant lags to use  [default: 12]
%	ar_order - AR order                       [default: prompt user]
%	nfft     - fft length for bispectrum      [default = 64]
%	nsamp    - number of samples per record   [row dimension of y]
%	overlap  - percentage overlap              [default = 0]
%	flag     - 'biased' or 'unbiased'         [default = 'biased']
%	arvec    - estimated AR parameters
%	Bspec    - estimated bispectrum, an nfft/2 by nfft upper
%	         - triangular array, corresponding to the
%	         - 0 <= f2 <= f1 < 0.5

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.8 $
%  A. Swami   January 20, 1993.

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

% ----- Parameter checks -------------------------------

[ly,nrecord]    = size(y);
if (ly == 1) y = y'; ly = nrecord; nrecord = 1;   end   %vectors must be cols

if (exist('maxlag')   ~= 1)  maxlag   = 12;       end
if (exist('ar_order') ~= 1)  ar_order = 0;        end
if (exist('nfft') ~=1)       nfft     = 64;       end
if (nfft <= 0)               nfft     = 64;       end
if (exist('nsamp') ~= 1)     nsamp    = ly;       end
nsamp   = min(ly, max(nsamp,  0));
if (exist('overlap') ~= 1)   overlap  = 0; 	  end
if (nrecord ~= 1)            overlap  = 0;        end   %input is a matrix
overlap = min(99, max(overlap,0));
if (exist('flag') ~= 1)          flag = 'biased'; end
if (flag(1:1) ~= 'b')            flag = 'unbiased'; end

% ------- conventional power spectrum
  pxx = zeros(2*nfft,1);
  for k=1:nrecord,
      pxx = pxx + abs( fft( y(:,k), 2*nfft ) ).^2 ;
  end
  pxx = pxx(1:nfft);  pxx = pxx/max(pxx);
  pxx(1) = pxx(2);                          % dynamic range problems

  waxis = [0:nfft-1]/(2*nfft);
  clf, subplot(221),
  semilogy(waxis, pxx); grid on 
  title('Power Spectrum'),      xlabel('frequency')
  set(gcf,'Name','Hosa QPCTOR')

% ----- estimate third-order cumulants  -----------------------------
  ycum = cum3est(y(:),maxlag,nsamp,overlap, flag,0);   % C(m,0) slice
  ycum = ycum(2*maxlag+1:-1:1);                        % C(m,m) slice

% ------ Set up the linear system of equations -----------------------
%   \sum_{k=0}^{p} a(k) C(k-m,k-m) = beta   m=0
%                                  = 0      m=1,..,p
%  where a(0) = 1; order p = 6N (N = number of frequency triplets)
%    H(f) = 1/A(f);  C(f1,f2) = beta H(f1) H(f2) H(-f1-f2)

  Amat = toeplitz(ycum(maxlag+1:-1:1), ycum(maxlag+1:2*maxlag+1));

% ------ How many triplets (AR order) ? ---------------------------------

  s = svd(Amat);

  subplot(222),
  plot(s,'-'), hold on, plot(s,'o'), hold off, grid on
  title('Singular values')

  while (ar_order <= 0)
     txt = ['specify AR order (p) to use: [1,' int2str(length(s)) '] ---->'];
     ar_order = input(txt); 
     if (isempty(ar_order)) ar_order = 0; 	end
     disp(' ')
  end

% ------ Obtain the AR vector ---------------------------
   ar_vec = [1; -Amat(:,2:ar_order+1) \ Amat(:,1)];

% ----- compute the AR transfer function H(w) = 1/A(w)
% ----- then, bispectrum B(w1,w2) = H(w1)H(w2)conj(H(w1+w2))

  Hf   = fft(ar_vec, 2*nfft);
  Hf   = Hf(1:nfft);
  Hf   = ones(nfft,1) ./ Hf;
  nfft1 = fix(nfft/2) + rem(nfft,2);
  Hf1  = Hf(1:nfft1);
  ind  = [nfft1:nfft, 1:nfft1-1];

  Bspec = (Hf1 * Hf.' ) .* conj(hankel(Hf1, Hf(ind)));
  Bspec = triu(Bspec);


% ------- displays ----------------

  [y,j] = max(abs(Bspec));  [val,i] = max(y);
  f1    = waxis(i);
  f2    = waxis(j(i));
  disp(['Maximum of bispectrum:  B(',num2str(f1),',',num2str(f2),') = ', ...
       num2str(val)])

  subplot(212),
%  contour(abs(Bspec),8,waxis,waxis(1:nfft/2)), grid
  contour(waxis,waxis(1:nfft/2), abs(Bspec),8), grid on
  hold on
  plot(waxis(1:nfft/2),waxis(1:nfft/2),'--')     % the diagonal line
  hold off
  xlabel('f1'), ylabel('f2')
  title('Estimated Parametric Bispectrum')

return
