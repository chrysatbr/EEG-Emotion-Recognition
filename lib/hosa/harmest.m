function [Pxx,ar1,ar2] = harmest(y,maxlag,p_order,flag,nfft,norder)
%HARMEST Frequencies of harmonics in colored Gaussian noise.
%	[Pxx,ar1,ar2] = harmest(y,maxlag,p_order,flag,nfft,norder)
%	      y - data  matrix [nsamp x nrecs]
%	 maxlag - number of cumulant lags to compute [default = nsamp/12]
%	p_order - order to use (dimension of signal subspace)
%	          user will be prompted if p_order <= 0
%	   flag - 'biased' or 'unbiased'
%	   nfft - fft length for spectra [default = 256]
%	 norder - cumulant order to use: 2 or 4 [default = 4]
%
%	   Pxx  - a  nfft x 7 matrix of spectral estimates
%	   ar1  - estimated parameters for the AR method.
%	   ar2  - estimated parameters for the min-norm method
%
%	The seven columns of Pxx contain the spectral estimates based on
%	the Eigenvector, Music, Pisarenko, ML, AR, periodogram methods,
%	and the minimum-norm method.
%

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.6 $
%  A. Swami   January 20, 1993;  August 8, 1994

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

% --------- Parameter checks ---------------------------------

[nsamp,nrecs] = size(y);
if (nsamp == 1) nsamp = nrecs; nrecs = 1;  y = y(:); end
if (exist('p_order') ~=1)    p_order = 0;            end
if (exist('flag') ~= 1)      flag    = 'biased';     end
if (flag(1) ~= 'b')          flag    = 'unbiased';   end
if (exist('maxlag') ~= 1)    maxlag  = nsamp/12;     end
if (exist('nfft') ~= 1)      nfft    = 256;          end
if (exist('norder') ~= 1)    norder  = 4;            end
if (norder ~= 2 & norder ~= 4)
   error('cumulant order - norder - should be 2 or 4')
end

Pxx = zeros(nfft,7);

% ------------ estimate second/fourth-order cumulants --------------
  if (norder == 4)
    M = maxlag;
    cum_4y = cumest(y,4,M,nsamp,0,flag,0,0);
    cum_4y = (cum_4y(M+1:2*M+1) + conj(cum_4y(M+1:-1:1)) )/2;
    cum_4y = cum_4y/cum_4y(1);
    Amat   = toeplitz(conj(cum_4y),cum_4y);
  else
    Amat = zeros(maxlag+1,maxlag+1);
    ind = [0:-1:-maxlag];
    for n=maxlag+1:nsamp
        xf = y(n+ind,:);
        xb = conj(flipud(xf));
        Amat = Amat + xf * xf' + xb * xb';
    end
    mu = mean(mean(y));
    Amat = Amat - mu*mu';                   % remove mean
    Amat = Amat / ( nrecs*(nsamp-maxlag) );
  end
  [U,S,V] = svd(Amat); sval = diag(S);

% -----------  how many harmonics ? -------------------------

  hold off, clf,
  set(gcf,'Name', ...
      ['Hosa HARMEST -  cum-order=', int2str(norder)] )
  ls = length(sval);
  plot(1:ls, sval,'-', 1:ls,sval,'go'),  grid on 
  title(['sval: cum-',int2str(norder)])
  drawnow
  p = p_order;

  while (p <= 0 | p > ls )
        p = input([' enter order to use [1,' int2str(ls) '] ---> ']);
        if (isempty(p)) p = 0; end
  end

 % ---------- Noise subspace methods --------------

   M    = maxlag+1;
   nvec = M - p;  mfft = nfft/2;
   wte  = ones(nvec,1) ./ sval(p+1:M);

   Pvm  = V(:,p+1:M) * V(:,p+1:M)';
   Pve  = V(:,p+1:M) * diag(wte) * V(:,p+1:M)';
   psum = zeros(nfft,2);
   psum(1,1) = sum(diag(Pve));
   psum(1,2) = sum(diag(Pvm));
   for k=1:maxlag
       psum(k+1,1)      = sum(diag(Pve,k));
       psum(nfft-k+1,1) = sum(diag(Pve,-k));
       psum(k+1,2)      = sum(diag(Pvm,k));
       psum(nfft-k+1,2) = sum(diag(Pvm,-k));
   end
   Pxx(:,1:2)  = ones(nfft,2) ./ real(fft(psum,nfft));

   [u1,s1,v1] = svd(Amat(1:p+1,1:p+1));  v1 = v1(:,p+1);
   Pxx(:,3)   = ones(nfft,1) ./ abs(fft(conj(v1),nfft)).^2;

% ------------ signal subspace method -----------------------------
% mlcapon
   for k=1:p
      Pxx(:,4) = Pxx(:,4) + abs(fft(conj(V(:,k))/sqrt(sval(k)),nfft)) .^ 2;
   end

% --------- AR method -----------------------------------------------
% --------- SVD low-rank approximation (a la Cadzow) ----------------
   Amat = U(:,1:p) * diag(sval(1:p)) * (V(:,1:p))';
   Ahat = [];
   [m,n] = size(U);
   for k=p+1:n
       Ahat = [Ahat; Amat(:,k-p:k)];
   end

%                                      force unity modulus solution
   avec = [1; -Ahat(:,2:p+1) \ Ahat(:,1)];
   avec = conj(avec);
if (exist('debug') == 1)
   ncol = floor((p+2)/2);
   nflip = floor((p+1)/2);

   Ahat(:,1:nflip) = Ahat(:,1:nflip) + Ahat(:,p+1:-1:p+2-nflip);
   avec = Ahat(:,2:ncol) \ Ahat(:,1);
   avec = [1; -avec];
   avec = [avec; avec(nflip:-1:1)];
end
ar1 = avec;
   Pxx(:,5) = ones(nfft,1) ./ abs(fft(avec,nfft)) .^2;

% ----------- The periodogram estimate  ----------------------

  ind = [1:nsamp]';
  for k=1:nrecs
      ys = y(ind);
      ys = ys - mean(ys);
      Yf = fft(ys,nfft) / nsamp;
      Pxx(:,6) = Pxx(:,6) + Yf .* conj(Yf);
      ind = ind + nsamp;
  end
  Pxx(:,6) = Pxx(:,6) / nrecs;
  Pxx(1,6) = Pxx(2,6);                 % dynamic range problems, else

% ----------------- Minimum-norm method -------------------

 gvec  = U(1,1:p).';
 gmat  = U(2:length(sval),1:p);
 ar2   = [1; - conj(gmat) * gvec / (1 - gvec' * gvec)];

 Pxx(:,7) = fft(ar2, nfft);
 Pxx(:,7) = ones(nfft,1) ./ ( Pxx(:,7) .* conj(Pxx(:,7)) ) ;

% ----------- Display estimates ---------------------

 waxis = [-mfft:mfft-1]/nfft;
 Pxx   = Pxx([mfft+1:nfft,1:mfft],:);

% - scale to max abs of unity for plots only

   spmax = max(Pxx);
   spmax = ones(1,length(spmax)) ./ spmax;

   clf,
   subplot(421),   plot(1:ls, sval,'-', 1:ls,sval,'g.'), grid on
   title(['sval: cum-',int2str(norder)])

subplot(422), plot(waxis,10*log10(Pxx(:,6)*spmax(6))), title('pxx'),grid on
subplot(423), plot(waxis,10*log10(Pxx(:,1)*spmax(1))), ylabel('eig'),grid on
subplot(424), plot(waxis,10*log10(Pxx(:,2)*spmax(2))), ylabel('music'),grid on
subplot(425), plot(waxis,10*log10(Pxx(:,3)*spmax(3))), ylabel('pisar'),grid on
subplot(426), plot(waxis,10*log10(Pxx(:,4)*spmax(4))), ylabel('ml'),grid on
subplot(427), plot(waxis,10*log10(Pxx(:,5)*spmax(5))), ylabel('ar'),grid on
subplot(428), plot(waxis,10*log10(Pxx(:,7)*spmax(7))), ylabel('min-norm'),
grid on

return






