function [sg, sl] = glstat(x,cparm,nfft)
%GLSTAT	Detection statistics for Hinich's Gaussianity and linearity tests
%	[sg, sl] = glstat(x,cparm,nfft)
%	x     - time series
%	cparm - resolution (smoothing) parameter c;
%	        0.5 < c < 1.0 [default = 0.51]
%	        increasing c reduces the variance,
%	        but increases the bias (window bandwidth).
%	nfft  - fft length [default = 128]
%	sg    - statistic for the Gaussianity test
%	        [observed chi_sq value, dg, PFA]
%	sl    - statistic for the linearity test
%	        [R_estimated, lambda, R_theoretical]
%
%	Under the assumption of Gaussianity, the test statistics S is
%	      chi-squared distributed with df degrees of freedom.
%	Under the assumption of non-zero skewness and linearity, the test
%	      statistic R should be approximately equal to the inter-quartile
%	      range of a chi-squared distribution with two degrees of freedom,
%	      and non-centrality parameter lambda.
%

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.4 $
%  A. Swami   January 20, 1993.
%  A. Swami,  October 14, 1994, Revised.
%  A. Swami,  February 23, 1995, added `N' printout in linearity test

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

% ---------------- parameter checks: ----------------------

   if (min(size(x)) ~= 1)
      error('x: should be a vector')
   end
   nsamp = length(x);
   if (exist('nfft') ~= 1) nfft = 128; end
   nrecs = floor(nsamp/nfft);
   nrecs = max(nrecs,1);
   ksamp = min(nfft,nsamp);
   if (nfft  > nsamp)
      disp([' glstat results are unreliable if zero-padding is severe: '])
      disp(['    fft length= ',int2str(nfft),' data length=',int2str(nsamp)])
   end

   if (exist ('cparm') ~= 1) cparm = 0.51; end
   if (cparm >= 1.0)
      error('cparm  cannot be greater than or equal to 1.')
   end
   if (cparm <=0.5 & nrecs == 1)
      error('cparm: for single segments: allowed range is (0.5,1.0) ')
   end

   M = fix(nfft^cparm);  M = M + 1 - rem(M,2);
   Msmuth = M;
% ----------------- estimate raw bispectrum ------------------------

   mrow = fix(nfft/3)+1; ncol = fix(nfft/2)+1;
   F = zeros(mrow,ncol);
   S = zeros(nfft,1);

   mask  = hankel ([1:mrow],[mrow:mrow+ncol-1]);
   ind   = 1:ksamp;
   for k=1:nrecs
       y = x(ind);
       xf = fft(y-mean(y), nfft);
       xc = conj(xf);
       S  = S + xf .* xc;                          % power spectrum
       F = F + xf(1:mrow) * xf(1:ncol).' .* ...
           reshape (xc(mask), mrow, ncol) ;
      ind = ind + ksamp;
   end
   F = F /(nfft*nrecs);                   % `Raw' bispectrum F as in Eq (2.3)
   S = S /(nfft*nrecs); 	          % `Raw' power spectrum

% ---------- zero out area outside principal domain ---------------

   ind = (0:(mrow-1))';
   F(1:mrow,1:mrow) = triu(F(1:mrow,1:mrow));
   Q = ones(mrow, ncol);
   Q(1:mrow,1:mrow) = triu(Q(1:mrow,1:mrow)) + eye(mrow);

% the 2f1 + f2 = 1 line:

   r = ( rem(nfft,3) == 2);        % ans oct 94: also -r next line
   for k=mrow+1:ncol-r
       j = k-mrow;
       Q(mrow-2*j+2:mrow, k) = zeros(2*j-1,1);
       F(:, k) = F(:,k) .* Q(:,k);
       Q(mrow-2*j+1, k) = 2;                      % factor 2 on boundary
   end

   F = F(2:mrow,2:ncol);      % in principal domain, no dc terms
   Q = Q(2:mrow,2:ncol);
   mrow = mrow-1;
   ncol = ncol-1;

% -------- smooth the estimated spectrum and bispectrum ------------
%  Msmuth x Msmuth box-car smoother
%  only every Msmuth term in the smoothed output is retained
%  "independent" estimates

   M = Msmuth;
   mby2 = (M+1)/2;
   m1 = rem(mrow,M); m2=rem(ncol,M);
   m1 =  - m1 + M * (m1 ~= 0);
   m2 =  - m2 + M * (m2 ~= 0);
   F  = [F, zeros(mrow,m2); zeros(m1,ncol+m2)];
   Q  = [Q, zeros(mrow,m2); zeros(m1,ncol+m2)];

   k  = size(F)/Msmuth;
   k1 = k(1); k2 = k(2);

% Apply the box car smoother
%  can replace the next ten lines or so with
% B =kron(eye(k1),ones(1,Msmuth)) * F * kron(eye(k2),ones(Msmuth,1))/Msmuth^2;
% Q =kron(eye(k1),ones(1,Msmuth)) * Q * kron(eye(k2),ones(Msmuth,1));

ind = 1:Msmuth;
B = zeros(k1,k2); Q1 = B;
for i=1:k1
    for j=1:k2
        t = F( (i-1)*Msmuth+ind, (j-1)*Msmuth+ind );
        B(i,j) = mean(t(:));
        t = Q( (i-1)*Msmuth+ind, (j-1)*Msmuth+ind );
        Q1(i,j) = sum(t(:));
    end
end
Q = Q1;

% --------------
% At this point B corresponds to B  in eq (2.6) of Hinich's paper
%    and Q corresponds to the definition following eq (2.8)
% --------------
   M   = Msmuth;
   S   = conv(S, ones(M,1))/M;
   S   = S(M+1:M:M+nfft-1);
   S1  = B*0;
%   S1  = S(1:k1) * S(1:k2)' .* hankel(S(1:k1),S(k1:k1+k2-1));
   S1  = S(1:k1) * S(1:k2)' .* hankel(S(2:k1+1),S(k1+1:k1+k2)); 
   S1  = ones(k1,k2) ./ S1;
   ind = find(Q > 0);
   Q   = Q(ind);
   Bic = S1(ind) .* abs(B(ind)).^2;          % squared bicoherence

%-------------------------------------------------------------------
%----------------- Gaussianity (non-skewness) test -----------------
%-------------------------------------------------------------------

   scale      = 2 * (Msmuth^4) / nfft;
   Xstat      = scale * Bic ./ Q;       % 2 |X(m,n)|^2 in (2.9)

   df      = length(Q) * 2;             % degrees of freedom
                                        % asymptotic value is (nfft/M)^2 / 6
   chi_val = sum(Xstat);                % observed value of chi-square r.v.

% Reference:
% Sankaran's approximation given in [7.9.4], p 197, of
% J.K. Patel and C.B. Read,   M. Dekker, New York, 1982.
% `Handbook of the Normal Distribution'
% the cdf of chi(df,lam,x) is approximated by phi(y) [normal distro cdf]
% where y = [ (x/(df+lam))^h - a ] / b;
% and h, a, b are defined below

% For a given df, and chi-sq value, the PFA is the probability that
% we will be wrong in accepting the alternate hypothesis; i.e.,
% data are non-skewed.

   h = 1/3;    b = sqrt(2/df) / 3;    a = (1-b^2);
   y = (chi_val/df)^(1/3);   y = (y-a)/b;
   pfa = 0.5 * erfc(y/sqrt(2));
   pfa = round(pfa*10000)/10000;     % round to five decimal places

   disp(['Test statistic for Gaussianity is ', num2str(chi_val), ...
         ' with df = ',int2str(df), ', Pfa = ',num2str(pfa)])
   sg  = [chi_val, df, pfa];

% -----------------------------------------------------------------
% ---------------------- linearity test ---------------------------
% -----------------------------------------------------------------

% ---- find the first and third quartiles of 2|Xmn|^2 = Xstat

   ind   = find(Q == M^2);
   rtest = Xstat(ind);
   rtest = sort(rtest);
   l1    = length(rtest);
 if (l1 < 4)
    fprintf('glstat: cparm (%g) is too large or nfft (%g) is too small', ...
    	     cparm,nfft);
    fprintf('\n   estimated interquartile range (R) set to NaN \n');
    Rhat = NaN;
 else
   lo1 = fix(l1/4);    lo2 = fix(l1/4+0.8);
   hi1 = fix(3*l1/4);  hi2 = fix(3*l1/4+0.8);
   r1  = (rtest(lo1)+rtest(lo2))/2;
   r3  = (rtest(hi1)+rtest(hi2))/2;
   Rhat  = r3 - r1;                 % sample inter-quartile range
end

% ----- estimated lambda: a scaled version of mean FD skewness
%   lam = 2/(df * Msmuth^2) * sum ( Q.*Xstat - 2);    % Hinich's (4.2)
    lam = mean(Bic) * 2 * Msmuth^2 / nfft;            % equivalent

%-----  find the `theoretical inter-quartile' range of chi-sq(df=2,lam)
% See reference cited earlier.

    df = 2;                                 % for our problem, df=2
    h  = 1/3 + (2/3) * ( lam / (df + 2*lam) )^2;
    a  = 1 + h * (h-1) * (df + 2*lam) / ( (df+lam)^2 ) ...
        - h * (h-1) * (2-h) * (1-3*h) * (df + 2*lam)^2 / (2 * (df+lam)^4);
    b  = h * sqrt(2 * (df+2*lam) ) / (df+lam)  ...
        * (  1 - (1-h) * (1-3*h) * (df+2*lam) / ( 2*(df+lam)^2) );

    yn   = sqrt(2) * erfinv(2*[-.25,.25]);       % quartiles of Phi
    xc   = (a+b*yn).^(1/h) * (df + lam);         % quartiles of Chi
    Rth  = xc(2)-xc(1);                          % IQ range

 sl      = [Rhat, lam, Rth];
 disp(['Linearity test:  R (estimated) = ',num2str(Rhat), ...
        ', lambda = ', num2str(lam), ', R (theory) = ', num2str(Rth), ...
	', N = ',int2str(length(rtest)) ])
return
