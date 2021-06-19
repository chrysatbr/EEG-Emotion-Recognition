function   cmat = cumtrue (ma, ar, norder, nlags, k)
%CUMTRUE Theoretical cumulants of an ARMA model
%	cmat = cumtrue (ma, ar, norder, nlags, k)
%	    ma - MA parameter vector
%	    ar - AR parameter vector;  defaults to [1] (pure MA model)
%	norder - cumulant order: should be 2, 3 or 4;  defaults to 3
%	 nlags - maximum number of cumulant lags to compute.
%	         default value is p+q, where p and q are the AR and MA orders
%	     k - if norder=4, k specifies the 3rd lag of the cumulant;
%	         default value is 0.
%	  cmat - computed cumulant vector or matrix.
%	If norder=2, cmat is a column vector of length 2*nlags + 1,
%	     and consists of C2(m), m=-nlags, .... , nlags
%	If norder=3, cmat is a (2*nlags + 1) by (2*nlags + 1) matrix;
%	      C3(i,j) is returned in cmat(i+nlags+1,j+nlags+1),
%	              i,j = -nlags, ... , nlags
%	      note that the axes point down and right; origin is at center
%	If norder=4, cmat is a (2*nlags + 1) by (2*nlags + 1) matrix;
%	      C4(i,j,k) is returned in cmat(i+nlags+1,j+nlags+1),
%	              i,j = -nlags, ... , nlags
%	      note that the axes point down and right; origin is at center

%
%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.3 $
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

% ----------------- Parameter checks -------------------------
   if (exist('ma') ~= 1)
        error('ma parameter vector must be supplied')
   end
   if (min(size(ma)) ~= 1)
        error('variable ma must be a vector')
   end
   q = length(ma) - 1;
   if (exist('ar') ~= 1) ar = [1]; end
   if (min(size(ar)) ~= 1)
        error('variable ar must be a vector')
   end
   p = length(ar) - 1;

   if (exist('norder') ~= 1) norder = 3; end
   if (norder < 2 | norder > 4)
      error ('norder should be 2, 3 or 4')
   end
   if (exist('nlags') ~= 1) nlags = p+q; end
   if (p == 0) nlags = min(nlags, q); end
   if (exist('k') ~= 1) k = 0; end
   klag3 = k;


% ------------ compute theoretical cumulants ---------------
% c2(i)     = sum_{n} h(n) h(n+i)
% c3(i,j)   = sum_{n} h(n) h(n+i) h(n+j)
% c4(i,j,k) = sum_{n} h(n) h(n+i) h(n+j) h(n+k)

   if (p == 0) h = ma(:);
   else,
       rpoles = abs(roots(ar));
       if (any(rpoles >= 1) )
          error('unstable AR polynomial passed')
       end
       rho = max(rpoles);
       nsamp = max(2*nlags,round(log(0.001)/log(rho)));
       h = filter(ma,ar,[1;zeros(nsamp,1)]);
   end

   maxlag = nlags;
   nlags  = 2*maxlag + 1;
   nlag1  = maxlag + 1;
   nsamp  = length(h);

   if (norder == 2) cmat = zeros(nlags,1);
   else, cmat = zeros(nlags, nlags); end

   if (norder == 2)
      cmat(nlag1) = h'*h;
      for n=1:maxlag
	cmat(nlag1+n) = h(1:nsamp-n)' * h(n+1:nsamp);
        cmat(nlag1-n) = cmat(nlag1+n);
      end
      return
   end

   if (norder == 3)
      for n=-maxlag:maxlag
          z = h*0;
          if (n >= 0) z(1:nsamp-n)  = h(1:nsamp-n) .* h(n+1:nsamp);
          else        z(-n+1:nsamp) = h(-n+1:nsamp) .* h(1:nsamp+n);
          end
          cmat(n+nlag1,nlag1) = h'*z;
          for k=1:maxlag
              cmat(n+nlag1,nlag1+k) = z(1:nsamp-k)' * h(k+1:nsamp);
              cmat(n+nlag1,nlag1-k) = h(1:nsamp-k)' * z(k+1:nsamp);
          end
      end
      return
   end

   if (norder == 4)
      z = h * 0;
      z2 = z;
      k = klag3;
      if (k >= 0) z2(1:nsamp-k) = h(1:nsamp-k) .* h(k+1:nsamp);
      else        z(-k+1:nsamp) = h(-k+1:nsamp) .* h(1:nsamp+k);
      end

      for n=-maxlag:maxlag
          z =  h*0;
          if (n >= 0) z(1:nsamp-n) = z2(1:nsamp-n) .* h(n+1:nsamp);
          else        z(-n+1:nsamp) = z2(-n+1:nsamp) .* h(1:nsamp+n);
          end

          cmat(n+nlag1,nlag1) = h'*z;
          for k=1:maxlag
              cmat(n+nlag1,nlag1+k) = z(1:nsamp-k)' * h(k+1:nsamp);
              cmat(n+nlag1,nlag1-k) = h(1:nsamp-k)' * z(k+1:nsamp);
          end
      end
      return
   end

return
