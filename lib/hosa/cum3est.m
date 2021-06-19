function   y_cum = cum3est (y, maxlag, nsamp, overlap, flag, k1)
%CUM3EST Third-order cumulants.
%	Should be invoked via "CUMEST" for proper parameter checks
%	y_cum = cum3est (y, maxlag, samp_seg, overlap, flag, k1)

%	y_cum = cum3est (y, maxlag, samp_seg, overlap, flag, k1)
%	       y: input data vector (column)
%	  maxlag: maximum lag to be computed
%	samp_seg: samples per segment
%	 overlap: percentage overlap of segments
%	   flag : 'biased', biased estimates are computed  [default]
%	          'unbiased', unbiased estimates are computed.
%	      k1: the fixed lag in c3(m,k1): see below
%	   y_cum:  estimated third-order cumulant,
%	           C3(m,k1)  -maxlag <= m <= maxlag

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.4 $
%  A. Swami   January 20, 1993

% Modified Jan 20, 94 to handle complex case properly.
%  c3(i,j) := E x^*(n) x(n+i) x(n+j)   (x assumed zero mean)

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

%  c3(i,j) := E x^*(n) x(n+i) x(n+j)   (x assumed zero mean)

%---------------- Parameter checks done by CUMEST --------------
   [n1,n2]  = size(y);
   N        = n1*n2;
   minlag   = -maxlag;
   overlap  = fix(overlap/100 * nsamp);
   nrecord  = fix( (N - overlap)/(nsamp - overlap) );
   nadvance = nsamp - overlap;

   y_cum = zeros(maxlag-minlag+1,1);

   ind = (1:nsamp)';
   nlags = 2 * maxlag + 1;
   zlag  = 1 + maxlag;
   if (flag(1) == 'b' | flag(1) == 'B')
   	scale = ones(nlags,1)/nsamp;
   else
       lsamp = nsamp - abs(k1);
       scale = [lsamp-maxlag:lsamp,lsamp-1:-1:lsamp-maxlag]';
       [m2,n2] = size(scale);
       scale = ones(m2,n2) ./ scale;
   end

   for i=1:nrecord
       x = y(ind); x = x(:) - mean(x);     % make sure we have a col vec
       cx = conj(x);
       z = x*0;

%                     create the "IV" matrix: offset for second lag

       if (k1 >= 0) z(1:nsamp-k1)  = x(1:nsamp-k1,:) .* cx(k1+1: nsamp,:);
       else         z(-k1+1:nsamp) = x(-k1+1:nsamp)  .* cx(1:nsamp+k1);
       end

%                     compute third-order cumulants

       y_cum(zlag)  =  y_cum(zlag) + z' * x;

       for k = 1:maxlag
           y_cum(zlag-k) = y_cum(zlag-k) + z([k+1:nsamp])' * x([1:nsamp-k]);
           y_cum(zlag+k) = y_cum(zlag+k) + z([1:nsamp-k])' * x([k+1:nsamp]);
       end

       ind = ind + nadvance;
   end

   y_cum = y_cum .* scale / nrecord;

return
