function avec = arrcest(y,p,q, norder,maxlag,samp_seg,overlap,flag)
%ARRCEST AR parameter estimates based on cumulants.
%	avec = arrcest (y,p,q, norder,maxlag,samp_seg,overlap,flag)
%	    y : time-series (vector or matrix)
%	    p : AR order
%	    q : MA order
%	norder: cumulant order: 2, 3 or 4         [default = 2]
%	         -3 : use correlations and 3rd order cumulants
%	         -4 : use correlations and 4th order cumulants
%	  maxlag:  maximum cumulant lag to be used  [default = p+q]
%	samp_seg: samples per segment for estimating cumulants
%	                           [default = length of y]
%	overlap : percentage overlap of segments    [default = 0]
%	   flag : 'biased' or 'unbiased'            [default = 'biased']
%	avec    : estimated AR parameter vector

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.4 $
%  A. Swami, January 20, 1992

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

% ------------------- parameter checks --------------------
 if (nargin < 2)
    error('insufficient number of parameters')
 end

 [nsamp, nrecs] = size(y);
 if (nsamp == 1) nsamp = nrecs; nrecs = 1; y = y.';  end

 if (p <= 0)                  avec = 1;    return,    end
 if (exist('q') ~= 1)           q = 0;                end
 if (q < 0)
    error('MA order q must be non-negative'),
 end
 if (exist('norder') ~= 1)      norder = 2;           end
 if ( norder ~=2 & abs (norder) ~=3 & abs(norder) ~= 4)
    error('norder must be 2, 3, 4, -3 or -4')
 end
 if (exist('maxlag') ~=1 )       maxlag = p + q;      end
 if (maxlag < p+q),
    disp(['ARRCEST: maxlag changed from ', int2str(maxlag), ...
      ' to ',int2str(p+q)])
    maxlag = p + q;
 end
 if (exist('samp_seg') ~= 1) samp_seg = nsamp; end
 if (exist('overlap') ~= 1)   overlap = 0;     end
 overlap = max(0, min(overlap, 99));
 if (exist('flag') ~= 1)      flag = 'biased'; end

 if (nrecs > 1) overlap = 0; samp_seg = nsamp; end

  minlag = -maxlag;
  nlags = maxlag - minlag + 1;
  Amat  = [];
  rvec  = [];

% --------------------  estimate cumulants --------------------

  if (norder ~= 2)
    kslice1 = [q-p, q];
    kslice2 = [0,0];

    kslice = (kslice1(2) - kslice1(1) + 1) * (kslice2(2) - kslice2(1) + 1);
    cum_y = zeros(nlags,kslice);

    morder = abs(norder);
    kloc = 0;
    for k1 = kslice1(1) : kslice1(2)
       for k2 = kslice2(1) : kslice2(2)
           kloc = kloc + 1;
           cum_y(:,kloc) = cumest(y, morder, maxlag, samp_seg, overlap,  ...
                                  flag, k1, k2);
       end
     end

%--------------- set up cumulant-based `normal' equations -------------
       Amat = hankel(cum_y(q-p+1-minlag+1:nlags-p, 1), ...
                  cum_y(nlags-p:nlags-1,1) );
       rvec = cum_y(q+1-minlag+1:nlags, 1);
       for k=2:kslice
           Amat = [Amat; hankel(cum_y(q-p+1-minlag+1:nlags-p, k), ...
                          cum_y(nlags-p:nlags-1,k) ) ];
           rvec = [rvec; cum_y(q+1-minlag+1:nlags, k) ];
       end
       rvec = -rvec;
   end

% ------- append correlation-based normal equations ------------
   if (norder == 2 | norder < 0)
       cor_y = cumest(y,2,maxlag,samp_seg,overlap,flag);
       AR = hankel(cor_y(q-p+1-minlag+1:nlags-p, 1), ...
                  cor_y(nlags-p:nlags-1,1) );
       br = -cor_y(q+1-minlag+1:nlags, 1);

       Amat = [Amat; AR];  rvec = [rvec; br];
   end

% ------------- compute LS estimate -------------------
   avec = Amat \rvec;
   avec = [1; avec(p:-1:1)];

return
