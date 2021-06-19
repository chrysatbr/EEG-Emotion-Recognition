function [avec, bvec] = armarts(y,p,q, norder,maxlag,samp_seg,overlap,flag)
%ARMARTS Estimates ARMA parameters via the residual time series method.
%	[avec, bvec] = armarts(y,p,q, norder,maxlag,samp_seg,overlap,flag)
%	      y : time-series (vector or matrix)
%	      p : AR order
%	      q : MA order
%	  norder: cumulant order: -4, -3, 3 or 4      [default = 3 ]
%	          -3 : use correlations and 3rd order cumulants
%	          -4 : use correlations and 4th order cumulants
%	  maxlag:  maximum cumulant lag to be used    [default = p+q]
%	samp_seg: samples per segment for estimating cumulants
%	                                              [default = length of y]
%	overlap : percentage overlap of segments      [default = 0]
%	   flag : 'biased' or 'unbiased'              [default = 'biased']
%	   avec : estimated AR parameter vector
%	   bvec : estimated MA parameter vector

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc.
%       $Revision: 1.4 $
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

%--------------- parameter checks ---------------------

 if (nargin < 3)
    error('insufficient number of parameters')
 end

 [nsamp, nrecs] = size(y);
 if (nsamp == 1) nsamp = nrecs; nrecs = 1; y = y.'; end

 if (p < 0)
    error('AR order cannot be negative')
 end

 if (q < 0)
    error('MA order cannot be negative')
 end

 if (exist('norder') ~= 1) norder = 3; end
 if ( abs (norder) ~=3 & abs(norder) ~= 4)
    error('norder must be  3, 4, -3 or -4')
 end

 if (exist('samp_seg') ~= 1)    samp_seg = nsamp; end
 if (exist('overlap') ~=1)       overlap = 0;     end
 overlap = max(0,min(overlap,99));
 if (exist('flag')  ~= 1)        flag = 'biased'; end
 maxlag0 = p + q;
 if (exist('maxlag') ~= 1)      maxlag = maxlag0; end
 if (maxlag < maxlag0)
    disp(['maxlag changed from ',int2str(maxlag), ' to ', int2str(maxlag0)])
    maxlag = maxlag0;
 end

 if (nrecs > 1)  overlap = 0; samp_seg = nsamp; end

%-------- estimate AR parameters, compute residual TS --------------
  avec  = [1];
  if ( p > 0)
     avec = arrcest(y,p,q, norder,maxlag,samp_seg,overlap,flag);
     y  = filter(avec,[1],y);
  end

%-------- estimate MA parameters ------------------------------------


  bvec = [1];
  if (q > 0)
     bvec = maest(y,q, abs(norder),samp_seg,overlap,flag);
  end

return
