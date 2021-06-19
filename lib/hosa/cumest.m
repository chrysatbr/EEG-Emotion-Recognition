function y_cum = cumest(y,norder,maxlag,nsamp,overlap,flag,k1,k2)
%CUMEST	Second-, third- or fourth-order cumulants.
%	y_cum = cumest (y, norder, maxlag, samp_seg, overlap, flag, k1, k2)
%	       y - time-series  - should be a vector
%	  norder - cumulant order: 2, 3 or 4 [default = 2]
%	  maxlag - maximum cumulant lag to compute [default = 0]
%	samp_seg - samples per segment  [default = data_length]
%	 overlap - percentage overlap of segments [default = 0]
%	           overlap is clipped to the allowed range of [0,99].
%	   flag  - 'biased' or 'unbiased'  [default = 'biased']
%	  k1,k2  - specify the slice of 3rd or 4th order cumulants
%	  y_cum  - C2(m) or C3(m,k1) or C4(m,k1,k2),  -maxlag <= m <= maxlag
%	           depending upon the cumulant order selected

% Functions cumN_est (N=2,3,4) should be invoked
%  via this routine, since this routine does extensive parameter checking.

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.3 $
%  A. Swami, January 20, 1993;  August 6, 1994.

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

% ----------- Parameter checks ------------------

 [ksamp,nrecs] = size(y);
 if (ksamp == 1)
     ksamp = nrecs; nrecs = 1;
 end
 if (exist('norder') ~= 1) norder = 2; end
 if (norder < 2 | norder > 4)
    error('cumulant order must be 2, 3 or 4')
 end
 if (exist('maxlag') ~= 1) maxlag = 0; end
 if (maxlag < 0)
    error ('"maxlag" must be non-negative ')
 end

 if (exist('nsamp') ~= 1) nsamp = ksamp; end
 if (nrecs > 1)           nsamp = ksamp; end
 if (nsamp <= 0 | nsamp > ksamp)
    nsamp = ksamp;
 end

 if (exist('overlap') ~= 1) overlap = 0; end
 if (nrecs > 1)             overlap = 0; end
 overlap = max(0,min(overlap,99));

 if (exist('flag') ~= 1) flag = 'biased'; end
 if (flag(1:1) ~= 'b' & flag(1:1) ~= 'B')   flag = 'unbiased';
 else, flag = 'biased';
 end

 if (exist('k1') ~= 1) k1 = 0; end
 if (exist('k2') ~= 1) k2 = 0; end

% ---------- go, estimate the cumulants ------------------------
 if (norder == 2)
    y_cum = cum2est (y, maxlag, nsamp, overlap, flag);

 elseif (norder == 3)
    y_cum = cum3est (y, maxlag, nsamp, overlap, flag, k1);

 elseif (norder == 4)
    y_cum = cum4est (y, maxlag, nsamp, overlap, flag, k1, k2);

 end

return
