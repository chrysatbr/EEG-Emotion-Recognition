function bvec = maest(y,q, norder,samp_seg,overlap,flag)
%MAEST  MA parameter estimation via the GM-RCLS algorithm, with Tugnait's fix
%	bvec = maest (y, q, norder, samp_seg, overlap, flag)
%	      y  - time-series (vector or matrix)
%	      q  - MA order
%	  norder - cumulant-order to use  [default = 3]
%	samp_seg - samples per segment for cumulant estimation
%	           [default: length of y]
%	 overlap - percentage overlap of segments  [default = 0]
%	    flag - 'biased' or 'unbiased'          [default = 'biased']
%	    bvec - estimated MA parameter vector

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

% ------ parameter checks -----------------------------

 if (nargin < 2)
    error('insufficient number of parameters')
 end

 [nsamp, nrecs] = size(y);
 if (nsamp == 1) nsamp = nrecs; nrecs = 1; y = y.'; end

 if (q <= 0)        bvec=1; return, end
 if (exist('norder') ~= 1) norder = 3; end
 if (norder ~= 3 & norder ~=4)
    error('cumulant order must be 3 or 4')
 end
 if (exist('samp_seg') ~= 1) samp_seg = nsamp; end
 if (exist('overlap') ~= 1) overlap = 0; end
 overlap = max(0, min(overlap,99));
 if (exist('flag') ~= 1) flag = 'biased'; end
 if (nrecs > 1)  samp_seg = nsamp;  overlap = 0; end

 % ---------- estimate cumulants and covariances --------------

  c2 = cumest(y,2,q, samp_seg, overlap, flag);
  c2 = [c2; zeros(q,1)];

  cumd = cumest(y,norder,q,samp_seg,overlap,flag,0,0);
  cumq = cumest(y,norder,q,samp_seg,overlap,flag,q,0);

  cumd = cumd(2*q+1:-1:1);
  cumd = [cumd; zeros(q,1)];
  cumq(1:q) = zeros(q,1);

% ----------- The GM-RCLS algorithm --------------------

  cmat = toeplitz(cumd, [cumd(1),zeros(1,q)]);
  rmat = toeplitz(c2,   [c2(1),  zeros(1,q)]);
  amat0 = [cmat, -rmat(:,2:q+1)];
  rvec0 = c2;

% ----------- The Tugnait fix ----------------------------

  cumq = [cumq(2*q+1:-1:q+1); zeros(q,1)];
  cmat4 = toeplitz(cumq, [cumq(1),zeros(1,q)]);
  c3   = cumd(1:2*q+1);
  amat0 = [amat0, zeros(3*q+1,1); zeros(2*q+1,q+1), cmat4(:,2:q+1), -c3];
  rvec0 = [rvec0; -cmat4(:,1)];

% ------------ Get rid of R(0) terms ----------------------

  row_sel = [1:q, 2*q+2:3*q+1, 3*q+2:4*q+1, 4*q+3:5*q+2];
  amat0 = amat0(row_sel,:);
  rvec0 = rvec0(row_sel);

% ------------ Solve for MA parms -------------------------

  bvec  = amat0 \ rvec0;
  b1    = bvec(2:q+1)/bvec(1);
  b2    = bvec(q+2:2*q+1);
  if (norder == 3)
     if (all(b2 > 0) )
        b1 = sign(b1) .* sqrt(0.5*(b1 .^2 + b2));
     else
        disp('MAEST: alternative solution b1 used')
     end
  else
     if (sign(b2) == sign(b1))
        b1 = sign(b1) .* (abs(b1) + abs(b2) .^(1/3) ) /2;
     else
        disp('MAEST: alternative solution b1 used')
     end
  end

  bvec = [1; b1];

return
