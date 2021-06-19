function [avec, bvec] = armaqs(y,p,q, norder,maxlag,samp_seg,overlap,flag)
%ARMAQS	Estimates ARMA parameters via the q-slice algorithm.
%	[avec, bvec] = armaqs(y,p,q, norder,maxlag,samp_seg,overlap,flag)
%	      y : time-series (vector or matrix)
%	      p : AR order
%	      q : MA order
%	  norder: cumulant order:  3 or 4         [default = 3 ]
%	  maxlag: maximum cumulant lag to be used [default = p + q]
%	samp_seg: samples per segment for estimating cumulants
%	                                          [default = length of y]
%	overlap : percentage overlap of segments  [default = 0]
%	   flag : 'biased' or 'unbiased'          [default = 'biased']
%	   avec : estimated AR parameter vector
%	   bvec : estimated MA parameter vector

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc.
%       $Revision: 1.4 $
%  A. Swami   January 20, 1993

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

% ------------------- parameter checks ----------------

 if (nargin < 3)
    error('insufficient number of parameters')
 end

 [nsamp, nrecs] = size(y);
 if (nsamp == 1) nsamp = nrecs; nrecs = 1; y = y.';  end

 if (p < 0)
    error('AR order cannot be negative')
 end
 if (p == 0)
    error('please use MAEST for the pure MA (p=0) case')
 end
 if (q < 0)
    error('MA order cannot be negative')
 end

 if ~exist('norder') norder = 3; end
 if (norder ~= 3 & norder ~= 4)
    error('norder must be 3, or 4')
 end

 maxlag0 = q + p;
 if (exist('maxlag') ~= 1) maxlag = maxlag0; end
 if (maxlag < maxlag0)
    disp(['ARMAQS: maxlag changed from ',int2str(maxlag), ...
            ' to ',int2str(maxlag0)])
    maxlag = maxlag0;
 end

 if (exist('samp_seg') ~= 1)  samp_seg = nsamp; end
 if (exist('overlap') ~=1)     overlap = 0;     end
 overlap = max(0, min(99, overlap) );
 if (exist('flag')  ~= 1)      flag = 'biased'; end

 if (nrecs > 1)  overlap = 0; samp_seg = nsamp; end

%-----------------------------------------------------------

% simultaneous AR and IR:
% first, the IR part:
% the q-slice IR equations are of the form
%     [I Ac][eh(0), ... ,eh(q), a(p), ... ,a(1)]' = -[bc; b];
% hence, concatenating,
%     |I Ac| |eh| = -|bc|
%     |0 A | |a | =  |b |

  ma_order = q;
  ar_order = p;
  zlag = max([p, p-q]);

  zlag1 = 1 + zlag + ma_order - ar_order;
  kloc = 0; k2 = 0;  cum_y = zeros(ar_order+1,ma_order+1);
  for k1 = 0:ma_order
       kloc = kloc + 1;
       alpha = cumest(y,norder,zlag,samp_seg,overlap,flag,k1,k2);
       cum_y(:,kloc) = alpha(zlag1:zlag1+ar_order);
  end
  cum_y = cum_y.';
  Acs   = cum_y(:,1:ar_order);  bcs = cum_y(:,ar_order+1);


%------------- Now for the AR part -----------------------------

  if ~rem(maxlag,2), maxlag = maxlag + 1; end
  nlags = 2*maxlag + 1;
  cum_y = zeros(nlags,p+1);
  for k = 1:p+1,
      cum_y(:,k) = cumest(y,norder,maxlag,samp_seg,overlap,flag,k+q-p-1,0);
  end
  pmax = (maxlag+1)/2;
  zlag = maxlag;
  zmax = zlag + maxlag - pmax + 1;
  AS = hankel(cum_y(zlag+1:zmax,1), ...
              cum_y(zmax:zlag+maxlag,1));
  bs = cum_y(zlag+pmax+1:zlag+maxlag+1,1);
  for k=2:p+1
      AS = [AS; hankel(cum_y(zlag+1:zmax,k), ...
                       cum_y(zmax:zlag+maxlag,k))];
      bs = [bs; cum_y(zlag+pmax+1:zlag+maxlag+1,k)];
  end


  [U,S,V] = svd([AS,bs]); V = V';
  Ahat = U(:,1:ar_order) * S(1:ar_order,1:ar_order) * V(1:ar_order,:);
  AS = Ahat(:,1:ar_order); bs=Ahat(:,ar_order+1);
  for k=2:p+1-ar_order
      AS = [AS; Ahat(:,k:k+ar_order-1)];
      bs = [bs; Ahat(:,k+ar_order)];
  end
  bs = -bs;

% --------------- Simultaneous TLS solution: -------------

  [mrow,ncol] = size(AS);
  Asvd = [ eye(q+1), Acs; zeros(mrow,q+1), AS];
  bsvd =  [ -bcs; bs];

  arma_order = ar_order + ma_order + 1;

  arma_vec = tls(Asvd,bsvd);
  avec = [1; arma_vec(arma_order:-1:ma_order+2)];
  hvec = arma_vec(1:ma_order+1);
  hvec = hvec/hvec(1);
  bvec = filter(avec,[1],hvec);

return
