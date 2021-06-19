function   y_cum = cum4x (w,x,y,z, maxlag, nsamp, overlap, flag, k1, k2)
%CUM4X 	Fourth-order cross-cumulants.
%       y_cum = cum4x (w,x,y,z, maxlag, samp_seg, overlap, flag, k1, k2)
%	 w,x,y,z  - data vectors/matrices with identical dimensions
%	           if w,x,y,z are matrices, rather than vectors, columns are
%	           assumed to correspond to independent realizations,
%	           overlap is set to 0, and samp_seg to the row dimension.
%	  maxlag - maximum lag to be computed    [default = 0]
%	samp_seg - samples per segment  [default = data_length]
%	 overlap - percentage overlap of segments [default = 0]
%	           overlap is clipped to the allowed range of [0,99].
%	   flag : 'biased', biased estimates are computed  [default]
%	          'unbiased', unbiased estimates are computed.
%	  k1,k2 : the fixed lags in C4(m,k1,k2); defaults to 0
%	   y_cum:  estimated fourth-order cross cumulant,
%          c4(t1,t2,t3) := cum( w^*(t), x(t+t1), y(t+t2), z^*(t+t3) )

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.4 $
%  A. Swami   January 20, 1995

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


% c4(t1,t2,t3) := cum( x^*(t), w(t+t1), y(t+t2), z^*(t+t3) )
%  cum(w,x,y,z) := E(wxyz) - E(wx)E(yz) - E(wy)E(xz) - E(wz)E(xy)
%  and, w,x,y,z are assumed to be zero-mean.


% ---- Parameter checks are done in CUMEST ----------------------
 [lx,nrecs]  = size(w);
 if ([lx,nrecs] ~= size(x) | [lx,nrecs] ~= size(y) | [lx,nrecs] ~= size(z))
      error('w,x,y,z should have identical dimensions')
end

if (lx == 1)
  lx = nrecs;  nrecs = 1;
end

if (exist('maxlag') ~= 1) maxlag = 0;         end
if (maxlag < 0)
    error ('"maxlag" must be non-negative ')
end

if (exist('nsamp') ~= 1)     nsamp = lx; end
if (nrecs > 1)               nsamp = lx; end
if (nsamp <= 0 | nsamp > lx) nsamp = lx; end

if (exist('overlap') ~= 1) overlap = 0; end
if (nrecs > 1)             overlap = 0; end
overlap = max(0,min(overlap,99));

if (exist('flag')   ~= 1)   flag = 'biased' ; end
if (flag(1:1) ~= 'b' & flag(1:1) ~= 'B')
        flag = 'unbiased';
else,   flag = 'biased';
end

overlap0 = overlap;
overlap  = fix(overlap/100 * nsamp);
nadvance = nsamp - overlap;
if (nrecs == 1)
   nrecs  = fix( (lx - overlap)/nadvance ) ;
end

if (exist('k1')     ~= 1)     k1 = 0;         end
if (exist('k2')     ~= 1)     k2 = 0;         end


% ------ scale factors for unbiased estimates --------------------

nlags = 2 * maxlag + 1;
zlag  = 1 + maxlag;

tmp   = zeros(nlags,1);
if (flag(1:1) == 'b'  | flag(1:1) == 'B')
     scale = ones(nlags,1) / nsamp;
     sc1 = 1/nsamp;
     sc2 = sc1; sc12 = sc1;
else
     ind   = [-maxlag:maxlag]';
     kmin  = min(0,min(k1,k2));
     kmax  = max(0,max(k1,k2));
     scale = nsamp - max(ind,kmax) + min(ind,kmin);
     scale = ones(nlags,1) ./ scale;
     sc1  = 1/(nsamp-abs(k1));
     sc2  = 1/(nsamp-abs(k2));
     sc12 = 1/(nsamp-abs(k1-k2));
end

% ----------- estimate second- and fourth-order moments; combine ------

 y_cum  = zeros(2*maxlag+1,1) ;
 rind = -maxlag:maxlag;
 ind   = 1:nsamp;

 for i=1:nrecs
     tmp = y_cum * 0;
     R_zy   = 0;    R_wy = 0;  M_wz = 0;
     ws = w(ind);  ws = ws(:) - mean(ws);
     xs = x(ind);  xs = xs(:) - mean(xs);
     ys = y(ind);  ys = ys(:) - mean(ys);  cys = conj(ys);
     zs = z(ind);  zs = zs(:) - mean(zs);

     ziv = xs * 0;
%                     create the "IV" matrix: offset for second lag

     if (k1 >= 0)
 	ziv(1:nsamp-k1)  = ws(1:nsamp-k1) .* cys(k1+1: nsamp);
	R_wy = R_wy + ws(1:nsamp-k1)' * ys(k1+1:nsamp);
     else
       	ziv(-k1+1:nsamp) = ws(-k1+1:nsamp)  .* cys(1:nsamp+k1);
	R_wy = R_wy + ws(-k1+1:nsamp)' * ys(1:nsamp+k1);
     end

%                     create the "IV" matrix: offset for third lag

     if (k2 >= 0)
        ziv(1:nsamp-k2) = ziv(1:nsamp-k2) .* zs(k2+1: nsamp);
        ziv(nsamp-k2+1:nsamp) = zeros(k2,1);
	M_wz = M_wz + ws(1:nsamp-k2).' * zs(k2+1:nsamp);
     else
        ziv(-k2+1:nsamp) = ziv(-k2+1:nsamp) .* zs(1:nsamp+k2);
        ziv(1:-k2)    = zeros(-k2,1);
	M_wz = M_wz + ws(-k2+1:nsamp).' * zs(1:nsamp-k2);
     end

     if (k1-k2 >= 0)
       R_zy = R_zy + zs(1:nsamp-k1+k2)' * ys(k1-k2+1:nsamp);
     else
       R_zy = R_zy + zs(-k1+k2+1:nsamp)' * ys(1:nsamp-k2+k1);
     end

     tmp(zlag)  =  tmp(zlag) + ziv' * xs ;
     for k = 1:maxlag
         tmp(zlag-k) = tmp(zlag-k) + ziv([k+1:nsamp])' * xs([1:nsamp-k]);
         tmp(zlag+k) = tmp(zlag+k) + ziv([1:nsamp-k])' * xs([k+1:nsamp]);
     end

     y_cum = y_cum + tmp .* scale ;     % fourth-order moment estimates done

     R_wx = cum2x(ws,      xs, maxlag,         nsamp, overlap0, flag);
     R_zx = cum2x(zs,      xs, maxlag+abs(k2), nsamp, overlap0, flag);
     M_yx = cum2x(cys,     xs, maxlag+abs(k1), nsamp, overlap0, flag);


     y_cum  = y_cum - R_zy * R_wx * sc12    ...
	       - R_wy * R_zx(rind - k2 + maxlag+abs(k2)+ 1) * sc1 ...
	       - M_wz'* M_yx(rind - k1 + maxlag+abs(k1)+ 1) * sc2 ;

     ind = ind + nadvance ;
end

y_cum = y_cum / nrecs;
return


% mods Sept 7, 1998
% fixed conjugation errors in computing R_wy, R_zy, M_yx 

