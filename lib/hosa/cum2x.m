function   y_cum = cum2x (x,y, maxlag, nsamp, overlap, flag)
%CUM2X Cross-covariance
%	y_cum  = cum2x (x,y,maxlag, samp_seg, overlap, flag)
%	  x,y    - data vectors/matrices with identical dimensions
%	           if x,y are matrices, rather than vectors, columns are
%	           assumed to correspond to independent realizations,
%	           overlap is set to 0, and samp_seg to the row dimension.
%	  maxlag - maximum lag to be computed    [default = 0]
%	samp_seg - samples per segment  [default = data_length]
%	 overlap - percentage overlap of segments [default = 0]
%	           overlap is clipped to the allowed range of [0,99].
%	   flag  - 'biased', biased estimates are computed  [default]
%	          'unbiased', unbiased estimates are computed.
%	   y_cum - estimated cross-covariance
%	           E x^*(n)y(n+m),   -maxlag <= m <= maxlag


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

%---------------- Parameter checks ---------------------------

[lx, nrecs] = size(x);
if ([lx,nrecs] ~= size(y))
   error('x,y should have identical dimensions')
end

if (lx == 1)
   lx = nrecs; nrecs = 1;
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

overlap  = fix(overlap/100 * nsamp);
nadvance = nsamp - overlap;
if (nrecs == 1)
   nrecs  = fix( (lx - overlap)/nadvance ) ;
end

nlags = 2*maxlag+1;
zlag  = maxlag + 1;
y_cum = zeros(nlags,1);

if (flag(1) == 'b' | flag(1) == 'B')
    scale = ones(nlags,1)/nsamp;
else
    scale = [lx-maxlag:lx,lx-1:-1:lx-maxlag]';
    scale = ones(2*maxlag+1,1) ./ scale;
end

% ----------------------------------------------------------

ind = (1:nsamp)';
for k=1:nrecs
    xs = x(ind);  xs = xs(:) - mean(xs);
    ys = y(ind);  ys = ys(:) - mean(ys);

    y_cum(zlag)  =  y_cum(zlag) + xs' * ys;

    for m = 1:maxlag
        y_cum(zlag-m)  = y_cum(zlag-m) + xs([m+1:nsamp])' * ys([1:nsamp-m]) ;
        y_cum(zlag+m)  = y_cum(zlag+m) + xs([1:nsamp-m])' * ys([m+1:nsamp]);
    end

    ind = ind + nadvance;
end

y_cum = y_cum .* scale / nrecs;

return
