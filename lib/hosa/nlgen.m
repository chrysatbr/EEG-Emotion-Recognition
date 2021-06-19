function y = nlgen (x, h, q)
%NLGEN	generates the output of a second order Volterra system
%	y = nlgen (x, h, q)
%	y(n) = sum_{k} h(k) x(n-k) + sum_{k} sum_{l} q(k,l) x(n-k)x(n-l)
%            the filters are assumed to be causal
%	 x is the input time series;  if x is a matrix, columns are assumed
%	   to correspond to different realizations.
%	h is the linear part of the transfer function; must be real
%	q is the quadratic part of the transfer function; must be real
%         and symmetric;  the origin is at (1,1) with  axes pointing
%         downwards and to the right
%	y is the output of the Volterra system (same dimensions as x).
%

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.5 $
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

% Parameter Checks ------------------------------------------------
[mrow,ncol]  = size(q);

if (mrow ~= ncol)
   error ('q must be a square matrix')
end

if (q ~= q.')
   error(' q must be real and symmetric')
end

if (any(imag(h) ~= 0))
   error (' h must be real ')
end

if (any(any(imag(x) ~= 0)))
   error('x must be a real-valued time-series (vector/matrix)')
end

[nsamp, nreal] = size(x);

if (nsamp == 1)
   trflag = 1;
   x = x.';
   nsamp = nreal;
   nreal = 1;
else
   trflag = 0;
end

% ---------- Generate --------
y = zeros(nsamp,nreal);

for m=1:nreal
    xm = x(:,m);
    ym = filter (h, 1, xm);

    for k=1:mrow
        yk = filter (q(k,:), 1, xm);   %  y(n,k) = sum_{l} q(k,l) x(n-l)
        x1 = [zeros(k-1,1); xm(1:nsamp-k+1)];
        ym = ym + yk .* x1(1:nsamp);     % accumulate for each k,
                                         %   y(n,k)x(n-k)
    end
    y(:,m) = ym;
end

if (trflag)
   x = x.';
end

return
