function z = ivcal(y, morder, lam)
%IVCAL	Instrumental variables
%	z = ivcal (y, morder, lam)
%	 y     - input time series (vector or matrix)
%	morder - cumulant order;  default = 3;
%	lam    - forgetting factor: range (0,1]; default = 1;
%
%	 z     - instrumental variable, computed as follows:
%	 morder       z(n)
%	 q <= 0      y(n+q)       (delay)
%	   1         sign (y(n)
%	   2         y(n)
%	   3         y^2(n)
%	   4         y^3(n) - 3 s(n) * y(n)
%	where s(n) is the estimated variance, at sample n.

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.4 $
%  A. Swami   April 15, 1993.

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


% ----- Parameter checks ----------------------

[nsamp,nrecs] = size(y);
if (exist('morder') ~= 1) morder = 3; end
if (exist('lam') ~= 1) lam = 1; end
if (lam <= 0 | lam > 1)
   error(' lambda should lie in range (0,1]')
end

if (nsamp == 1)
    y = y(:); nsamp = nrecs; nrecs = 1; flipit = 1;
else
    flipit = 0;
end

% ------- Compute IV ------------------------
z = zeros(nsamp, nrecs);
if     (morder <= 0) z(-morder+1:nsamp,:) =  y(1:nsamp+morder,:);
elseif (morder == 1) z = sign(y);
elseif (morder == 2) z = y;
elseif (morder == 3) z = y.^2;
elseif (morder == 4)
       if (lam == 1)
          c2 = cumsum(y.^2) ./ ([1:nsamp]' * ones(1,nrecs));
          z = y.^3 - 3 * y .* c2;
       else
          var = y * 0;
	  for k=1:nrecs,
	    var(:,k) = filter(1, [1,-lam], y(:,k).^2);
	  end
          z = y.^3 - 3 * (1-lam) * var .* y ;
       end
else error(' morder should be less than 5')
end

if (flipit) z = z.'; end
return
