function [arvec, fpe, wt] = rivtr(y,morder,arorder, lambda, delta, nsmuth)
%RIVTR	 Recursive instrumental algorithm using the transversal structure.
%	[arvec, fpe, wt] = rivtr(y,morder,arorder,lambda,delta,nsmuth)
%	      y - time series
%	 morder - cumulant order          [default = 4]
%	arorder - AR order                [default = 2]
%	 lambda - forgetting factor       [default = 0.998]
%	 delta  - initialization constant [sign(R_zy(0))]
%	          where z = ivcal(y,morder,1), & R_zy(0) = E ( y(n) z(n) )
%	 nsmuth - smoothing window for AR estimation
%	          The default value is  min(nsamp/4,50), where
%	          nsamp is the length of the time-series y.
%	  arvec - steady-state AR parameter vector
%	   fpe  - prediction error
%	    wt  - weight matrix (as a function of time)
%	          each row corresponds to a time sample

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

% ---------------- Parameter Checks ------------------

 y=y(:); nsamp = length(y);
 if (exist('morder') ~= 1) morder = 4; end
 if (morder ~= 2 & morder ~= 3 & morder ~= 4)
    error(' morder should be 2, 3 or 4')
 end
 if (exist('arorder') ~= 1) arorder = 2; end
 if (exist('lambda') ~= 1) lambda = 0.998; end
 if (lambda <= 0 | lambda > 1)
    error(' lambda should be in the range (0,1]')
 end
 if (exist('nsmuth') ~= 1) nsmuth = min(nsamp/4, 50); end

 z = ivcal(y, morder, lambda);

 if (exist('delta') ~= 1)
   delta = sign(mean(y'*ivcal(y,morder,1)));
%   disp([' Initialization value delta set to ', int2str(delta)])
 end

 sswt = zeros(arorder,1);
 pmat = eye(arorder) * delta;
 z = [zeros(arorder,1); z];
 y = [zeros(arorder,1); y];

 laminv = 1/lambda;

%%

 for n = arorder+1:nsamp+arorder            % loop over time
   pmat = laminv * pmat;
   v    = z(n-1:-1:n-arorder);              % extract current v vector
   u    = y(n-1:-1:n-arorder);              % extract current u vector
   kvec = pmat * v / (1 + u' * pmat * v);   % the gain vector
   alph = y(n) - sswt' * u;                 % a priori prediction error
   sswt = sswt + alph * kvec;               % updated weight vector
   pmat = pmat - kvec * (u' * pmat);        % updated matrix inverse
   fpe(n-arorder) = y(n) - sswt' *u;        % prediction error
   wt(n,:) = sswt';                         % save time-varying weights
 end
 wt = wt(arorder+1:arorder+nsamp, :);
 arvec = [1; -mean(wt(nsamp-nsmuth+1:nsamp,:))'];


return
