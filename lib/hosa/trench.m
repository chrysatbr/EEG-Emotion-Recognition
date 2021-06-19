function [amat,cmat,pf,gamf,gamb] = trench(c,r)
%TRENCH	 Trench recursion (non-symmetric Toeplitz matrix)
%	[amat,cmat, pfe,gamf,gamb] = trench(c,r)
%	c     - first column of Toeplitz matrix, r(k), k=0,...,M
%	r     - first row    of Toeplitz matrix, r(-k), k=0,...,M
%	          defaults to c'.
%	amat  - estimated forward AR predictors, of orders 1 through M
%	        the k-th column contains the AR(k) coefficients
%	cmat  - estimated backward AR predictors, of orders 1 through M
%	        the k-th column contains the AR(k) coefficients
%	pfe   - final prediction error variance for orders 1 through M
%	        the k-th element contains the value for order k
%	gamf   - estimated forward reflection coefficients
%	gamb   - estimated backward reflection coefficients
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

% ---------------------------------------------------------
c = c(:);
if (exist('r') ~= 1) r = c'; end
if (c(1) ~= r(1))
   disp(['r(1) ~= c(1);  setting r(1) = c(1) ']);
   r(1) = c(1);
end
m = min(length(c), length(r)) - 1;
c = c(1:m+1); r = r(1:m+1);
r = r(:);

pf = zeros(m+1,1);
df = zeros(m,1);
db = zeros(m,1);
gamf = zeros(m,1);
gamb = zeros(m,1);

pf(1) = c(1);
df(1) = c(2);
db(1) = r(2);

avec = [1];
cvec = [1];

amat = zeros(m+1,m);
cmat = zeros(m+1,m);
for k=1:m
    gamf(k) = -df(k) / pf(k);
    gamb(k) = -db(k) / pf(k);
    if ( gamf(k)*gamb(k)  > 1)
       disp(['gamf * gamb > 1 at order',int2str(k)])
       disp(['stopping iterations'])
       break
    end
    pf(k+1)  = pf(k) * (1 - gamf(k)*gamb(k));
    aold  = avec;
    avec  = [avec; 0] + gamf(k) * [0; cvec];
    cvec  = [0; cvec] + gamb(k) * [aold ;0];
    if (k < m)
        df(k+1) = c(k+2:-1:2).' * avec;
        db(k+1) = r(2:1:k+2).'  * cvec;
    end
    amat(1:k+1,k)   = avec;
    cmat(m+1-k:m+1,k) = cvec;
end
pf = pf(2:m+1);

return
