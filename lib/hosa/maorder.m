function qopt = maorder(x,qmin,qmax,pfa, flag)
%MAORDER MA order  determination
%	qopt = maorder(x, qmin, qmax,pfa, flag)
%	x    - time series (must be a vector)
%	qmin - minimum MA order    [default = 0]
%	qmax - maximum MA order    [default = 10]
%	pfa  - probability of false alarm  [default = 0.05]
%	flag - if non-zero (default is 1), various values are displayed
%	qopt - estimated MA order

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

% -------------------- Parameter checks --------------------------

if (exist('x') ~= 1)
   error(' vector x must be specified')
end

[m,n] = size(x);
if (min(m,n) ~= 1)
   error(' x must be a vector, not a matrix')
end
if (m==1)
   x = x(:); m=n; n=1;
end
nsamp = m;

if (exist('qmin') ~= 1) qmin = 0; end
if (exist('qmax') ~= 1) qmax = 10; end
if (qmax < qmin)  qmax = qmin + 1; end
if (exist('pfa') ~= 1) pfa = 0.05; end
if (exist('flag') ~= 1) flag = 1;  end

cvec = cumest (x, 3, qmax+1,nsamp,0,'unbiased');    % the c(m,0) slice
cvec = cvec(qmax+3:2*qmax+3);                        % c(m,0) m=1, ... ,qmax+1

var   = zeros(qmax+1,1);
ord   = zeros(qmax+1,1);
thres = zeros(qmax+1,1);

xsq    = x.^2;
z      = zeros(nsamp,1);
inveps = erfinv(1-pfa);

% ----------------------------------------------------------------

for q=qmin:qmax
    m=q+1;
    z   = xsq(1:nsamp-m) .* x(m+1:nsamp);         % z(n) = y^2(n) * y(n+m)
    rzz = cumest (z,2,2*qmax,nsamp-m,0,'unbiased'); % covariance of z(n)
    rzz = rzz + cvec(m)^2; 	                  % correlation of z(n)
    var(m) = sum(rzz) / (nsamp-m);                % variance of c(q+1,0)
    thres(m) = inveps * sqrt(2*var(m));           % (1-pfa) sigma bounds
     ord(m) = (abs(cvec(m)) < thres(m));          % less ? then MA(q)
end
q    = max(find(ord == 0));          % max m for which MA(m) hypothesis fails
ind  = find( ord(q+1:qmax+1) > 0);   % now find where hypothesis holds
qopt = min(ind) + q - 1;             % fix the index
if (isempty(qopt)) qopt = qmax+1; end

if (flag)
   disp(['   q      var(cqk)          cqk            thres      result '])
   for q=qmin:qmax
       m = q+1;
       fprintf(' %3g %15.5e %15.5e',q, var(m), cvec(m))
       fprintf(' %15.5e     %1g \n', thres(m), ord(m))
   end
   disp([' Estimated MA order is ', int2str(qopt) ])
end
return

