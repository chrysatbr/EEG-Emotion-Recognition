function [delay,avec] = tde(x,y, max_delay, nsamp,svdflag)
%TDE	Time-delay estimation using cross-cumulant method.
%	[delay,avec] = tde(x,y, max_delay,nsamp,svdflag)
%	x        - data at sensor 1
%	y        - data at sensor 2
%	max_delay - maximum delay
%	nsamp     - samples/record for computing cumulants
%	svdflag   - if its value is non-zero, the SVD will be computed and
%	            you will be asked to choose an order for the low-rank
%	            approximation;  default value is 0
%	delay     - estimated delay (positive means that y lags x)
%	avec      - the "AR vector" in the parametric method


%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.9 $
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

%------------ Parameter Checks -----------------------------
  [m1,n1] = size(x);
  [lx,nrecs] = size(y);
  if (m1 ~= lx | n1 ~= nrecs)
     error('matrices x and y should have the same dimensions')
  end
  if (lx == 1), lx = nrecs; nrecs = 1; x = x(:); y = y(:); end
  if (nrecs > 1) nsamp = lx; end
  if (exist('nsamp') ~= 1) nsamp = lx; end
  if (nsamp <= 0 | nsamp > lx) nsamp = lx; end
  if (nrecs == 1)
     nrecs = fix(lx/nsamp); lx = nsamp * nrecs;
  end
  if (exist('svdflag') ~= 1) svdflag = 0; end

% -------- compute the cross-cumulants ---------------

  ind = 1:nsamp;  ind1=2:nsamp;  ind2=1:nsamp-1;
  zn = zeros(1,nsamp); zp = zn; z = zn;
  MD = max_delay;
  nlag1 = max_delay + 1;
  nlag2 = 2 * max_delay + 1;
  Rxyx = zeros(2*max_delay+1,3);
  Rxxx = zeros(4*max_delay+1,3);

  jind2 = nsamp-max_delay:nsamp+max_delay;
  jind4 = nsamp-2*max_delay:nsamp+2*max_delay;

  for k=1:nrecs
    x1 = x(ind) - mean(x(ind));
    y1 = y(ind) - mean(y(ind));
    z  = x1 .* x1;                                 % x(n) x(n)
    zn(2:nsamp) = x1(ind1) .* x1(ind2);             % x(n)x(n-1)
    zp(1:nsamp-1) = zn(2:nsamp);                    % x(n)x(n+1)

    t1 = conv(zn',flipud(y1));
    t2 = conv(z,  flipud(y1));
    t3 = conv(zp',flipud(y1));

    Rxyx(:,1) = Rxyx(:,1) + t1(jind2);
    Rxyx(:,2) = Rxyx(:,2) + t2(jind2);
    Rxyx(:,3) = Rxyx(:,3) + t3(jind2);

    t1 = conv(zn',flipud(x1));
    t2 = conv(z,  flipud(x1));
    t3 = conv(zp',flipud(x1));

    Rxxx(:,1) = Rxxx(:,1) + t1(jind4);
    Rxxx(:,2) = Rxxx(:,2) + t2(jind4);
    Rxxx(:,3) = Rxxx(:,3) + t3(jind4);

   ind = ind + nsamp;
end
clear t1 t2 t3 jind2 jind4

Rxxx  = Rxxx/length(x(:));
Rxyx  = Rxyx/length(x(:));

% -----  set up the system of linear equations -------------------

 Amat = zeros(2*MD+1,2*MD+1);  rvec = zeros(2*MD+1,1);
    for ksl = 1:3
      tmp = toeplitz(Rxxx(nlag2:4*MD+1,ksl), Rxxx(nlag2:-1:1,ksl));
      Amat = Amat + tmp' * tmp;
      rvec = rvec + tmp' * Rxyx(:,ksl);
    end

% -----  SVD low-rank approximation (optional) ---------------------

set(gcf,'Name','Hosa TDE')
if (svdflag)
   clf
   [umat,smat,vmat] = svd(Amat);
   s = diag(smat); ns =length(s); mord = ns - 2;

   clf,subplot(211),plot(s),title('singular values'), grid on 
   subplot(212)                 % for next plot
   txt1 = ['Order for SVD low-rank approximation ---> [',int2str(mord),']'];
   mord = input(txt1);
   if (isempty(mord)), mord = ns -2; end
   if (mord > ns | mord <= 0) mord = ns-2; end
   if (mord == 1) smat = 1/s(1);
      
   % else,  smat = diag(ones(s(1:mord)) ./ s(1:mord));
   %  Modified by A. Swami - September 11, 1998.
   else, smat = diag(ones(mord,1) ./ s(1:mord))
   end
   avec = umat(:,1:mord) * smat * vmat(:,1:mord)' * rvec;
else
    avec = Amat \ rvec;
end
    avec = avec(2*MD+1:-1:1);

[val, delay] = max(abs(avec));
delay = -MD - 1 + delay;

disp(['Estimated delay= ',num2str(delay)])

plot(-max_delay:max_delay,avec, -max_delay:max_delay,avec,'o'), grid on 
title(['TDE: parameter vector, delay = ' int2str(delay)])
xlabel('Delay in samples')

return
