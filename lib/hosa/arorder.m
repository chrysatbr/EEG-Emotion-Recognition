function  p = arorder(y,norder, pmax,qmax, flag)
%ARORDER  estimates AR order
%	p   = arorder (y, norder, pmax, qmax, flag)
%	y      - time-series
%	norder - cumulant order to use, should be 2,3,4,-3 or -4 [default = 3]
%	         a value of -3 (-4) indicates that both correlation
%	         and third- (fourth-) order cumulants should be used
%	pmax   - maximum AR order  [default = 10]
%	qmax   - maximum MA order  [default = 10]
%	flag   - if 1, the internally chosen AR order is returned in 'p',
%	         otherwise, the plot of singular values is displayed,
%	         and the user is prompted to choose the order.
%	p      - estimated AR order

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.6 $
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

%---------- Parameter checks ---------------------------------------

if (exist('y') ~= 1)
   error(' vector y must be specified')
end

[m,n] = size(y);
if (min(m,n) ~= 1)
   error(' x must be a vector, not a matrix')
end
if (m==1)
   y = y(:); m=n; n=1;
end
nsamp = m;

if (exist('qmax') ~= 1) qmax = 10; end
if (exist('pmax') ~= 1) pmax = 10; end

if (exist('norder') ~= 1) norder = 3; end
if (norder ~= 2 & abs(norder) ~= 3 & abs(norder) ~= 4)
   error(' norder must be 2, 3, -3, 4 or -4')
end

if (exist('flag') ~= 1) flag = 1; end

maxlag = qmax + pmax;
minlag = -maxlag;
nlags  = maxlag - minlag + 1;

% cumulant estimates --------------------------------------------------

  if (norder ~= 2)
    kslice1 = [-pmax, qmax];
    kslice2 = [0,0];
    kslice = (kslice1(2) - kslice1(1) + 1) * (kslice2(2) - kslice2(1) + 1);
    cum_y = zeros(nlags,kslice);

    morder = abs(norder);
    kloc = 0;
    for k1 = kslice1(1) : kslice1(2)
       for k2 = kslice2(1) : kslice2(2)
           kloc = kloc + 1;
           cum_y(:,kloc) = cumest(y, morder, maxlag, nsamp, 0,  ...
                                  'biased', k1, k2);
       end
     end

% ----------------- set up the Hankel matrix ------------
       q = qmax; p = pmax;
       Amat = hankel(cum_y(q-p+1-minlag+1:nlags-p, 1), ...
                  cum_y(nlags-p:nlags-1,1) );
       rvec = cum_y(q+1-minlag+1:nlags, 1);
       for k=2:kslice
           Amat = [Amat; hankel(cum_y(q-p+1-minlag+1:nlags-p, k), ...
                          cum_y(nlags-p:nlags-1,k) ) ];
           rvec = [rvec; cum_y(q+1-minlag+1:nlags, k) ];
       end
       rvec = -rvec;
   end
%%
   if (norder == 2 | norder < 0)
       q = qmax; p = pmax;
       cor_y = cumest(y,2,maxlag,nsamp);
       AR = hankel(cor_y(q-p+1-minlag+1:nlags-p, 1), ...
                  cor_y(nlags-p:nlags-1,1) );
       br = -cor_y(q+1-minlag+1:nlags, 1);

       Amat = [Amat; AR];  rvec = [rvec; br];
   end

%%-------------- compute svd -----------------------

   [U,S,V] = svd(Amat);
   sd = diag(S); sd = sd/sd(1);
   sdd = -diff(sd);
   [rmax,popt] = max(sdd);
   plot(sdd); grid on 
   title(['difference in singular values, cumulant order=',int2str(norder)])
   set(gcf, 'Name', 'Hosa ARORDER')
   if (flag == 1)
     p = popt;
     return
   end
   disp(['The singular values of the cumulant matrix are displayed on',...
          'the graphics window'])
%   plot(sd); grid, title('singular values')
%   ss = cumsum(sd) / sum(sd);
%   plot(ss); grid, title('cumulative singular values')

   disp(['The estimated AR order is ', int2str(popt)])
   rtxt = ['choose AR order (0 to ',int2str(pmax),' ) ---> [', ...
   		int2str(popt),'] '];
   p = -1;
   while (p < 0 | p > pmax)
      p = input(rtxt);
      if (isempty(p)) p = popt; end
   end
   disp(' ')
