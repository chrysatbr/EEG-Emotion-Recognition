function [arvec, fref, bref,fpe] =  ...
            rivdl (y,morder,arorder,lam,delta,thres, nsmuth)
%RIVDL	Recursive instrumental variable algorithm using the double lattice
%	[arvec, fref, bref, fpe] = rivdl(y,morder,p,lam,delta,thres,nsmuth)
%	      y - data vector
%	 morder - cumulant order (2, 3 or 4;  default = 4)
%	      p - number of stages                          (default = 2)
%	    lam - forgetting factor:  0 < lambda =< 1;      (default = 0.998 )
%	  delta - initialization value for F(n=0) and B(n=0) (default = 0.01)
%	  thres - threshold check for division by "zero"    (default = 0.0001)
%	 nsmuth - smoothing window for AR estimation
%	          the default value is  min(nsamp/4,50),
%	          where nsamp is the length of the time-series y.
%	 arvec - AR parameters corresponding to the smoothed final reflection
%	         coefficient estimates
%	fref, bref -  forward and backward reflection coefficients of the top
%	              lattice; rows correspond to time;
%	              columns s correspond to stages.
%	 fpe   - final prediction error (upper leg of lattice)

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.4 $
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

% ------- Parameter checks ---------------------------
    [nsamp, nrecs] = size(y);
    if (nsamp == 1) y = y.'; nsamp = nrecs; nrecs = 1; end
    if (nrecs > 1)
       error('y must be a vector, not a matrix')
    end

    if (exist('morder') ~= 1) morder = 4; end
    if (morder ~= 2 & morder ~= 3 & morder ~= 4)
       error(' morder should be 2, 3 or 4')
    end
    if (exist('lam') ~= 1) lam = 0.998; end
    if (lam <= 0 | lam > 1)
       error(' lambda should be in the range (0,1]')
    end
    if (exist('delta') ~= 1) delta = 0.01; end
    if (exist('thres') ~= 1) thres = 0.0001; end
    if (exist('arorder')~= 1) arorder = 2; end
    if (exist('nsmuth') ~= 1) nsmuth = min(nsamp/4, 50); end

    z = ivcal (y, morder, lam);      % the instrumental variable

%  Initialize quantities at time 0

    delf_o = zeros(arorder,1);                  % (5.83 a)
    delb_o = delf_o;                            % (5.83 a)
    fcor   = delta * ones(arorder+1,1);         % (5.83 b)
    bcor_o = fcor;                              % (5.83 b)

    b_o    = zeros(arorder+1,1);                % all internal variables are 0
    bt_o   = b_o;                               % at time 0

% Loop over time
    for n=1:nsamp

        f(1) = y(n);  b(1) = y(n);                % (5.83 c)
        ft(1) = z(n);  bt(1) = z(n);              % (5.83 d)
        fcor(1) = lam * fcor(1) + y(n) * z(n);    % (5.83 e)
        bcor(1) = fcor(1);                        % (5.83 e)
        gam(1) = 1;                               % (5.83 f)

        for m=1:arorder                           % loop over order

            gfac = 1;  bfac = 1;  ffac = 1;       % divide-by-zero checks
            if (abs(gam(m))    > thres) gfac = 1/gam(m);    end
            if (abs(bcor_o(m)) > thres) bfac = 1/bcor_o(m); end
            if (abs(fcor(m))   > thres) ffac = 1/fcor(m);   end

            delf(m) = lam * delf_o(m) + f(m) * bt_o(m) * gfac;      % (5.84 a)
            delb(m) = lam * delb_o(m) + ft(m) * b_o(m) * gfac;      % (5.84 b)

            fref(n,m) = - delf(m) * bfac;                           % (5.85 a)
            bref(n,m) = - delb(m) * ffac;                           % (5.85 b)

            ftref     = - delb(m) * bfac;                          % (5.85 c')
            btref     = - delf(m) * ffac;                          % (5.85 c")


            f(m+1)    = f(m) + fref(n,m) * b_o(m);                 % (5.85 d)
            b(m+1)    = b_o(m) + bref(n,m) * f(m);                 % (5.85 e)

            ft(m+1)  = ft(m) + ftref * bt_o(m);                    % (5.85 f)
            bt(m+1)  = bt_o(m) + btref * ft(m);                    % (5.85 g)

            fcor(m+1) = fcor(m) + fref(n,m) * delb(m);             % (5.85 h)
            bcor(m+1) = bcor_o(m) + bref(n,m) * delf(m);           % (5.85 i)

            gam(m+1)  = gam(m) - bt_o(m) * b_o(m) * bfac;          % (5.85 j)
       end

       delf_o = delf; delb_o = delb;
       bcor_o = bcor;
       b_o    = b;  bt_o = bt;
       fpe(n) = f(arorder+1);

     end

% convert to AR coefficients

     fr = mean(fref(nsamp-nsmuth:nsamp,:));
     br = mean(bref(nsamp-nsmuth:nsamp,:));


     a = zeros(arorder+1, arorder+1);
     a(:,1) = ones(arorder+1,1);
     c = diag(ones(arorder+1,1));

     for m = 1: arorder
         a(m+1,m+1) = fr(m);
         c(m+1,1)   = br(m);
         for k = 2: m
             a(m+1,k) = a(m,k) + fr(m) * c(m,k-1);
             c(m+1,k) = c(m,k-1)   + br(m) * a(m,k);
         end
     end

    arvec = a(m+1,:)';

return
