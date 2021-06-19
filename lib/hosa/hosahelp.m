function hosahelp
%HOSAHELP Higher-Order Spectral Analysis Toolbox 
%         -- Version 2.0.3 (R12 compliant)  27 Dec 2000
%	  Quick one line help
%           [avec, bvec] = armaqs   (y,p,q,norder,maxlag,sampSeg,overlap,flag)
%           [avec, bvec] = armarts  (y,p,q,norder,maxlag,sampSeg,overlap,flag)
%                   zmat = armasyn  (default)
%                      p = arorder  (y,norder, pmax,qmax, flag)
%                   avec = arrcest  (y,p,q,norder,maxlag,sampSeg,overlap,flag)
%[he,ceps,A,B,minh,maxh] = biceps   (y,p,q,nsamp,overlap,flag,lh)
%            [hest,ceps] = bicepsf  (y,nlag,nsamp, overlap,flag, nfft, wind)
%            [bic,waxis] = bicoher  (y,  nfft, wind, nsamp, overlap)
%	     [bic,waxis] = bicoherx (w,x,y,  nfft, wind, nsamp, overlap)
%          [Bspec,waxis] = bispecd  (y,  nfft, wind, nsamp, overlap)
%          [Bspec,waxis] = bispecdx (x, y, z,nfft,wind,nsamp,overlap,plotflag)
%          [Bspec,waxis] = bispeci  (y,nlag,nsamp, overlap,flag, nfft, wind)
%          [Bspec,waxis] = bispect  (ma, ar, nfft)
%                  y_cum = cum2est  (y, maxlag, nsamp, overlap, flag)
%                  y_cum = cum3est  (y, maxlag, nsamp, overlap, flag, k1)
%                  y_cum = cum4est  (y, maxlag, nsamp, overlap, flag, k1, k2)
%                  y_cum = cum2x    (x,y, maxlag, nsamp, overlap, flag)
%                  y_cum = cum3x    (x,y,z, maxlag, nsamp, overlap, flag, k1)
%                  y_cum = cum4x    (w,x,y,z,maxlag,nsamp,overlap,flag,k1,k2)
%                  y_cum = cumest   (y,norder,maxlag,nsamp,overlap,flag,k1,k2)
%                   cmat = cumtrue  (ma, ar, norder, nlags, k)
%   [spec,theta,bearing] = doa      (ymat, dspace, dtheta,nsource,order,delta)
%                   smat = doagen   (default)
%               [sg, sl] = glstat   (x,cparm,nfft)
%          [Pxx,ar1,ar2] = harmest  (y,maxlag,p_order,flag,nfft,norder)
%                   zmat = harmgen  (default)
%	                   hosademo (figure_number)
%			   hosahelp
%     [a,theta,alpha,fr] = hprony   (x,p)
%                      z = ivcal    (y, morder, lam)
%                   bvec = maest    (y,q, norder,samp_seg,overlap,flag)
%                   qopt = maorder  (x,qmin,qmax,pfa, flag)
%                   hest = matul    (bisp)
%                      y = nlgen    (x, h, q)
%                  [h,q] = nlpow    (x,y,nfft)
%                [ht,qt] = nltick   (x,y,nfft, wind,segsamp,overlap)
%	       [loc,val] = pickpeak (spec, npicks, rdiff)
%                   zdat = qpcgen   (default)
%         [ar_vec,Bspec] = qpctor   (y,maxlag,arOrder,nfft,nsamp,overlap,flag)
%[arvec, fref, bref,fpe] = rivdl    (y,morder,arorder,lam,delta,thres, nsmuth)
%       [arvec, fpe, wt] = rivtr    (y,morder,arorder, lambda, delta, nsmuth)
%                      u = rpiid    (nsamp, in_type,p_spike)
%           [delay,avec] = tde      (x,y,max_delay,  nsamp,svdflag)
%           [delay,ctau] = tdeb     (x,y,max_delay,  nfft,wind,nsamp,overlap)
%               [s1, s2] = tdegen   (default)
%            [delay,rxy] = tder     (x,y,max_delay,  segsamp,overlap,nfft)
%               [x,flag] = tls      (A,b)
%   [amat,cmat,pf,gf,gb] = trench   (c,r)
%          [Tspec,waxis] = trispect (ma, ar, nfft, f3)
%            [wx, waxis] = wig2     (x0,nfft,		flag)
%            [wx, waxis] = wig2c    (x0,nfft,sigma,	flag)
%            [wx, waxis] = wig3     (x0,nfft,		flag)
%            [wx, waxis] = wig3c    (x0,nfft,sigma,	flag)
%            [wx, waxis] = wig4     (x0,nfft,		flag)
%            [wx, waxis] = wig4c    (x0,nfft,sigma,	flag)
%

help hosahelp

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
