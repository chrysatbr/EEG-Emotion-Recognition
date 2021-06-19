function [bic,waxis] = bicoherx (w,x,y,  nfft, wind, nsamp, overlap)
%BICOHERX - Direct (FD) method for estimating cross-bicoherence
%	[bic,waxis] = bicoherx (w,x,y,  nfft, wind, segsamp, overlap)
%	w,x,y - data vector or time-series
%	      - should have identical dimensions
%	nfft - fft length [default = power of two > nsamp]
%	       actual size used is power of two greater than 'nsamp'
%	wind - specifies the time-domain window to be applied to each
%	       data segment; should be of length 'segsamp' (see below);
%		otherwise, the default Hanning window is used.
%	segsamp - samples per segment [default: such that we have 8 segments]
%	        - if x is a matrix, segsamp is set to the number of rows
%	overlap - percentage overlap, 0 to 99  [default = 50]
%	        - if y is a matrix, overlap is set to 0.
%	bic     - estimated cross-bicoherence: an nfft x nfft array, with
%	          origin at center, and axes pointing down and to the right.
%	waxis   - vector of frequencies associated with the rows and columns
%	          of bic;  sampling frequency is assumed to be 1.

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc.
%       $Revision: 1.7 $
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


% --------------------- parameter checks -----------------------------

    if (size(w) ~= size(x) | size(x) ~= size(y) )
       error(' w, x, and y should have identical dimensions')
    end
    [ly, nrecs] = size(y);
    if (ly == 1)
    	ly = nrecs; nrecs = 1;
         w = w(:);  x = x(:);  y = y(:);
    end

    if (exist('nfft') ~= 1)            nfft = 128; end
    if (exist('overlap') ~= 1)      overlap = 50;  end
    overlap = max(0,min(overlap,99));
    if (nrecs > 1)                  overlap = 0;   end
    if (exist('nsamp') ~= 1)          nsamp = 0;  end
    if (nrecs > 1)                    nsamp = ly;  end
    if (nrecs == 1 & nsamp <= 0)
       nsamp = fix(ly/ (8 - 7 * overlap/100));
    end
    if (nfft  < nsamp)   nfft = 2^nextpow2(nsamp); end

    overlap  = fix(overlap/100  * nsamp);
    nadvance = nsamp - overlap;
    nrecs    = fix ( (ly*nrecs - overlap) / nadvance);

% ----------------------------------------------------------------------
    if (exist('wind') ~= 1) wind = hanning(nsamp); end
    [rw,cw] = size(wind);
    if (min(rw,cw) ~= 1 | max(rw,cw) ~= nsamp)
	   disp(['Segment size  is ',int2str(nsamp)])
	   disp(['"wind" array  is ',int2str(rw),' by ',int2str(cw)])
	   disp(['Using default window'])
	   wind = hanning(nsamp);
    end
    wind = wind(:);
% ---------------- accumulate triple products ----------------------

    bic  = zeros(nfft,nfft);
    Pyy  = zeros(nfft,1);         Pww = Pyy; Pxx = Pyy;

    mask = hankel([1:nfft],[nfft,1:nfft-1] );   % the hankel mask (faster)
    Yf12 = zeros(nfft,nfft);
    ind  = [1:nsamp]';

    for k = 1:nrecs
        ws  = w(ind); ws = (ws-mean(ws)).* wind;
        Wf  = fft(ws,nfft)  / nsamp;  CWf = conj(Wf);
        Pww = Pww + Wf .* CWf;

        xs  = x(ind); xs = (xs-mean(xs)).* wind;
        Xf  = fft(xs,nfft)  / nsamp;  CXf = conj(Xf);
        Pxx = Pxx + Xf .* CXf;

        ys  = y(ind);  	ys = (ys(:) - mean(ys)) .* wind;
	Yf  = fft(ys,nfft)  / nsamp;  CYf = conj(Yf);
	Pyy = Pyy + Yf .* CYf;

        Yf12(:)  = CYf(mask);
        bic = bic + (Wf * Xf.') .* Yf12;

        ind = ind + nadvance;
    end

    bic     = bic / nrecs;
    Pyy     = Pyy  / nrecs;   Pww = Pww / nrecs;  Pxx = Pxx / nrecs;
    mask(:) = Pyy(mask);
    bic = abs(bic).^2 ./ (Pww * Pxx.' .* mask);
    bic = fftshift(bic) ;

% ------------ contour plot of magnitude bispectum --------------------

   if (rem(nfft,2) == 0)
       waxis = [-nfft/2:(nfft/2-1)]'/nfft;
   else
       waxis = [-(nfft-1)/2:(nfft-1)/2]'/nfft;
   end

   hold off, clf
%   contour(bic,4,waxis,waxis),grid
   contour(waxis,waxis,bic,4), grid on 
   title('Cross-Bicoherence')
   xlabel('f1'), ylabel('f2')
   set(gcf,'Name','Hosa BICOHERX')

   [colmax,row] = max(bic)  ;
   [maxval,col] = max(colmax);
   row = row(col);
   disp(['Max: bic(',num2str(waxis(row)),',',num2str(waxis(col)),') = ', ...
          num2str(maxval) ])

return
