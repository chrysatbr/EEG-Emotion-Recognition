function [hest,ceps] = bicepsf (y,nlag,nsamp, overlap,flag, nfft, wind)
%BICEPSF Non-parametric IR estimation using the bicesptrum (FFT method)
%	[hest, ceps] = bicepsf (y,nlag,segsamp,overlap,flag,nfft, wind)
%	y       - data vector or time-series
%	nlag    - number of lags to compute [must be specified]
%	segsamp - samples per segment    [default: row dimension of y]
%	overlap - percentage overlap     [default = 0]
%	flag    - 'biased' or 'unbiased' [default is 'unbiased']
%	nfft    - FFT length to use      [default = 128]
%	wind    - window function to apply:
%	      if wind=0, the Parzen window is applied (default);
%	      otherwise the hexagonal window with unity values is applied.
%       hest    - estimated impulse response
%       ceps    - estimated complex cepstrum

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



% ------------- parameter checks -------------------------------

    [ly, nrecs] = size (y);
    if (ly == 1) y=y(:);   ly = nrecs; nrecs = 1;      end
    if (exist('overlap') ~= 1)   overlap = 0;          end
    overlap = min(99, max(overlap,0));
    if (nrecs > 1)               overlap = 0;          end
    if (exist('nsamp') ~= 1)     nsamp   = ly;         end
    if (nsamp > ly | nsamp <= 0) nsamp   = ly;         end
    if (exist('flag') ~= 1)      flag    = 'biased';   end
    if (flag(1) ~= 'b')          flag    = 'unbiased'; end
    if (exist('nfft') ~= 1)      nfft    = 128;        end
    if (nfft <= 0)               nfft    = 128;        end
    if (exist('wind') ~= 1)      wind    = 0;          end

    nlag = min(nlag, nsamp-1);
    if (nfft  < 2*nlag+1)   nfft = 2^nextpow2(nsamp); end

% ---------------- create the lag window --------------------
    Bspec = zeros(nfft,nfft) ;

    if (wind == 0)
         indx = (1:nlag)';
         window = [1; sin(pi*indx/nlag) ./ (pi*indx/nlag)];
     else
         window = ones(nlag+1,1);
     end
     window = [window; zeros(nlag,1)];

% ---------------- cumulants in non-redundant region -----------------
% define cum(i,j) = E conj(x(n)) x(n+i) x(n+j)
% for a complex process, we only have cum(i,j) = cum(j,i)
%

     overlap  = fix(nsamp * overlap / 100);
     nadvance = nsamp - overlap;
     nrecord  = fix ( (ly*nrecs - overlap) / nadvance );

     c3 = zeros(nlag+1,nlag+1);
     ind = [1:nsamp]';
     for k=1:nrecord,
         x = y(ind); x = x - mean(x);
         ind = ind + nadvance;
         for j=0:nlag
             z = x(1:nsamp-j) .* x(j+1:nsamp);
             for i=j:nlag
                 sum = z(1:nsamp-i)' * x(i+1:nsamp);
                 if (flag(1) == 'b'), sum = sum/nsamp;
                 else, sum = sum / (nsamp-i);
                 end
                 c3(i+1,j+1) = c3(i+1,j+1) + sum;
             end
         end
     end
     c3 = c3 / nrecord;

% cumulants elsewhere by symmetry  ------------------------------------------
     c3 = c3 + tril(c3,-1)';           % complete I quadrant
     c31 = c3(2:nlag+1,2:nlag+1);
     c32 = zeros(nlag,nlag);  c33 = c32;  c34 = c32;
     for i=1:nlag,
         x = c31(i:nlag,i);
         c32(nlag+1-i,1:nlag+1-i) = x';
         c34(1:nlag+1-i,nlag+1-i) = x;
         if (i < nlag)
           x = flipud(x(2:length(x)));
           c33 = c33 + diag(x,i) + diag(x,-i);
         end
     end

     cmat = [ [c33, c32, zeros(nlag,1)]; [ [c34; zeros(1,nlag)] , c3 ] ];


% ----------- apply lag-domain window  -----------------------------

wcmat = cmat;
mcmat = cmat;
if (wind ~= -1)
     indx = [-nlag:nlag]';
     for k=-nlag:nlag
         wcmat(:,k+nlag+1) = cmat(:,k+nlag+1) ...
            .* window(abs(indx-k)+1) .* window(abs(indx)+1)  ...
             * window(abs(k)+1);
         mcmat(:,k+nlag+1) = k*wcmat(:,k+nlag+1);       % k R(k,n)
     end
end

% ------ compute the cepstrum and the impulse response ---------

    Bspec = fft2(mcmat, nfft, nfft) ./ fft2(wcmat, nfft,nfft);

%    Bceps = fftshift(ifft2(Bspec, nfft, nfft));
%    om = [-nfft/2:nfft/2-1] / nfft;
%    clf, contour(abs(Bceps),8,om,om), grid,
%    title('bicepstrum'), pause

    ceps  = fftshift(diag(real(ifft2(Bspec, nfft,nfft))));
    ceps  = ceps ./ [-nfft/2:-1,1,1:nfft/2-1]' ;
    ceps(nfft/2+1) = 1;
    hest = fftshift(real(ifft(exp(fft(fftshift(ceps))))));
    clf, subplot(211)
    plot(-nfft/2:nfft/2-1, ceps),title('complex cepstrum'), grid on
    subplot(212)
    plot(-nfft/2:nfft/2-1, hest),title('impulse response'), grid on
    set (gcf, 'Name','Hosa BICEPSF')

return
