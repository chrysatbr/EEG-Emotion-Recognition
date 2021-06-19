function [hest] = matul (bisp)
%MATUL	Impulse response estimation using Matsuoka-Ulrych algorithm
%	hest = matul(bisp)
%	bisp  - the estimated bispectrum
%	        (e.g., as computed by bispecd or bispeci).
%	 hest - estimated impulse response

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.7 $
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

% ----- parameter checks ----------------
  [nfft,mfft]=size(bisp);
  if (rem(nfft,2) == 0),
      N = nfft/2 + 1;
      bisp =[bisp(nfft,N:nfft), bisp(nfft,1);
            bisp(1:N-1,N:nfft), bisp(1:N-1,1)];
   else,
      N = (nfft+1)/2;
      bisp = bisp(1:N,N:nfft);
   end

% -------------------------------
% terms at origin and nyquist
  B00  = bisp(1,1);

  psi1   = angle(bisp(1,:));
  mshift = (psi1(1)-psi1(2)) * nfft/(2*pi);
  x      = exp(sqrt(-1)*mshift*2*pi/nfft *[0:N-1]);
  bisp   = bisp .* (x.' * x);

% Step 0: extract the psi(1,:) slice  and obtain raw estimate
  psi1 = angle(bisp(2,2:N-1))';     % psi(1,k) k=1,...,M/2 -1
  phi1 = [0; -cumsum(psi1)];        % phi(1)=0; phi(2),...,phi(M/2).

  phi2 = phi1(1:N-2) - (1:N-2)'/(nfft/2) * phi1(N-1);
  phi2(N-1) = phi1(N-1);
  phi2 = rem(phi2,pi);

% set up the matrix for the matsuoka-ulrych algorithm

K=fix(N/2); L=K*(N-K);
amat = zeros(L,N); rpsi = zeros(L,1); rmag = zeros(L,1);

loc = 0;
for i=1:fix(N/2)
    j = N+1-2*i;
    ind = loc+1:loc+j;
    amat(ind,i)     =  ones(j,1);
    amat(ind,i:N-i) = amat(ind,i:N-i) + eye(j);
    amat(ind,i+(i:N-i)) = amat(ind,i+(i:N-i)) - eye(j);
    rpsi(ind) = angle(bisp(i+1,i+1:N-i+1))';
    rmag(ind) = log(abs(bisp(i+1,i+1:N-i+1)))';
    loc = loc + j;
end

% --- find the 2pi corrections due to Rangoussi-Giannakis ---

kwrap = fix( (amat(:,1:N-1)*phi2 - rpsi(:)) /(2*pi));
phi = amat(:,1:N-1) \ (rpsi(:) + 2*pi*kwrap);
mag = abs(amat(:,1:N-1)) \ rmag(:);

% --- combine to get tf; and hence, ir.

phi0 = pi * (B00 < 0);
mag0 = abs(B00)^(1/3);
phi(N-1)=0;
phz = [phi0; phi; -flipud(phi(1:N-2))];
mag = exp(mag);
mag = [mag0; mag; flipud(mag(1:N-2))];

hf = mag .* exp(sqrt(-1)*phz);
hest  = real(fftshift(ifft(hf)));

hest = hest/max(abs(hest));
clf, plot(-(N-2):(N-1),hest), grid on 
title('Estimated impulse response')
xlabel('sample number')
set(gcf, 'Name','Hosa MATUL')

return
