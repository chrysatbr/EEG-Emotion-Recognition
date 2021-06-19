function [h,q] = nlpow(x,y,nfft)
%NLPOW	Second-order Volterra System Identification, arbitrary input
%	[h, q] = nlpow (x, y, nfft);
%	x   - input to the Volterra system
%	y   - output of the Volterra system
%	     x,y must have identical dimensions;
%	     if matrices, columns correspond to realizations.
%	nfft - FFT length to use
%	h   - estimated IR of the linear part
%	q   - estimated IR of the quadratic part
%	    origin at (1,1), with axes pointing right and downwards

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

% Parameter Checks ------------------------------------------------
if (exist('x') ~= 1  | exist('y') ~= 1)
   error('both x and y must be specified')
end

if (size(x) ~= size(y))
   error('x and y must have the same dimensions')
end

[nsamp,nreal] = size(x);
if (nsamp == 1)
   x = x(:);  y = y(:);  nsamp = nreal; nreal = 1;
end

if (exist('nfft') ~= 1) nfft = 0; end
lfft =  2^nextpow2(nsamp);
if (lfft/2 == nsamp) lfft = nsamp; end
nfft = max(nfft,lfft);
x = [x; zeros(nfft-nsamp,nreal)];
y = [y; zeros(nfft-nsamp,nreal)];

% Compute the FT's --------------------------------------
xf = fft(x,nfft);
yf = fft(y,nfft);

hf = zeros(nfft/2+1,1);
qf = zeros(nfft/2,nfft/2);

% Set up equations for each sum frequency ----------------
for m=0:nfft/2
   oddflag = (rem(m,2) == 1);
   if (m == 0)
       ind1 = 0:(nfft/4 -1); ind2 =-ind1;
   elseif (oddflag)
      ind1 = [(m+1)/2 : nfft/4]';
      ind2 = [(m-1)/2 :-1:m-nfft/4]';
   else
      ind1 = m/2 : nfft/4;
      ind2 = m/2 : -1 : m -nfft/4;
   end
   iloc = ind1(1);   indx = ind2;
   ind1 = ind1  + (ind1 < 0) * nfft + 1;
   ind2 = ind2  + (ind2 < 0) * nfft + 1;

   xmat = [xf(m+1,:) ; xf(ind1,:) .* xf(ind2, :) ];
   [l1,l2] = size(xmat);
   Amat = conj (xmat * xmat') / nreal;
   lvec = conj (xmat * yf(m+1,:)' ) / nreal;
   hvec = Amat \ lvec;
   hf(m+1) = hvec(1);


   for i = 1:l1-1
       i1 = nfft/4-1 + i + iloc;            j1 = nfft/4+1 - indx(i);
       i2 = nfft/4-1 + i - iloc + oddflag;  j2 = nfft/4+1 +  iloc + i - 1;
       qf(j1, i1) = hvec(i+1);
       if (i ~= l1-1)
            qf(j2, i2) = hvec(i+1)';
       end
   end
end


% Fill up entire region by symmetry
qf = fliplr(qf);
qdiag = diag(qf);
qf  = fliplr (qf + qf.' - diag(qdiag) );


% For the quadratic part, frequencies are in the range [-0.25,0.25]
% Pack the FD matrix properly for the inverse

N = nfft;
L = nfft/2;
K = nfft/4;

qff = zeros(N,N);
qff (1:K+1, 1:K+1) = qf(K+1:-1:1, K:L);
qff (1:K+1, N-K+2 : N) = qf(K+1:-1:1, 1:K-1);
qff (N-K+2 : N, 1:K+1) = qf(L:-1: K+2, K:L);
qff (N-K+2 : N, N-K+2 : N) = qf (L:-1:K+2, 1:K-1);

q = ifft2(qff);
q = real(q);
q = q(1:K,1:K);      % assumed to be causal

hff = [hf ;  conj(hf(L:-1:2))];
h   = real(ifft(hff));

w1 = [0:nfft/2]/nfft;
w2 = [-nfft/4+1:nfft/4]/nfft;


clf,
subplot(221)
semilogy(w1,abs(hf)), title('linear TF'), xlabel('f'), grid on 
subplot(222)
%contour(flipud(abs(qf)), 6, w2, w2), title('quadratic TF')
contour(w2,w2,flipud(abs(qf)), 6), title('quadratic TF'), grid on 
xlabel('f1'), ylabel('f2')

subplot(223)
plot(h),  title('linear part: IR'), xlabel('t1'), grid on 
subplot(224)
contour(q), title('quadratic part: IR'), grid on 
xlabel('t1'), ylabel('t2')
set(gcf, 'Name','Hosa NLPOW')
