function  [a,theta,alpha,fr] = hprony (x,p)
%HPRONY	Prony's method for modeling transients
%	[a,theta,alpha,fr] = hprony (x,p)
%	x  - transient signal
%	p  - order    [default = length(x)/10]
%	HPRONY fits the model:
%	     x(n) = sum_{1 to p} A(k) z(k)^n
%	where A(k) = a(k) exp (j * theta(k))
%	and   z(k) = exp( alpha(k) + j * 2 * pi * fr (k) )
%	The estimated values of a(k), theta(k), alpha(k) and fr(k)
%	   k=1,...,p, are returned in the vectors, a, theta, alpha and fr.


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

%------------- parameter checks --------
[nsamp,ncol] = size(x);
if (min(nsamp,ncol) ~= 1)
   error('hprony: x must be a vector');
end

if (nsamp == 1)
   nsamp = ncol; ncol=1; x = x(:);
end

ind = find(abs(x)>0);           %  leading/trailing zeros cause problems
x   = x(min(ind):max(ind));	%  with prony
nsamp = length(x);
xone  = x(1);
x     = x/x(1);                 % to handle v3.5/v4.0 toeplitz incompatibility

if (exist('p') ~= 1) p = fix(nsamp/10); end

% --------------------------
[ma,ar] = prony (x, p-1, p);            % fit an ARMA model

% junk = toeplitz([1; j]);      % to handle v3.5/v4.0 toeplitz incompatibility
			      % prony calls toeplitz
% if (junk(2,1) ~= j )
%    ma = conj(ma);  ar = conj(ar);
% end

ma = ma * xone;

[res,poles,delay] = residue(ma,ar);   % convert to poles-residues

a = abs(res);                         % required conversions
theta = angle(res);
alpha = real(log(poles));
fr    = imag(log(poles)) / (2*pi);


% ----------------------- plot the data and the model ---------------
xest = ( (a .* exp(j*theta)) * ones(1,nsamp) ) .*  ...
	exp( (j * 2*pi* fr + alpha) * [0:(nsamp-1)] );
if (p > 1) xest = sum(xest).';
else       xest = xest(:);
end
x    = x * xone;

if (exist('debug') ==  1)
  if (all(imag(x) == 0))
    plot([x,xest]), 
    title('True (-) and estimated (--)  data')
  else
    clf, subplot(211)
    plot(real([x,xest])), 
    title('True (-) and estimated (--)  data: real part')
    subplot(212)
    plot(imag([x,xest])), 
    title('True (-) and estimated (--)  data: imaginary part')
  end
  mse = (x-xest)'  * (x-xest) / length(x);
  disp(['Mean squared error = ',num2str(mse)])
end


% ---------------------- display parameter estimates ------------

if (all(imag(x)==0))
   disp('AR parameters are:')
   disp(ar')
   disp('MA parameters are:')
   disp(ma')
   disp(['delay = ',int2str(delay)])
else
   format compact
   disp('amplitudes, thetas, alphas, frequencies')
   disp(a.'), disp(theta.'), disp(alpha.'), disp(fr.')
end

return
