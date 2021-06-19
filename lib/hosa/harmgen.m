function zmat = harmgen(default);
%HARMGEN Synthetics for the harmonic retrieval problem
%	zmat = harmgen(default)
%	If parameter 'default' is passed, the default values are used;
%	otherwise, the user is prompted for all parameters

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.3 $
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

%------- prompt user for parms ---------------------

if (exist('default') ~= 1) default = 0; else default = 1; end
if (default)
   nsamp = 128;      nrealize = 64;
   nsin   = 2;       fsamp = 1;
   nvar   = 0.5;     mavec = [1];  arvec = [1  -1.058  0.81];
   theta  = [0.1   1  0;    0.2   1  0] ;
   rand('seed',0),   randn('seed',0)
else
  disp(' Synthetics For Harmonic Retrieval Problems ')
  nsamp = 0;
  while (nsamp <= 0)
     nsamp   = input(' samples per realization     ---> [64] ');
     if (isempty(nsamp)) nsamp = 64; end
  end
  nrealize = 0;
  while (nrealize <= 0)
    nrealize = input(' number of realizations      ---> [64] ');
    if (isempty(nrealize)) nrealize = 64; end
  end

  nsin   = input('how many harmonics ?         ---> [2] ');
  if (isempty(nsin)) nsin = 2; end
  fsamp  = input('sampling frequency           ---> [1] ');
  if (isempty(fsamp)) fsamp = 1; end
  nvar   = input('noise variance               ---> [0] ');
  if (isempty(nvar)) nvar = 0; end
  if (nvar > 0)
     disp('  ')
     disp('Remember to enclose vectors within [  ]')
     mavec = input('MA filter for gaussian noise ---> [1] ');
     if (isempty(mavec)) mavec = 1; end
     unstable = 1;
     while (unstable)
     	arvec = input('AR filter for gaussian noise ---> [1] ');
        if (isempty(arvec)) arvec = 1; end
	unstable = any(abs(roots(arvec)) >= 1);
	if (unstable)
	   disp('Unstable AR polynomial: try again');
	end
     end
  end

  disp('  ')
  disp(['frequencies must be normalized to range (0,', ...
       num2str(fsamp/2),')'])
  theta = zeros(nsin,3);
  txtf = [' frequency for harmonic number '];
  for i = 1:nsin
    while (theta(i,1) >= fsamp/2 | theta(i,1) <= 0)
      txt1 = [txtf,int2str(i),' ---> [',num2str(fsamp*i/(5*nsin)),'] '];
      frq = input(txt1);
      if (isempty(frq)), theta(i,1) = fsamp*i/(5*nsin);
         else, theta(i,1) = frq; end
    end
    while (theta(i,2)  <= 0)
      txt2 = [' amplitude for harmonic number ',int2str(i),' ---> [1] '];
      amp = input(txt2);
      if (isempty(amp)), theta(i,2) = 1; else, theta(i,2) = amp; end
    end
  end

end    %  end of parameter selection

  twopi = 2. * pi;
  theta(:,1) = theta(:,1) * twopi / fsamp;
  tvec = (0:nsamp-1)'; alpha = ones(nsamp,1);
  zmat = zeros(nsamp,nrealize);
  y = zeros(nsamp,1); acgn = y;
  nstd = sqrt(nvar);

%
% ---------------- Generate realizations --------------------------
%

  for realize = 1:nrealize
      theta(:,3) = rpiid(nsin,'uni') * twopi;
      y = cos(tvec * (theta(:,1))' + alpha * (theta(:,3))') * theta(:,2);
      if (nvar > 0)
	 acgn = rpiid(nsamp,'nor');
         acgn = acgn - mean(acgn);
         acgn = filter(mavec,arvec,acgn);
         acgn = (acgn - mean(acgn)) * nstd / std(acgn);
         y = y + acgn;
      end
      zmat(:, realize) = y;
  end

return
