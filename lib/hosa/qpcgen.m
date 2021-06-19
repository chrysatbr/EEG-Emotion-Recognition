function zdat = qpcgen(default)
%QPCGEN	Generates synthetics for the quadratic phase-coupling problem.
%	zdat = qpcgen(default)
%	If parameter 'default' is passed, the default values are used;
%	otherwise, the user is prompted for all parameters

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

% --------- prompt for frequencies and amplitudes of coupled harmonics

if (exist('default') ~= 1) default = 0; else default = 1; end
if (default)
   nsamp = 64;  nrun = 64; nvar = 1.5;  mavec = 1; arvec = 1;
   fsamp = 1;   nqpc = 1;  frq = [0.1 0.15 0.25];  amq = [1 1 1];
   uqpc  = 1;   fru  = 0.4;  amu = 1;
   rand('seed',0),  randn('seed',0)
else

   disp(' Synthetics For Quadratic Phase Coupling ')
   nsamp = 0;
   while (nsamp <= 0)
     nsamp    = input('samples per realization       ---> [64] ');
     if (isempty(nsamp)) nsamp=64; end
   end
   nrun = 0;
   while (nrun <= 0)
     nrun    = input('number of realizations        ---> [64] ');
     if (isempty(nrun)) nrun=64; end
   end

   nvar   = input('noise variance               ---> [0] ');
   if (isempty(nvar)) nvar = 0; end

   if (nvar > 0)
      disp(' ')
      disp('Remember to enclose vectors within [  ]')
      mavec = input('MA filter for gaussian noise ---> [1] ');
      if (isempty(mavec)) mavec = 1; end
      unstable = 1;
      while (unstable)
         arvec = input('AR filter for gaussian noise ---> [1] ');
         if (isempty(arvec)) arvec = 1; end
  	 unstable = any(abs(roots(arvec)) >= 1.);
	 if (unstable)
	    disp('Unstable AR polynomial: try again')
	 end
      end
   end

   fsamp  = input('sampling frequency           ---> [1] ');
   if (isempty(fsamp)) fsamp = 1; end

   nqpc = input('number of phase-coupled frequency triplets ---> [1] ');
   if (isempty(nqpc)), nqpc = 1; end
   if (nqpc > 0)
      disp(' ')
      disp(['frequencies must be in the range (0,', num2str(fsamp/2),')'])
      disp('recall that f3 = f1 + f2 in the qpc problem')
      frq = zeros(nqpc,3);  amq = ones(nqpc,3);
      frq(:,1) = fsamp/(10*nqpc) * [1:nqpc]';
      frq(:,2) = 2 * fsamp/(10*nqpc) * [1:nqpc]';
   end

   txtf1 = 'frequency of first harmonic  ---> [';
   txtf2 = 'frequency of second harmonic ---> [';

   for k=1:nqpc
      disp(' ')
      disp(['Phase-coupled triplet number ',int2str(k)])
      while (frq(k,3) <= 0. | frq(k,3) >= fsamp/2 )
         txtf = [txtf1, num2str(frq(k,1)),'] '];
         fr1 = 0.;
         while (fr1 <= 0. | fr1 >= fsamp/2)
            fr1 = input(txtf);
            if (isempty(fr1)), fr1 = frq(k,1); end
         end
         frq(k,1) = fr1;

         txtf = [txtf2, num2str(frq(k,2)),'] '];
         fr2 = 0.;
         while (fr2 <= 0. | fr2 >= fsamp/2)
            fr2 = input(txtf);
            if (isempty(fr2)), fr2 = frq(k,2); end
         end
         frq(k,2) = fr2;

         frq(k,3) = frq(k,1) + frq(k,2);
      end

      amp  = input(' amplitude of  first harmonic ---> [1] ');
      if (~isempty(amp)),  amq(k,1)  = amp; end
      amp  = input(' amplitude of second harmonic ---> [1] ');
      if (~isempty(amp)),  amq(k,2)  = amp; end
      amp  = input(' amplitude of  third harmonic ---> [1] ');
      if (~isempty(amp)),  amq(k,3)  = amp; end
   end

% prompt for frequencies and amplitudes of uncoupled harmonics
   uqpc = input('number of uncoupled harmonics ---> [0] ');
   if (isempty(uqpc)), uqpc = 0; end

   if (nqpc + uqpc == 0)
      error('total number of harmonics = 0!!')
   end

   if (uqpc > 0)
      disp(' ')
      disp(['frequencies must be in the range (0,', num2str(fsamp/2),')'])
      fru = zeros(uqpc,1);  amu = ones(uqpc,1);
      fru = fsamp/(10*uqpc) * [1:uqpc]' + fsamp/20;
   end

   for k=1:uqpc
       disp(' ')
       disp(['Non Phase-coupled harmonic number ',int2str(k)])
       ufrq = 0;
       while (ufrq <= 0. | ufrq >= fsamp/2)
          txtf = [' Frequency for harmonic ',int2str(k), ...
                 ' ----> [',num2str(fru(k,1)),'] '];
          ufrq = input(txtf);
          if (~isempty(ufrq)),
             if (ufrq > 0. & ufrq < fsamp/2) fru(k,1) = ufrq; end
          end
       end
       amp  = input(' amplitude of  the harmonic   ---> [1] ');
       if (~isempty(amp)),  amu(k,1)  = amp; end
   end
end

% other parameters

  freqs = [frq(:); fru(:)]' * 2 * pi /fsamp;
  amps  = [amq(:); amu(:)];

% ------------ generate synthetics --------------------------

  tvec = (0:nsamp-1)'; alpha = ones(nsamp,1);
  zdat = zeros(nsamp,nrun);
  nstd = sqrt(nvar);

  for krun=1:nrun
      thq = rpiid(nqpc*3, 'uni') * 2 * pi;
      thq = reshape(thq, nqpc, 3);
      thu = rpiid(uqpc,'uni') * 2 * pi;

      if (nqpc > 0) thq(:,3) = thq(:,1) + thq(:,2); end
      theta = [thq(:); thu(:)]';
      y = cos(tvec * freqs + alpha * theta) * amps;

      if (nvar > 0)
         acgn = rpiid(nsamp,'nor');
         acgn = acgn - mean(acgn);
         acgn = filter(mavec,arvec,acgn);
         acgn = (acgn - mean(acgn))/std(acgn) * nstd;
         y = y + acgn;
      end
      zdat(:,krun) = y;
  end

return
