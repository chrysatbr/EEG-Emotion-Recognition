function  [s1, s2] = tdegen (default)
%TDEGEN	Synthetics for time-delay estimation
%	[s1, s2] = tdegen(default)
%	 s1, s2  noisy signals at the two sensors
%	If the parameter 'default' is passed, the default settings are
%	used; otherwise, the user is prompted for all parameters.

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

% ------------------------------------------------------------------------

if (exist('default') ~= 1) default = 0; else default = 1; end
if (default)
   nsamp = 4096;    delay = 16;    stype = 'exp';
   nvar1 = 1;       nvar2 = 1;     ntype = 'nor';  amp2  = 1;
   ar1   = 1;       ma1   = 1;     ar2   = 1;      ma2   = [1:6,5:-1:1];
   rand('seed',0);  randn('seed',0);
else
   disp(' Synthetics for Time-Delay Estimation ')
   nsamp = input(' Number of samples  ----> [1024] ');
   if (isempty(nsamp)), nsamp = 1024; end

   disp(' ')
   disp('Note: fractional delays will be rounded off ')
   delay = input(' Signal delay       ----> [16] ');
   if (isempty(delay)) delay = 16; end
   delay = round(delay);

   amp2 = input('Amplitude gain for signal at sensor 2 ---> [1] ');
   if (isempty(amp2)) amp2 = 1; end

   stype = 'exp';

   nvar1 = input(' Noise variance at first sensor ----> [1] ');
   if (isempty(nvar1)), nvar1 = 1; end
   if (nvar1 > 0)
   	nvar2 = input(' Noise variance at second sensor ----> [1] ');
	if (isempty(nvar2)), nvar2 = 1; end
	ntype = input([' Noise type: ''uni'', ''lap'', ''nor'' ', ...
	              '----> [''nor''] ']);
        if (isempty(ntype)), ntype = 'nor'; end
   else
      nvar2 = 0;
   end

   if (nvar1 > 0)
      disp(' ')
      disp(['Specify ARMA parameters for the colored noise ', ...
            'at the first sensor:'])
      disp(' All vectors must be enclosed within [   ]')

      stab_flag = 1;
      while (stab_flag ==1)
      	ar1 = input('AR parameters for noise at sensor 1 ---> [1] ');
	if (isempty(ar1)) ar1=1; end
	stab_flag = any(abs(roots(ar1)) >= 1);
	if (stab_flag)
       	   disp('Unstable AR polynomial: try again ')
        end
      end
      ma1 = input('MA parameters for noise at sensor 1 ---> [1] ');
      if (isempty(ma1)), ma1 = 1; end
   end

   if (nvar2 > 0)
      disp(' ')
      disp('The noise at sensor 1 is passed through an ARMA filter ')
      disp(' to obtain the noise at sensor 2: ')
      disp(' All vectors must be enclosed within [   ]')

      stab_flag = 1;
      while (stab_flag ==1)
      	ar2 = input('AR parameters for noise transfer function ---> [1] ');
	if (isempty(ar2)) ar2=1; end
	stab_flag = any(abs(roots(ar2)) >= 1);
	if (stab_flag)
           disp('Unstable AR polynomial: try again')
	end
      end
      ma2 = input('MA parameters for noise transfer function ---> [1] ');
      if (isempty(ma2)), ma2 = 1; end
   end
end
% ----------

 g = rpiid(nsamp+abs(delay),'exp');

 if (delay >= 0)  s2 = g(1:nsamp);  s1 = g(delay+1:delay+nsamp);
 else,            s1 = g(1:nsamp);  s2 = g(-delay+1:-delay+nsamp);
 end
 s1 = s1/std(s1); s2 = s2/std(s2) * amp2;

 if (nvar1 > 0)
    g = rpiid(nsamp,ntype);
    g = filter(ma1,ar1,g);
    g = (g - mean(g))/std(g) * sqrt(nvar1);
    s1 = s1 + g;
    if (nvar2 > 0)
        g = filter(ma2,ar2,g);
        g = (g - mean(g))/std(g) * sqrt(nvar2);
        s2 = s2 + g;
    end
 end

return
