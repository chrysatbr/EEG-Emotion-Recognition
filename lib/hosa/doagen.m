function smat = doagen(default)
%DOAGEN	Synthetics for the DOA problem.
%	smat = doagen(default)
%	If parameter 'default' is passed, the default settings are used;
%	otherwise, the user is prompted for all parameters.
%
%	smat is the nsamp x msens matrix of noisy sensor signals
%	     nsamp is the number of samples,
%	     msens is the number of sensors.
%

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

% ----------------------------------------------------------------
  format compact

if (exist('default') ~= 1) default = 0;  else default = 1; end
def_msens = 8;
def_dspace = 0.5;
def_nsamp = 4096;
def_msource = 2;
% bearing  =  [-15 -25];   2/17 commented out: problems if #sources=1
in_type = 'lap';
amp = [1, 1];
nvar = 1;
ar = [1 0 0.81];
ma = 1;
if (default)
   bearing = [-15, -25];
   rand('seed',0); randn('seed',0);
   for j=1:def_msource
       smat(:,j) = rpiid(def_nsamp,in_type)*amp(j);
   end
   msens = def_msens; dspace = def_dspace; nsamp = def_nsamp;
end

if (default <= 0)
% --- Number of sensors
  txt = ['Number of sensors                      ---> [', ...
        int2str(def_msens),'] '];
  msens = input(txt);
  if (isempty(msens)), msens = def_msens; end
  if (msens <= 0)
     msens = def_msens;
     disp(['   Number of sensors set to ',int2str(def_msens)])
  end

% --- Sensor spacing
  txt = ['Sensor spacing in units of wavelengths ---> [', ...
         num2str(def_dspace),']' ];
  dspace = input(txt);
  if (isempty(dspace)), dspace = def_dspace; end
  if (dspace <= 0.)
      dspace = def_dspace;
      disp(['  Sensor spacing set to ', num2str(def_dspace)])
  end

% --- Number of snapshots
  txt = ['Sensor signal length in samples        ---> [', ...
        int2str(def_nsamp),'] '];
  nsamp = input(txt);
  if (isempty(nsamp)), nsamp = def_nsamp; end
  if (nsamp <= 0.)
     nsamp = def_nsamp;
     disp(['  Sensor signal length set to ',int2str(def_nsamp)])
  end

% --- Number of sources
  txt = ['Number of sources                      ---> [', ...
        int2str(def_msource),'] '];
  msource = input(txt);
  if (isempty(msource)) msource = def_msource; end
  if (msource <= 0)     msource = def_msource; end


% ---- Source parameters: Bearing, Amplitude, P.D.F.
 disp(' ')
 disp(' Note: standard beamwidth is  pi/msens');
 disp(['Allowed source signal pdfs are: ''exp'', ''lap'', ''uni'', ''nor'''])

  smat = zeros(nsamp,msource);
  for n = 1:msource
      disp(' ')
      disp(['---- Specify parameters for source number ', num2str(n),' :'])
      bear0 = -5 -10*n;
      bear = input(['Source bearing in degrees  ---> [',num2str(bear0),']']);
      if (isempty(bear)) bear = bear0; end
      bearing(n) = bear;

      amp     = input('source amplitude           ---> [1] ');
      if (isempty(amp)) amp = 1; end
      if (amp <= 0) amp = 1; end

      in_type = input('source signal pdf          ---> [''lap''] ');
      if (isempty(in_type)), in_type ='lap'; end
      if (~strcmp(in_type,'exp') & ~strcmp(in_type,'lap') & ...
          ~strcmp(in_type,'uni') & ~strcmp(in_type,'nor'))
         in_type = 'lap';
	 disp(['  Source signal pdf set to ',in_type])
      end
      smat(:,n) = rpiid(nsamp,in_type) * amp;
  end

% noise:
  nvar = input('noise variance                        ---> [0] ');
  if (isempty(nvar)) nvar = 0; end

  if (nvar > 0.)

% --- Noise color (spatial ARMA)
     disp(' ')
     disp(['the AR(2) vector [1,-2*r*a,r^2] will have peaks at theta'])
     disp(['    where a = cos(sin(theta)*pi) '])
     unstable = 1;
     while (unstable)
         ar = input('AR vector for the noise          ---> [1]');
         if (isempty(ar)) ar = 1; end
         unstable = any(abs(roots(ar)) >= 1);
	 if (unstable)
	    disp('Unstable AR polynomial: try again ')
	 end
     end
     ma = input('MA vector for the noise          ---> [1]');
     if (isempty(ma)) ma = 1; end
  end
end                            %----> end of parameter prompting

% the propagation matrix
  bearing = bearing * pi / 180;    % degrees to radians
  amat = exp(sqrt(-1) * [0:msens-1]' * (2*pi*dspace*sin(bearing)));

% the noisefree sensor signal
  smat = smat * amat.';
% --- Convert ARMA noise to equivalent MA form:
  if (nvar > 0)
     cvec = cumtrue(ma,ar,2,msens);
     q = (length(cvec) - 1) /2 ;
     if (q < msens)
        cvec = [zeros(msens-q,1); cvec; zeros(msens-q,1)];
     end
     cmat = toeplitz(cvec(msens+1:2*msens));
     amat = chol (cmat);

% --- Generate white noise, then color it:
     gmat = reshape(rpiid(nsamp*msens,'nor'),nsamp,msens);
     gmat = gmat * amat;
     lby21 = fix((nsamp-1)/2);
     mask  = [1; ones(lby21,1)*2;1; zeros(lby21,1)];
     for i=1:msens
        gmat(:,i) = ifft( fft(gmat(:,i)).*mask );        % complex noise
	g = gmat(:,i) - mean(gmat(:,i));
	q = sqrt( real(g'*g) / (nsamp-1) );
	gmat(:,i) = gmat(:,i) / q * sqrt(nvar);
     end
     smat = smat + gmat;
  end

return
