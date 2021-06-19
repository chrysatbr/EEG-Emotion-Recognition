function zmat = armasyn (default)
%ARMASYN Generates ARMA synthetics.
%	zmat = armasyn(default)
%	if default > 0, default settings are used, otherwise
%	the user is prompted for all parameters.
%	zmat is the  generated ARMA synthetics; each column corresponds
%               to a different realization

%   ar, ma:     arma parameters for the input noise
%   in_type:    input noise pdf (exp, lap, uni, nor or bga)
%   nvar:       noise variance
%   ar_n, ma_n: arma parameters for the additive noise
%   n_type:     additive noise pdf (exp, lap, uni, or nor)


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

% --------------------------------------------------------------

if (exist('default') ~= 1) default = 0; else default = 1; end
if (default)
  rand('seed',0), randn('seed',0)
  N = 1024;    nrun = 1;
  ar = [1, -1.5, 0.8];  ma = [1]; in_type = 'exp';
  nvar = 0.01; ar_n = [1];  ma_n = [1];    n_type = 'nor';
  p_spike = 0;
else
  N       = input(' samples per realization   [256] --->');
  if (isempty(N)) N = 256; end
  if (N <= 0) N = 256; end
  nrun    = input(' number of realizations      [1] --->');
  if (isempty(nrun)) nrun=1; end
  if (nrun <= 0) nrun = 1; end

% ensure that user-specified AR polynomial is stable
  stab_flag = 1;
  while (stab_flag ==1)
    ar      = input(' AR coefficients: signal     [1] --->');
    if (isempty(ar)) ar=1; end
    stab_flag = any(abs(roots(ar)) >= 1);
    if (stab_flag)
       disp('Unstable AR polynomial: try again ')
    end
  end
  ma        = input(' MA coefficients: signal     [1] --->');
  if (isempty(ma)) ma=1; end

  in_type = 'xyz';
  while (in_type ~= 'uni' & in_type ~= 'exp' & in_type ~= 'lap' ...
       & in_type ~= 'nor' & in_type ~= 'bga')
    in_type = ...
    input(' signal pdf: ''uni'' ''exp'' ''lap'' ''nor'' ''bga'' [''exp''] --->');
    if (isempty(in_type)) in_type ='exp'; end
  end

  p_spike = 0.;
  if (in_type == 'bga')
     while (p_spike <= 0. | p_spike >= 1)
        p_spike = input(' input: probability of spike for bga  --->');
        if (isempty(p_spike)) p_spike=0; end
     end
  end
  disp(' ')
  disp('The noise free signal will be normalized to unity variance')
  nvar  = input(' noise variance              [0] --->');
  if (isempty(nvar)) nvar =0; end

  ar_n = [1]; ma_n = [1];
  if (nvar > 0)
    stab_flag = 1;
%           ensure that user-specified AR polynomial is stable
    while (stab_flag ==1)
      ar_n    = input(' AR coefficients: noise      [1] --->');
      if (isempty(ar_n)) ar_n = 1; end
      stab_flag = any(abs(roots(ar_n)) >= 1);
      if (stab_flag)
         disp('Unstable AR polynomial: try again ')
      end
    end
    ma_n    = input(' MA coefficients: noise      [1] --->');
    if (isempty(ma_n)) ma_n = 1; end

    n_type = 'xyz';
    while (n_type ~= 'uni' & n_type ~= 'exp' & n_type ~= 'lap' ...
         & n_type ~= 'nor')
         n_type  = input(' noise pdf: ''uni'' ''exp'' ''lap'' ''nor'' [''nor''] --->');
         if (isempty(n_type)) n_type = 'nor'; end
    end
  end

end
zmat = zeros(N,nrun);

for i=1:nrun                   % generate realizations

    u = rpiid(N,in_type,p_spike);
    y = filter(ma,ar,u);
    y = y / std(y);

    if (nvar > 0)
       w   = rpiid(N,n_type);
       y_n = filter(ma_n,ar_n,w);
       y_n =  y_n * sqrt(nvar) / std(y_n);
       y   = y + y_n;
    end
    zmat(:,i) = y;
end

return
