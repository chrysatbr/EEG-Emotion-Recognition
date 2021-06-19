function  [loc,val] = pickpeak(spec,npicks,rdiff)
%PICKPEAK Picks peaks 
% [loc,val] = pickpeak(spec,npicks,rdiff)
%	spec   - data vector or matrix 
%	npicks - number of peaks desired              [default = 2]
%	rdiff  - minimum spacing between picked peaks [default = 5]
%       loc    - vector of locations (indices) of the picked peaks
%	val    - vector corresponding values 
%	A 0 in location (i,j) of array loc (or a NaN in array val)
%	indicates that the j-th data vector has less than i peaks
%	with a separation of rdiff or more. 

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.4 $
%  A. Swami Jan 20, 1995. 

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


% ---- parameter checks  -------------------------------------------

if (exist('rdiff') ~= 1)  rdiff =  5;                  end 
if (exist('npicks') ~= 1) npicks = 2;                  end 

% ---- convert row vectors to col vectors  -------------------------

[mrows,ncols]  = size(spec);
if (mrows==1) mrows=ncols; ncols=1; spec = spec(:);   end

% ---- edit out NaNs and Infs ---------------------------------------

good = find (finite(spec)); 
rmin = min(spec(good)) - 1; 
bad  = find(~finite(spec));
if (~isempty(bad)) 
   spec(bad) = ones(size(bad)) * rmin; 
end 

% ---- find a peak, zero out the data around the peak, and repeat 

val =  ones(npicks,ncols) * NaN ; 
loc =  zeros(npicks,ncols) ; 

for k=1:ncols
                                           % Find all local peaks: 
    dx = diff([rmin; spec(:,k); rmin]);    % for a local peak at either end 
    lp = find(dx(1:mrows)   >= 0 ... 
            & dx(2:mrows+1) <=0);          % peak locations 
    vp = spec(lp,k);                       % peak values

    for p=1:npicks
       [v,l] = max(vp);                   % find current maximum
       val(p,k) = v;  loc(p,k) = lp(l);   % save value and location 

       ind = find(abs(lp(l)-lp) > rdiff);  % find peaks which are far away

       if (isempty(ind)) 
           break                           % no more local peaks to pick
       end 
       vp  = vp(ind); 			   % shrink peak value array
       lp  = lp(ind);                      % shrink peak location array 
    end
end 
