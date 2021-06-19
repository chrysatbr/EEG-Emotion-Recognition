function [hest,ceps,A,B,minh,maxh] = biceps(y,p,q,nsamp,overlap,flag,lh) 
%BICEPS	Non-parametric system identification using the bicepstrum. 
%	BICEPS estimates the complex cepstrum, using third-order
%	cumulants, and then reconstructs the impulse response. 
%
%	[hest,ceps,A,B,minh,maxh] = biceps(y,p,q,sampseg,overlap,flag,lh)  
%	      y - data vector or matrix
%	    p,q - cepstral orders (truncation points)   [MUST be specified]
%	sampseg - samples per record                    [default = row size] 
%	overlap - percentage overlap of records         [default = 0; max=99] 
%	  flag  - 'biased' or 'unbiased'                [default = 'biased'] 
%	    lh  - hest indices will range from -lh to lh [default = 2(p+q)]
% 
%	if y is a matrix, each column is assumed to be a different realization
%	or record; in this case, overlap is set to 0, and record size is set 
%	to the row length of the matrix. 
% 
%	hest - estimated impulse response, h(n):  n= -lh, ... ,lh  
%	ceps - estimated complex cepstrum 
%	A    - estimated A(k)'s, k=1,...,p
%	B    - estimated B(k)'s, k=1,...,q 
%       minh - minimum phase component of hest; length lh+1
%	maxh - maximum phase component of hest; length lh+1 

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.7 $
%  A. Swami   January 20, 1993.   Revised Jan 20, 1995. 

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

%--------------- parameter checks ----------------------- 

  [ly, nrecord]  = size(y); 
  if (ly == 1) y = y(:); ly = nrecord; nrecord = 1;       end 

  if (exist('p') ~=1) error('order p must be specified'); end 
  if (exist('q') ~=1) error('order q must be specified'); end 

  if (nrecord > 1)   nsamp = ly;  overlap = 0;            end 
  if (exist('nsamp') ~= 1)          nsamp = 0;            end 
  if (nrecord > 1)                  nsamp = ly;           end 
  if (exist('overlap') ~= 1)      overlap = 0;            end
  if (nrecord > 1)                overlap = 0;            end 
  if (exist('flag') ~= 1)            flag = 'biased';     end
  if (exist('lh')   ~= 1)              lh = 2*(p+q);      end 

  overlap = min(99,(max(overlap,0))); 
  nsamp   = min(ly,max(nsamp,0)); 

% ------ estimate cumulants ---------------------------------------
  
  w     = max(p,q);  
  w     = w + rem(w,2);  
  wby2  = w/2; 
  cum_y = zeros(w+p+q+1,4*w+1); 
  M     = 2 * w; 
  n_loc = wby2+w+1; 

  for k = -wby2-w:wby2+w
     cum_y(n_loc+k,:) = cumest(y(:),3,M,nsamp,overlap,flag,k).';     
  end 

% -------- set up system of linear equations ----------------------
%    \sum_{i=-q}^{p} i ch(i) [ C(m-i,n) - C(m+1,n+i) ] = m C(m,n) 
%           with, |m| <= max(p,q) and |n| < max(p,q)/2


  xorg = M + 1; 
  yorg = wby2 + w + 1;        % c3(0,0) is in cum_y(xorg,yorg) 
  ind1 = (-w+q:w+q) + xorg; 
  ind2 = (-w+q:-1:-w-p) + xorg; 
  ind3 = yorg + (q:-1:-p);
  ind4 = xorg + (-w:w); 

  Amat = []; 
    for n = -wby2:wby2 
        alpha = cum_y(yorg+n,:);                     % picks up (.,n) slice 
        tmp   = toeplitz(alpha(ind1), alpha(ind2));  % I term b(m-k,n) m=-w:w
        tmp1  = cum_y(ind3-n,ind4-n).';                        
        Amat  = [Amat; tmp-tmp1]; 
    end 

% ------ solve and convert to complex cepstral coefficients -----------

    Amat      = Amat(:,[1:q,q+2:p+q+1]); 
    rvec      = cum_y(yorg+ (-wby2:wby2), ind4) * diag(-w:w); 
    rvec      = rvec'; 
    rvec      = rvec(:); 
    dcepstrum = Amat\rvec; 
    cepstrum  = dcepstrum ./ [-q:-1,1:p]'; 
    ceps      = [cepstrum(1:q); 0; cepstrum(q+1:q+p)]; 


% ------- estimate causal/anti-causal IR's directly: 

    maxh   = real( ifft( exp(  fft( [1;cepstrum(q:-1:1) ], lh+1 )  )));
    minh   = real( ifft( exp(  fft( [1;cepstrum(q+1:q+p)], lh+1 )  )));
    hest   = conv(flipud(maxh),minh); 
    A      = - dcepstrum(q+1:q+p);
    B      =   dcepstrum(q:-1:1); 

% ------- or estimate impulse response using Oppenheim and Schafer's method 
opschaf = 0; 
if (opschaf == 1) 
    cpos   = cepstrum(q+1:q+p) .* [1:p]'; 
    cneg   = cepstrum(q:-1:1)  .* [1:q]';  

    ir_len = lh; 
    cpos   = [cpos; zeros(ir_len-p,1)];
    cneg   = [cneg; zeros(ir_len-q,1)];
   
    ha     = zeros(ir_len+1,1); hc = ha;  ha(1) = 1; hc(1) = 1; 
    for n=1:ir_len 
        ha(n+1) = cneg(1:n).' * ha(n:-1:1) /n; 
        hc(n+1) = cpos(1:n).' * hc(n:-1:1) /n; 
    end
    ha   = ha(ir_len+1:-1:1);
    hest = conv(ha,hc); 
    minh = hc; maxh = ha; 
    hest = hest/hest(2*(p+q)+1); 
    A    = -cpos(1:p);
    B    = -cneg(1:q); 
end

% --------- display estimates --------------------


    clf, subplot(211)
    plot([-q:p]', ceps),  grid on
    title('complex cepstrum')
    xlabel('sample number')

    subplot(212) 
    taxis = [-lh:lh]'; 
    plot(taxis, hest), grid on 
    title('impulse response ') 
    xlabel('sample number') 
    set (gcf, 'Name','Hosa BICEPS')
return

