function [spec,theta,bearing] = doa(ymat, dspace, dtheta,nsource,order,delta)
%DOA	Direction-of-arrival estimation. 
%	[spec,theta] = doa(ymat, dspace, dtheta,nsource,order,delta)
%	 ymat  - sensor array data: each column corresponds to a different 
%	         sensor; rows correspond to time samples (snapshots)
%	dspace - array spacing "d" in units of wavelengths    default = 0.5
%	dtheta - angular resolution in degrees                default = 2. 
%       nsource - number of sources: default value is 0; the singular values
%                of the cross-cumulant or cross-correlation matrix 
%		 (see `order' below) will be displayed, and the user will
%	         be prompted to enter the number of sources. 
%	order   - cumulant order to use; should be 2 or 4;    default = 4
%	delta  - displacement (in number of elements) for ESPRIT [default = 1]
%	  spec - is the array of estimated ``spectra'';  the columns
%	         correspond to estimates based on the Eigenvector, Music, 
%	           Pisarenko, ML, AR, min-norm and beamformer methods. 
%	 theta -  Vector of  bearings corresponding to the rows of spec
%	bearing - Vector of source bearings (in degrees) estimated by ESPRIT. 

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.6 $
%  A. Swami   January 20, 1993; May 1995 

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

%----------- Parameter checks ------------------------

[nsamp, msens] = size(ymat); 
if (msens == 1 | nsamp == 1)
   error('ymat - should be a matrix')
end

if (exist('dspace')  ~= 1) dspace  = 1/2; end
if (exist('dtheta')  ~= 1) dtheta  = 2;   end 
if (exist('delta')  ~= 1) delta  = 1; end 
if (exist('nsource') ~= 1) nsource = 0; end 
if (exist('order')   ~= 1) order = 4; end
if (order ~= 2 & order ~= 4) 
   error(' cumulant order should be 2 or 4')
end 

%----- estimate the spatial cross-cumulant matrix ----
cmat = zeros(msens,msens); 

ymat  = ymat - ones(nsamp,1) * mean(ymat);   % remove mean 

if (order == 4)
   zmat = ymat .* ymat .* conj(ymat);    % the iv for the fourth moments
   cmat = conj(zmat' * ymat) /nsamp;     % fourth-order cross-moments 
   r1   = conj(ymat' * ymat) / nsamp;    % correlation matrix 
   r2   = conj(ymat.'* ymat) / nsamp;    % moment matrix 

                                         % cross-cumulant matrix: 
   cmat = cmat - 2 * diag(diag(r1)) * r1 - diag(conj(diag(r2))) * r2; 
else 
   cmat = conj(ymat' * ymat) / nsamp;    % correlation matrix 
end 

%----- determine number of sources --------------------

 [umat, smat, vmat] = svd(cmat);
 svec  = diag(smat);

hold off, clf, 
set(gcf,'name','Hosa DOA')
stem(svec),  title('svals-doa'), grid on 

while (nsource <= 0 | nsource > msens) 
   txt = ['specify number of sources (order p) : [1,' int2str(msens) ']'];
   nsource = input([txt '  ---> ']);
   if (isempty(nsource)) nsource = 0; end 
end 



% ----------- esprit ------------

Ksens = msens - delta; 
ind1  = 1:Ksens;
ind2  = ind1 + delta; 
uyy   = umat(ind1,1:nsource);
vzz   = umat(ind2,1:nsource);
dvec  = eig( (uyy') * vzz);    % generalized lambda(ryy,ryz) 

dvec = angle(dvec)/(2*pi*dspace*delta); 
bearing = asin(dvec) * 180/pi; 

format compact 
disp(['Number of sensors = ',int2str(msens)])
disp(['Displacement for ESPRIT  = ',int2str(delta),' sensors'])
disp(['bearings (in degrees) estimated by ESPRIT ']) 
disp(bearing) 

%-------- estimate bearing spectra ----------------------------------

prmat = zeros(msens,msens); 
for i=1:msens 
    if (svec(i) > 0) 
       prmat = prmat + vmat(:,i) * vmat(:,i)' / svec(i); 
    end 
end
prvec = prmat(:,1); 

theta = [-90:dtheta:90]' * (pi/180); 
omega = 2*pi*dspace*sin(theta);
mth   = length(theta); 

mlc   = zeros(mth,1);  smus = mlc;   spis = mlc; seig  = mlc; par   = mlc; 
beam  = mlc;  smin = mlc; 

gvec = umat(1,1:nsource).'; 
gmat = umat(2:msens,1:nsource); 
arvec = [1 ; - conj(gmat) * gvec / (1 - gvec'*gvec) ]; 

tmat  = vmat(:,nsource+1:msens);
wvmat = tmat * diag(ones(msens-nsource,1) ./ svec(nsource+1:msens)) * tmat';
tmat  = tmat * tmat'; 

for k =1:mth
  steer   = exp(sqrt(-1) * omega(k) * [0:msens-1]');
  beam(k) = abs(steer' * cmat * steer); 
  smin(k) = 1./ abs(arvec.' * steer)^2; 
  mlc (k) = 1./ abs(steer' * prmat * steer);
  par (k) = 1./ (abs(prvec' * steer)^2); 

  smus(k) = 1./ (steer'  * tmat  * steer); 
  seig(k) = 1./ (steer'  * wvmat * steer); 
  spis(k) = 1./ abs(steer' * vmat(:,msens) )^2; 
end 

%--------- display the angular spectra -----

theta = theta*180/pi;
spec  = abs( [seig, smus, spis, mlc, par, smin, beam] ); 

% ------ normalize to abs-max of unity for display only -----
spmax = max(spec); 
spmax = ones(1,7) ./ spmax;
spec  = spec * diag(spmax);

tname = ['svals of c',int2str(order)]; 
subplot(421), stem(svec),  title(tname), grid on 
subplot(422), plot(theta, 10*log10(spec(:,3))), ylabel('pisar'),grid on
subplot(423), plot(theta, 10*log10(spec(:,1))), ylabel('eig'),grid on
subplot(424), plot(theta, 10*log10(spec(:,2))), ylabel('music'),grid on
subplot(425), plot(theta, 10*log10(spec(:,4))), ylabel('ml'),grid on
subplot(426), plot(theta, 10*log10(spec(:,5))), ylabel('ar'),grid on
subplot(427), plot(theta, 10*log10(spec(:,6))), ylabel('min'),grid on
subplot(428), plot(theta, 10*log10(spec(:,7))), ylabel('beam'),grid on

tname = ['Hosa DOA c' int2str(order)]; 
set(gcf,'name',tname) 

return

