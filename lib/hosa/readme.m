% README file for the HOSA Toolbox.
% Version 2.0.3 (R12 Compliant) 27 Dec 2000
%
% Note: There have been no changes in Toolbox functionality.
% *********************************************************************
% Bug fixes for version 2.0.3:
%  
% 1) CUM4X -  Corrected conjugation errors related to the computation of 
%             R_wy, R_zy and M_yx.  
%
% 2) TDE -    Corrected a size error when the 'svdflag' input argument was
%             used.
%
% *********************************************************************
% The HOSA Manual:
%
% The classification example in the Case Studies section of the Higher 
% Order Spectral Analysis Toolbox manual does not define x and y. 
% Users cannot run the example because they do not have the two 
% underwater acoustic signals.  

% *********************************************************************
%   Changes to the Toolbox for version 2.0.2:
%  
%   1) There have been no changes in Toolbox functionality.  Several 
%      of the demo M-files have been modified.  In particular, each demo 
%      can be invoked separately without going through the HOSADEM or 
%      HOSADEMO functions.
%    
%   2) Command-line demos will now plot to only one figure window.
%   
%   3) Case-study demo figures will be closed upon completion of each
%      case study.
%
%   4) The matul function has been changed to correct a bug.  
%      The function now produces the correct coefficient matrix 
%      based on the Matsuoka-Ulrych paper.
%   
%   5) The harmest, doa, tde, and qpctor functions have been modified
%      to validate user-entered order. 
%
%   6) The hprony function has been modified to fix incompatibilities
%      with the toeplitz function.
%
%   7) Several other minor changes (i.e. add grid lines) 
%
% *********************************************************************
% The HOSA Manual:
%
% The classification example in the Case Studies section of the Higher 
% Order Spectral Analysis Toolbox manual does not define x and y. 
% Users cannot run the example because they do not have the two 
% underwater acoustic signals.  

%**********************************************************************
%
% Changes to the HOSA Manual
%  A bug in routine glstat.m has been fixed;  this leads to
%  changes in the output of glstat.m;  HOSA manual pages should 
%  be corrected as shown below (none of the interpretations change)
%
%  On p 1-21,1-22: 
%
%  glstat(g,0.51,256)
% Test statistic for Gaussianity is 22.179 with df = 48, Pfa = 0.9995
% Linearity test:  
% R (estimated) = 0.88819, lambda = 0.68932, R (theory) = 2.9288, N = 14
%
%  glstat(u,0.51,256)
% Test statistic for Gaussianity is 17.4885 with df = 48, Pfa = 1
% Linearity test:  
% R (estimated) = 0.72383, lambda = 0.51704, R (theory) = 2.7453, N = 14
% 
%   glstat(e,0.51,256)
% Test statistic for Gaussianity is 253.3529 with df = 48, Pfa = 0
% Linearity test:  
% R (estimated) = 7.8894, lambda = 9.4555, R (theory) = 8.4655, N = 14
% 
%   glstat(x,0.51,256)
% Test statistic for Gaussianity is 277.5194 with df = 48, Pfa = 0
% Linearity test:  
% R (estimated) = 6.7513, lambda = 10.6519, R (theory) = 8.968, N = 14
% 
%   glstat(z,0.51,256)
% Test statistic for Gaussianity is 12640.0657 with df = 48, Pfa = 0
% Linearity test:  
% R (estimated) = 606.9323, lambda = 492.5759, R (theory) = 59.9088, N = 14
% 
%   glstat(l,0.51,256)
% Test statistic for Gaussianity is 49.931 with df = 48, Pfa = 0.3965
% Linearity test:  
% R (estimated) = 2.6047, lambda = 1.8124, R (theory) = 4.0038, N = 14
% 
%    p 1-96 (sunspot data) 
% 
% Test statistic for Gaussianity is 357.4639 with df = 60, Pfa = 0
% Linearity test:  
% R (estimated) = 14.8592, lambda = 11.0332, R (theory) = 9.1222, N = 16
% 
%    p 1-98 (sunspot data, differenced)
% 
% Test statistic for Gaussianity is 250.1965 with df = 70, Pfa = 0
% Linearity test:  
% R (estimated) = 13.5335, lambda = 6.4449, R (theory) = 7.0433, N = 16
% 
% 
%    p 1-101 (canadian lynx data)
% 
% Test statistic for Gaussianity is 196.752 with df = 28, Pfa = 0
% Linearity test:  
% R (estimated) = 6.8468, lambda = 11.299, R (theory) = 9.2282, N = 5
% 
%    p 1-109 (laughter data)
% 
% Test statistic for Gaussianity is 71.3231 with df = 48, Pfa = 0.0161
% Linearity test:  
% R (estimated) = 2.3216, lambda = 2.376, R (theory) = 4.472, N = 14
%
%**********************************************************************
%
% The classification example in the Case Studies section of the Higher 
% Order Spectral Analysis Toolbox manual does not define x and y. 
% Users cannot run the example because they do not have the two 
% underwater acoustic signals.  
%
%
% The laughter example code in the Case Studies section of the doc
% has an error: Change the two lines following
% % --------------------- power spectra and cum-4 spectra 
% to
% figure(3), [px2,a21,a22] = harmest(sp,25,12,'biased',512,2);
% figure(4), [px4,a41,a42] = harmest(sp,25, 8,'biased',512,4);
%
%
% The example for estimating cumulants in the Polyspectra and Linear 
% Processes section of the manual has a syntax error in its use of the 
% contour command.  The arguments were in the wrong order. Change the 
% line to read as the following and it will work fine.
%   subplot(122), contour(-n:n,-n:n,cmat,8)
%
%**********************************************************************

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.10 $
%  A. Swami  November 21, 1997

% RESTRICTED RIGHTS LEGEND
% Use, duplication, or disclosure by the Government is subject to
% restrictions as set forth in subparagraph (c) (1) (ii) of the
% Rights in Technical Data and Computer Software clause of DFARS
% 252.227-7013.
% Manufacturer: United Signals & Systems, Inc., P.O. Box 2374,
% Culver City, California 90231.
%
%  This material may be reproduced by or for the U.S. Government pursuant
%  to the copyright license under the clause at DFARS 252.227-7013.

disp('HOSA Toolbox Version 2.0.3 (R12 compliant) 27 Dec 2000')
disp('Press any key to see readme file'),pause
clc, help readme
