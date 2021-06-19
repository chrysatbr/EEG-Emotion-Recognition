% Higher-Order Spectral Analysis Toolbox.
% Version 2.0.3 (R12 compliant)  27 Dec 2000
%
% New Features.
%   Readme      - Important release information about the HOSA Toolbox 
%                 (double click on Readme, or type "whatsnew directoryname", 
%                  i.e. "whatsnew hosa" to display this file).
%
% Higher-Order Spectrum Estimation: conventional methods 
%   cum2x      - Estimates cross-covariance
%   cum3x      - Estimates third-order cross-cumulants
%   cum4x      - Estimates fourth-order cross-cumulants 
%   cumest     - Estimates auto-cumulants, orders two, three, or four 
%                (cum2est, cum3est, cum4est are sub-ordinate routines) 
%   bicoher    - Estimates bicoherence, direct method 
%   bicoherx   - Estimates cross-bicoherence, direct method 
%   bispecd    - Bispectrum estimation (direct method) 
%   bispecdx   - Cross-Bispectrum estimation (direct method) 
%   bispeci    - Bispectrum estimation (indirect method) 
%   glstat     - Gaussianity-Linearity detection statistics (Hinich test) 
%
% Higher-Order Spectrum Estimation: parametric methods
%   armaqs     - Estimates ARMA parameters via q-slice algorithm 
%   armarts    - Estimates ARMA parameters via residual time-series algorithm 
%   armasyn    - Generates ARMA synthetics 
%   arorder    - Estimates AR order 
%   arrcest    - Estimates AR parameters using correlation &/or cumulants 
%   bispect    - Theoretical bispectrum of an ARMA process 
%   cumtrue    - Computes theoretical (true) cumulants of ARMA processes 
%   maest      - Estimates MA parameters (GM algorithm) 
%   maorder    - Estimates MA order 
%   rpiid      - Generates a sequence of i.i.d. random variables, various p.d.f.'s
%   trispect   - Computes 2-D slice of true trispectrum of ARMA process 
%
% Quadratic Phase Coupling (QPC) 
%   qpcgen     - Generates quadratically-phase coupled harmonics in noise 
%   qpctor     - Detection of quadratic phase coupling via the TOR method 
%   
% Second-Order Volterra Systems
%   nlgen      - Computes the output of a second-order Volterra system 
%   nlpow      - Power's method for parameters of 2nd-order Volterra system
%   nltick     - Tick 's method for parameters of 2nd-order Volterra system
%
% Harmonic Retrieval 
%   harmest    - Estimates frequencies of harmonics  
%   harmgen    - Generates harmonics in Gaussian (colored) noise 
%
% Time-Delay Estimation (TDE) 
%   tde        - Time-delay estimation using third-order cross cumulants 
%   tdeb       - Time-delay estimation using third-order cross bispectrum
%   tdegen     - Generates synthetics for time-delay estimation 
%   tder       - Time-delay estimation using cross-correlation 
% 
% Array Processing - Direction of Arrival (DOA) estimation 
%   doa        - Estimates number of sources and their bearings (cum2 or cum4) 
%   doagen     - Synthetic generator for the DOA problem 
%
% Adaptive Linear Prediction 
%   ivcal      - Computes instrumental variables  (used by rivtr and rivdl)
%   rivdl      - Recursive instrumental variable algorithm: double-lattice filter
%   rivtr      - Recursive instrumental variable algorithm: transversal filter 
%
% Impulse Response (IR), Magnitude and Phase Retrieval
%   biceps     - Estimates impulse response via the bicepstrum method 
%   bicepsf    - Estimates impulse response via the bicepstrum (FFT) method 
%   matul      - Magnitude and phase retrieval (Matsuoka-Ulrych algorithm) 
%
% Time-Frequency Distributions
%   wig2       - Wigner spectrum
%   wig2c      - Wigner spectrum, with Choi-Williams type filtering
%   wig3       - Wigner bispectrum, diagonal slice 
%   wig3c      - Wigner bispectrum, with Choi-Williams type filtering
%   wig4       - Wigner trispectrum, diagonal slice
%   wig4c      - Wigner trispectrum, with Choi-Williams type filtering
%
% Utilities 
%   hprony     - Prony's method for modeling transients 
%   pickpeak   - Picks peaks subject to a separation condition 
%   tls        - Total Least Squares Solution to a set of linear equations 
%   trench     - Trench recursion for non-symmetric Toeplitz matrices
%
% Demos and Quick Help
%   hosahelp   - One-line synopsis of all HOSA mfiles 
%   hosademo   - A demo of the HOSA Toolbox
%

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.15 $
%  A. Swami   October 9, 1996 


%     RESTRICTED RIGHTS LEGEND
% Use, duplication, or disclosure by the Government is subject to
% restrictions as set forth in subparagraph (c) (1) (ii) of the 
% Rights in Technical Data and Computer Software clause of DFARS
% 252.227-7013. 
% Manufacturer: United Signals & Systems, Inc., P.O. Box 2374, 
% Culver City, California 90231. 
%
% This material may be reproduced by or for the U.S. Government pursuant 
% to the copyright license under the clause at DFARS 252.227-7013. 


