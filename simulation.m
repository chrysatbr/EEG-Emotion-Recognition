%% VISUAL ANALYSIS
%
% Notes:
% 1. Since we are dealing with large files, instead of reading the samples all
% at once, we can work per participant extracting the features. 
%
% 2. Currently, the analysis is based on visualization plots, so only one
% EEG fullsignal is analyzed. The fullsignal is determined by the variables:
% 'idParticipant', 'idVideo', 'idChannel'
%
% 3. No usage of EEGLAB, TEAP toolboxes
%
% 4. Computing the bands using 'wavedec' and plot the bispectrum ('HOSA' 
% toolbox) for each one of them including the full fullsignal too.
% 
% 5. The original signal: 'fullsignal' is a struct with 1 member:
% 1) samples
% The bands: 'delta','theta','alpha','beta','gamma' are structs with 3 members:
% 1) coefficients, 2) samples in time domain, 3) level (DWT)

clc;
clear;
%close all;

% Replace the path of the database with your local path
databasePath = 'D:\Downloads\databases\deap\data_preprocessed_matlab\';
fprintf('The path of the database: "%s"\n', databasePath)

% setting the path
rootDir = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(rootDir));

numParticipants = 32;
numVideos = 40;
numSignals = 40;
numLabels = 4;
rate = 256;
numSamples = 8064;
fprintf('Finished set up .. \n\n')

% Choose your data! Participant,video,electrode (channel)
idParticipant = 2;
idVideo = 1;
% USAGE: "idChannel.Fp1" gets the Fp1 electrode data
electrodes = struct('Fp1',1,'AF3',2,'F3',3,'F7',4,'FC5',5,'FC1',6,'C3',7,'T7',8,'CP5',9, ...
    'CP1',10,'P3',11,'P7',12,'PO3',13,'O1',14,'Oz',15,'Pz',16,'Fp2',17,'AF4',18, ...
    'Fz',19,'F4',20,'F8',21,'FC6',22,'FC2',23,'Cz',24,'C4',25,'T8',26,'CP6',27, ...
    'CP2',28, 'P4',29,'P8',30,'PO4',31,'O2',32,'hEOG',33,'vEOG',34,'zEMG',35, ...
    'tEMG',36,'GSR',37,'RESP',38,'BVP',39,'HST',40);
electrodesString = string(fieldnames(electrodes));
% to change electrodes, change the member of the struct e.g. electrodes.F8
idChannel = electrodes.Fp1; 

fprintf('Loading participant: %d/32\n',idParticipant)
if idParticipant<10 
    matName = ['s0' int2str(idParticipant) '.mat'];
else
    matName = ['s'  int2str(idParticipant) '.mat'];
end
participant = load([databasePath matName]);
fullsignal.samples = squeeze(participant.data(idVideo,idChannel,:));
valence = participant.labels(idVideo,1);
arousal = participant.labels(idVideo,2);
fprintf("Video ID: %d\nChannel: %s\n",idVideo,electrodesString(idChannel));

% Classify emotional regions
if valence<5 && arousal<5 
    emotionQuarter = 'LVLA';
elseif valence<5 && arousal>5
    emotionQuarter = 'LVHA';
elseif valence>5 && arousal<5
    emotionQuarter = 'HVLA';
elseif valence>5 && arousal>5
    emotionQuarter = 'HVHA';
end
fprintf('Emotion Label: %s\n',emotionQuarter)
fprintf('Finished loading ...\n')

% Discrete Wavelet Decomposition
motherWavelet = 'db5';
numLevels = 5;

%[C,L] = wavedec(X,N,Lo_D,Hi_D) returns the decomposition structure as above, given the low- and high-pass decomposition filters you specify

[components,levels] = wavedec(fullsignal.samples,numLevels,motherWavelet);
[gamma.coeff,beta.coeff,alpha.coeff,theta.coeff] = detcoef(components,levels,[1 2 3 4]);
delta.coeff = appcoef(components,levels,motherWavelet);

delta.level = 5; theta.level = 5; alpha.level = 4; beta.level = 3; gamma.level = 2;
fprintf('Finished decomposition ...\n')

%% DWT Coefficients
figure
subplot(5,1,1)
plot(delta.coeff)
axis tight; title('Delta (0-4 Hz)')

subplot(5,1,2)
plot(theta.coeff)
axis tight; title('Theta (4-8 Hz)')

subplot(5,1,3)
plot(alpha.coeff)
axis tight; title('Alpha (8-16 Hz)')

subplot(5,1,4)
plot(beta.coeff)
axis tight; title('Beta (16-32 Hz)')

subplot(5,1,5)
plot(gamma.coeff)
axis tight; title('Gamma (32-64 Hz)')
suptitle(['Coefficients | Label: ' emotionQuarter])

%% Reconstruct coefficients in time domain
fprintf('Reconstructing wavelet coefficients ...\n')
step = 1/rate; finish = step*(numSamples-1);

figure
subplot(6,1,1)
plot(0:step:finish,fullsignal.samples)
axis tight; title('Original fullsignal')

subplot(6,1,2)
delta.samples = wrcoef('a',components,levels,motherWavelet,delta.level);
plot(0:step:finish,delta.samples)
axis tight; title('Delta (0-4 Hz)')

subplot(6,1,3)
theta.samples = wrcoef('d',components,levels,motherWavelet,theta.level);
plot(0:step:finish,theta.samples)
axis tight; title('Theta (4-8 Hz)')

subplot(6,1,4)
alpha.samples = wrcoef('d',components,levels,motherWavelet,alpha.level);
plot(0:step:finish,alpha.samples)
axis tight; title('Alpha (8-16 Hz)')

subplot(6,1,5)
beta.samples = wrcoef('d',components,levels,motherWavelet,beta.level);
plot(0:step:finish,beta.samples)
axis tight; title('Beta (16-32 Hz)')

subplot(6,1,6)
gamma.samples = wrcoef('d',components,levels,motherWavelet,gamma.level);
plot(0:step:finish,gamma.samples)
axis tight; title('Gamma (32-64 Hz)')
suptitle(['Time domain | Label: ' emotionQuarter])
fprintf('Finished reconstructing .. \n')

%% Bispectrum Direct
clc;
fprintf('Bispectrum Direct started ...\n')

% Available options to pass: 
% 'fullsignal','gamma','beta','alpha','theta','delta' 
signalToTest = fullsignal;

samples = signalToTest.samples;
M = fix(numSamples/16);
overlap = 10;
display = 1;

tic
[bispd, waxis] = bispecd(samples,'nfft',0,M,rate,overlap,display);
toc

tic
[bicod, waxis] = bicoher(samples,'nfft',0,M,0);
toc

fprintf('Bispectrum Direct finished ...\nWaiting for the plots ...\n')
