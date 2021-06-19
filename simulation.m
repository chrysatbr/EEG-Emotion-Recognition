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
% 5. The original signal: 'fullsignal' and the bands:
% 'delta','theta','alpha','beta','gamma' are structs with three members:
% 1) samples, 2) rate and 3) number of samples (numSamples)

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
fullsignal.rate = 256;
fullsignal.numSamples = 8064;
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

[components,levels] = wavedec(fullsignal.samples,numLevels,motherWavelet);
[gamma.samples,beta.samples,alpha.samples,theta.samples] = detcoef(components,levels,[1 2 3 4]);
delta.samples = appcoef(components,levels,motherWavelet);

gamma.rate = fullsignal.rate/2; gamma.numSamples = numel(gamma.samples);
beta.rate  = gamma.rate/2; beta.numSamples = numel(beta.samples);
alpha.rate = beta.rate/2; alpha.numSamples = numel(alpha.samples);
theta.rate = alpha.rate/2; theta.numSamples = numel(theta.samples);
delta.rate = theta.rate/2; delta.numSamples = numel(delta.samples);

fprintf('Finished decomposition ...\n')

%% DWT PLOTS
figure
subplot(6,1,1)
step = 1/fullsignal.rate; finish = step*(fullsignal.numSamples-1);
plot(0:step:finish,fullsignal.samples)
axis tight; title('Original fullsignal')

subplot(6,1,2)
step = 1/delta.rate; finish = step*(delta.numSamples-1);
plot(0:step:finish,delta.samples)
axis tight; title('Delta (0-4 Hz)')

subplot(6,1,3)
step = 1/theta.rate; finish = step*(theta.numSamples-1);
plot(0:step:finish,theta.samples)
axis tight; title('Theta (4-8 Hz)')

subplot(6,1,4)
step = 1/theta.rate; finish = step*(alpha.numSamples-1);
plot(0:step:finish,alpha.samples)
axis tight; title('Alpha (8-16 Hz)')

subplot(6,1,5)
step = 1/beta.rate; finish = step*(beta.numSamples-1);
plot(0:step:finish,beta.samples)
axis tight; title('Beta (16-32 Hz)')

subplot(6,1,6)
step = 1/beta.rate; finish = step*(gamma.numSamples-1);
plot(0:step:finish,gamma.samples)
axis tight; title('Gamma (32-64 Hz)')
suptitle(['Label: ' emotionQuarter])

%% Bispectrum Direct
clc;
fprintf('Bispectrum Direct started ...\n')

% Available options to pass: 
% 'fullsignal','gamma','beta','alpha','theta','delta' 
signalToTest = fullsignal;

samples = signalToTest.samples;
NFFT = nextpow2(signalToTest.numSamples);
M = min(512,signalToTest.numSamples);
rate = signalToTest.rate;
overlap = 10;
display = 1;

tic
[bispd, waxis] = bispecd(samples,NFFT,0,M,rate,overlap,display);
toc

tic
[bicod, waxis] = bicoher(samples,NFFT,0,M,0);
toc

fprintf('Bispectrum Direct finished ...\nWaiting for the plots ...\n')