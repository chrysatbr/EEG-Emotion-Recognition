function [fullsignal,bands] = loadAndDecomposeSEED(databasePath,idParticipant,idSession,idVideo,idChannel)

% loadAndDecomposeSEED Loading data and signal decomposition for SEED
%
% Description of input data structure:
% https://bcmi.sjtu.edu.cn/~seed/seed.html
%
% - 15 participants x 15 experiments (movies) x 3 times repeated the experiment (sesions)
% - 15 experiments (movies) x 3 times per experiment = 45 .mat files per participant
% - Each mat file contains data for these 15 experiments represented as a
% 2d array: 62 (number of channels) x <samples>
% - Number of samples per array isn't contant
% - Each experiment has an emotional tag in the 'label.mat'

clear participant
matFiles = dir([databasePath int2str(idParticipant) '_*.mat']);
sessions = {matFiles(1).name matFiles(2).name matFiles(3).name};

% trying to distinguish regex expressions 'djc_eeg1','djc_eeg11'
if idVideo >= 1  && idVideo <= 5  
    regExp = sprintf('g%d$',idVideo); 
else
    regExp = sprintf('%d$',idVideo);
end

clear participant fullsignal
participant = load([databasePath sessions{idSession}], '-regexp', regExp);
fileName    = string(fieldnames(participant));
participant = struct2cell(participant);
participant = participant{1}; % grab the 2d array
fullsignal.samples = participant(idChannel,:);

fprintf('Participant: %d/15\n',idParticipant)
fprintf('Session:     %d/3\n',idSession)
fprintf('Video:       %d/15\n',idVideo)
fprintf('Channel:     %d/62\n',idChannel)
fprintf('File:        %s\n',sessions{idSession})
fprintf('Array name:  %s\n',fileName)

% reading labels
labels = load([databasePath 'label.mat']);
labels = struct2cell(labels);
labels = labels{1};

idLabel = labels(idVideo);
if idLabel == -1
    emotionQuarter = 'Negative';
elseif idLabel == 0
    emotionQuarter = 'Neutral';
elseif idLabel == 1
    emotionQuarter = 'Positive';
else
    fprintf('Wrong parsed emotion')
    exit(-1)
end
fprintf('Emotion:     %s\n\n',emotionQuarter)

fullsignal.label = emotionQuarter;

% Discrete Wavelet Decomposition
motherWavelet = 'db5';
numLevels = 5;

% Coefficients
% TODO: Frequency rate for SEED is 200 and not 128. The DWT levels should be 
% revised or to find another signal decomposition technique
[c,l] = wavedec(fullsignal.samples,numLevels,motherWavelet);
[gamma.coeff,beta.coeff,alpha.coeff,theta.coeff] = detcoef(c,l,[2 3 4 5]);
delta.coeff = appcoef(c,l,motherWavelet);

delta.level = 5; theta.level = 5; alpha.level = 4; beta.level = 3; gamma.level = 2;

% Reconstructing coefficients to time domain
delta.samples = wrcoef('a',c,l,motherWavelet,delta.level);
theta.samples = wrcoef('d',c,l,motherWavelet,theta.level);
alpha.samples = wrcoef('d',c,l,motherWavelet,alpha.level);
beta.samples  = wrcoef('d',c,l,motherWavelet,beta.level);
gamma.samples = wrcoef('d',c,l,motherWavelet,gamma.level);

bands = [gamma beta alpha theta delta];
end
