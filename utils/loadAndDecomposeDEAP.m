function [fullsignal,bands] =  loadAndDecomposeDEAP(databasePath,idParticipant,idVideo,idChannel)

% loadAndDecomposeDEAP Loading data and signal decomposition for DEAP
% 
% Description of input data structure:
% https://www.eecs.qmul.ac.uk/mmv/datasets/deap/readme.html
%
% Loading data per pariticipant from the DEAP dataset and decomposing the 
% obtained signal, using wavelets, to frequency bands that correspond to 
% the brain rythms
% 
% Output:
% fullsignal - The obtained preprocessed signal represented as a struct
% with two members: a) samples and b) label which
% corresponds to the samples of the signal and the emotional label.
%
% bands - A 1D matrix containing the brain rythms:
% gamma,beta,alpha,theta,delta. Each rythm is a struct with three members:
% a) coeff, b) samples, c) level (DWT)
%
% NOTES:
% 1. Emotional labels are divided to:
% a) Low  Valence Low  Arousal  ('LVLA')
% b) Low  Valence High Arousal  ('LVHA')
% c) High Valence Low  Arousal  ('HVLA')
% d) High Valence High Arousal  ('HVHA')

if idParticipant<10 
    matName = ['s0' int2str(idParticipant) '.mat'];
else
    matName = ['s'  int2str(idParticipant) '.mat'];
end
clear participant fullsignal bands
participant = load([databasePath matName]);
fullsignal.samples = squeeze(participant.data(idVideo,idChannel,:));
valence = participant.labels(idVideo,1);
arousal = participant.labels(idVideo,2);

% Classify emotional regions
if valence<=5 && arousal<=5 
    emotionQuarter = 'LVLA';
elseif valence<=5 && arousal>=5
    emotionQuarter = 'LVHA';
elseif valence>=5 && arousal<=5
    emotionQuarter = 'HVLA';
elseif valence>=5 && arousal>=5
    emotionQuarter = 'HVHA';
else
    fprintf('ERROR: Invalid emotional region\n')
end
fullsignal.label = emotionQuarter;

% Discrete Wavelet Decomposition
motherWavelet = 'db5';
numLevels = 5;

% Coefficients
%[c,l] = wavedec(X,N,Lo_D,Hi_D) returns the decomposition structure as above, given the low- and high-pass decomposition filters you specify
[c,l] = wavedec(fullsignal.samples,numLevels,motherWavelet);
[gamma.coeff,beta.coeff,alpha.coeff,theta.coeff] = detcoef(c,l,[1 2 3 4]);
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