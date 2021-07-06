%% SETUP

clc;
clear;
%close all;

seedPath = 'D:\Downloads\databases\seed\SEED\Preprocessed_EEG\';
fprintf('The path of the database: "%s"\n', seedPath)

% add directories 'utils' and 'libs' to path
% setting the path
rootDir = fileparts(matlab.desktop.editor.getActiveFilename);
eval(['cd ' rootDir])
addpath(genpath([rootDir '/lib']));
addpath(genpath([rootDir '/utils']));

rate = 200;
% Number of samples is different per samples

% Electrodes 62 channels
% todo read csv file
%electrodes = struct()

%% SEED Load Single Analysis
idParticipant = 1;  % range 1-15
idSession = 1;      % range 1-3
idVideo = 3;        % range 1-15
idChannel = 1;      % range 1-62
[fullsignal,bands] = loadAndDecomposeSEED(seedPath,idParticipant,idSession,idVideo,idChannel);
gamma = bands(1); beta= bands(2); alpha = bands(3); theta = bands(4); delta = bands(5);

%% DWT Coefficients
titleDescription = {'Gamma (32-64 Hz)' 'Beta (16-32 Hz)' 'Alpha (8-16 Hz)' 'Theta (4-8 Hz)' 'Delta (0-4 Hz)'};
figure
for i = 1:1:5
    subplot(5,1,i)
    plot(bands(i).coeff)
    axis tight; title(titleDescription(i))
end
sgtitle(['Brain Rythms DWT Coefficients | Label: ' fullsignal.label])

%% Reconstruct coefficients in time domain
step = 1/rate;
titleDescription = {'Original full signal' 'Gamma (32-64 Hz)' 'Beta (16-32 Hz)' 'Alpha (8-16 Hz)' 'Theta (4-8 Hz)' 'Delta (0-4 Hz)'};
figure
for i = 1:1:6
    subplot(6,1,i)
    if i == 1 
        plot(0:step:step*(numel(fullsignal.samples)-1),fullsignal.samples)
    else
        plot(0:step:step*(numel(bands(i-1).samples)-1),bands(i-1).samples)
    end
    axis tight; title(titleDescription(i))
end
sgtitle(['Brain Rythms Time domain | Label: ' fullsignal.label])

%% Bispectrum Direct
clc;
fprintf('Bispectrum Direct started ...\n')

% Available options to pass: 
% 'fullsignal','gamma','beta','alpha','theta','delta' 
signalToTest = fullsignal;

samples = signalToTest.samples;
M = fix(numel(signalToTest.samples)/16);
overlap = 10;
display = 1;

tic
[bispd, waxis] = bispecd(samples,'nfft',0,M,rate,overlap,display);
toc

tic
%[bicod, waxis] = bicoher(samples,'nfft',0,M,0);
toc

fprintf('Bispectrum Direct finished ...\nWaiting for the plots ...\n')