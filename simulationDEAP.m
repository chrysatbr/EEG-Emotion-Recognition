%% SETUP

clc;
clear;
%close all;

% Replace the path of the database with your local path
deapPath = 'D:\Downloads\databases\deap\data_preprocessed_matlab\';
fprintf('The path of the database: "%s"\n', deapPath)

% add directoris 'utils' and 'libs' to path
% setting the path
rootDir = fileparts(matlab.desktop.editor.getActiveFilename);
eval(['cd ' rootDir])
addpath(genpath([rootDir '/lib']));
addpath(genpath([rootDir '/utils']));

rate = 128;
numSamples = 8064;

% Electrodes
electrodes = struct('Fp1',1,'AF3',2,'F3',3,'F7',4,'FC5',5,'FC1',6,'C3',7,'T7',8,'CP5',9, ...
    'CP1',10,'P3',11,'P7',12,'PO3',13,'O1',14,'Oz',15,'Pz',16,'Fp2',17,'AF4',18, ...
    'Fz',19,'F4',20,'F8',21,'FC6',22,'FC2',23,'Cz',24,'C4',25,'T8',26,'CP6',27, ...
    'CP2',28, 'P4',29,'P8',30,'PO4',31,'O2',32,'hEOG',33,'vEOG',34,'zEMG',35, ...
    'tEMG',36,'GSR',37,'RESP',38,'BVP',39,'HST',40);
electrodesString = string(fieldnames(electrodes));

%% DEAP Load Single Analysis
% Choose your data! Participant,video,electrode (channel)
idParticipant = 1;
idVideo = 1;
% to change electrodes, change the member of the struct e.g. 'electrodes.F8'
idChannel = electrodes.Fp1;

% Choose your data! Participant,video,electrode (channel)
fprintf('Start loading and wavelet decomposition... \n')
[fullsignal,bands] = loadAndDecomposeDEAP(deapPath,idParticipant,idVideo,idChannel);
gamma = bands{1}; beta= bands{2}; alpha = bands{3}; theta = bands{4}; delta = bands{5};
fprintf('Finished loading and decomposing ...\n\n')
fprintf('Participant ID: %d\nVideo ID: %d\nChannel: %s\n',idParticipant,idVideo,electrodesString(idChannel));
fprintf('Emotion Label: %s\n',fullsignal.label)


%% DWT Coefficients
titleDescription = {'Gamma (32-64 Hz)' 'Beta (16-32 Hz)' 'Alpha (8-16 Hz)' 'Theta (4-8 Hz)' 'Delta (0-4 Hz)'};
figure
for i = 1:1:5
    subplot(5,1,i)
    plot(bands{i}.coeff)
    axis tight; title(titleDescription(i))
end
sgtitle(['Brain Rythms DWT Coefficients | Label: ' fullsignal.label])

%% Reconstruct coefficients in time domain
step = 1/rate; finish = step*(numel(fullsignal.samples)-1);
titleDescription = {'Original full signal' 'Gamma (32-64 Hz)' 'Beta (16-32 Hz)' 'Alpha (8-16 Hz)' 'Theta (4-8 Hz)' 'Delta (0-4 Hz)'};
figure
for i = 1:1:6
    subplot(6,1,i)
    if i == 1 
        plot(0:step:finish,fullsignal.samples)
    else
        plot(0:step:finish,bands{i-1}.samples)
    end
    axis tight; title(titleDescription(i))
end
sgtitle(['Brain Rythms Time domain | Label: ' fullsignal.label])

%% Bispectrum Direct
%clc;
fprintf('Bispectrum Direct started ...\n')

% Available options to pass: 
% 'fullsignal','gamma','beta','alpha','theta','delta' 
signalToTest = fullsignal;

samples = signalToTest.samples;
M = fix(numel(samples)/16);
nfft = 2^nextpow2(M);
freqBins = rate/2^(nextpow2(M))
overlap = 10;
display = 1;

tic
[bispd, waxis] = bispecd(samples,nfft,0,M,rate,overlap,display);
toc

%tic
%[bicod, waxis] = bicoher(samples,nfft,0,M,rate,overlap,display);
%toc

% frequency spectrum (todo: welch method)
% nfft = 2^nextpow2(numel(samples))
% psd = abs(fftshift(fft(samples)));
% fshift = [-nfft/2:nfft/2-1]'*rate/nfft;
% figure
% plot(fshift,psd);


fprintf('Bispectrum Direct finished ...\nWaiting for the plots ...\n')

%% Bulk Analysis Per Participant Bispectrum Direct
% 
% NOTES:
% 1. This would probably take a long time! Here we will try to implement a
% brute force visualization. Trying to create a canvas of 40x32 plots,
% plotting the direct bispectrum of the EEG signals (full signal, not rythms)
% for all the channel per participant. The goal is to spot something interesting!
% 
% 2. In the root directory, a "plots" directory will be created that will
% contain 40 figures that each one depict the bispectrum for all the EEG
% channels
%
% 3. Avoid running this section. You could use the above single
% visualization methods.

clc;

idParticipant = 5;
fprintf('Bulk visualization per participant with id: %d ... \n',idParticipant)

% create folders to save plots
eval('mkdir plots')
participantDir = sprintf('participant_%d/bispectrum',idParticipant);
eval(['mkdir plots/' participantDir])
eval(['mkdir plots/' participantDir '/HVHA'])
eval(['mkdir plots/' participantDir '/LVHA'])
eval(['mkdir plots/' participantDir '/HVLA'])
eval(['mkdir plots/' participantDir '/LVLA'])

M = fix(numSamples/16);
freqBins = rate/2^(nextpow2(M))
overlap = 10

numVideo = 40;    % max 40
numChannels = 32; % max 32

tic
for idVideo = 1:numVideo
    %f = figure;
    f = figure('visible','off');
    set(gcf, 'Position', get(0, 'Screensize')); % maximize figure
    count = 1;
    fprintf('Channels:\n')
    for idChannel = 1:numChannels
        fprintf('%d,',idChannel)
        [fullsignal,bands] = loadAndDecomposeDEAP(deapPath,idParticipant,idVideo,idChannel);
        emotionLabel = fullsignal.label;
        clear waxis bands
        [bisp, waxis,zeroPos] = bispecd(fullsignal.samples,'nfft',0,M,rate,overlap,0);
        clear fullsignal magnBisp
        magnBisp = abs(bisp(zeroPos:end,zeroPos:end));
        clear bisp
        
        % plot
        subplot(4,8,count)
        meshTarget = mesh(waxis(zeroPos:end),waxis(zeroPos:end),magnBisp); grid on
        
        % find peak
        [row,col,value] = maxMatrix(magnBisp);
        x = waxis(row+zeroPos-1);
        y = waxis(col+zeroPos-1);
        z = magnBisp(row,col);
        
        titleHandler = title(sprintf('%d,%s(%d)\n(%0.2f,%0.2f,%0.2f)', ...
            idVideo,electrodesString(idChannel),idChannel,x,y,z));
        titleHandler.FontSize = 8;
        
        count = count+1;
    end
    fprintf('\nVideo ID: %d | Emotion Label: %s\n',idVideo,emotionLabel)
    text = sprintf('Participant ID: %d | Video ID: %d | Emotion Label: %s\n',idParticipant,idVideo,emotionLabel);
    sgtitle(text);
    
    fileName = sprintf('bispecd_video%d',idVideo);
    fullpath = ['plots/' participantDir '/' emotionLabel '/' fileName];
    print(f,fullpath,'-dpng','-r150')
    clf(f)
end
toc

%% Bulk visualization per channel


%% Histogram coupling frequencies per channel
clc;

M = fix(numSamples/16);
freqBins = rate/2^(nextpow2(M))
overlap = 10

numParticipants = 32;
numChannels = 32; % max 32
numVideo = 40;    % max 40
numLabels = 4;

% keep track number of occurences of coupling bispectrum frequencies (f1,f2)
maxFrequencies = numVideo*numParticipants; % max possible
freqChannel = -1. .* ones(numChannels,numLabels,2,maxFrequencies); % preallocate
count = zeros(numChannels,numLabels);

for idParticipant = 1:numParticipants
    fprintf('Participant ID: %d\n',idParticipant)
    for idVideo = 1:numVideo 
        for idChannel = 21:numChannels
            fprintf('%d,',idChannel)
            [fullsignal,bands] = loadAndDecomposeDEAP(deapPath,idParticipant,idVideo,idChannel);
            emotionLabel = fullsignal.label;
            clear waxis bands
            [bisp, waxis,zeroPos] = bispecd(fullsignal.samples,'nfft',0,M,rate,overlap,0);
            clear fullsignal magnBisp
            magnBisp = abs(bisp(zeroPos:end,zeroPos:end));
            clear bisp

            % find peak
            [row,col,value] = maxMatrix(magnBisp);
            x = waxis(row+zeroPos-1);
            y = waxis(col+zeroPos-1);
            z = magnBisp(row,col);

            if emotionLabel == "HVHA"
                label = 1;
            elseif emotionLabel == "LVHA"
                label = 2;
            elseif emotionLabel == "HVLA"
                label = 3;
            elseif emotionLabel == "LVLA"
                label = 4;
            end

            % count occurences of coupling frequencies
            count(idChannel,label) = count(idChannel,label) + 1;
            iter = count(idChannel,label);
            freqChannel(idChannel,label,1,iter) = x;
            freqChannel(idChannel,label,2,iter) = y;
        end
        fprintf('\nVideo ID: %d | Emotion Label: %s\n',idVideo,emotionLabel)
    end
end

labels = ["HVHA" "LVHA" "HVLA" "LVLA"];
eval('mkdir plots')
eval('mkdir plots/histo')

% calculate histogram and export plots
for idChannel = 21:numChannels
    channelName = sprintf('channel_%d',idChannel);
    eval(['mkdir plots/histo/' channelName])
    for idLabel = 1:4
        % remove redudant preallocated frequencies
        k = find(freqChannel(idChannel,idLabel,1,:) == -1,1);

        %f = figure;
        f = figure('visible','off');
        set(gcf, 'Position', get(0, 'Screensize')); % maximize figure

        x = squeeze(freqChannel(idChannel,idLabel,1,1:k-1));
        y = squeeze(freqChannel(idChannel,idLabel,2,1:k-1));
        freqRange = 0:0.25:15;
        hist3([x, y],'CDataMode','auto','Ctrs',{freqRange freqRange});
        h = hist3([x, y],'CDataMode','auto','Ctrs',{freqRange freqRange});
        %histogram2(x,y,0:0.25:15,0:0.25:15,'DisplayStyle','tile');
        
        % get the most frequent element
        [row,col,value] = maxMatrix(h);
        f1 = (row-1)*0.25;
        f2 = (col-1)*0.25;
        
        title(labels(idLabel))
        xticks(sort([0:1:15 f1]));yticks(sort([0:1:15 f2]))
        xlabel('f1');ylabel('f2')
        colorbar
        view(2)

        %print(fullpath,'-dpng','-r150')
        fullpath = sprintf('plots/histo/%s/%s.png',channelName,labels(idLabel));
        exportgraphics(gcf,fullpath,'ContentType','image')
        clf(f)
    end
end