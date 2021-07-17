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

% Electrodes
electrodes = struct('Fp1',1,'AF3',2,'F3',3,'F7',4,'FC5',5,'FC1',6,'C3',7,'T7',8,'CP5',9, ...
    'CP1',10,'P3',11,'P7',12,'PO3',13,'O1',14,'Oz',15,'Pz',16,'Fp2',17,'AF4',18, ...
    'Fz',19,'F4',20,'F8',21,'FC6',22,'FC2',23,'Cz',24,'C4',25,'T8',26,'CP6',27, ...
    'CP2',28, 'P4',29,'P8',30,'PO4',31,'O2',32,'hEOG',33,'vEOG',34,'zEMG',35, ...
    'tEMG',36,'GSR',37,'RESP',38,'BVP',39,'HST',40);
electrodesString = string(fieldnames(electrodes));

%% DEAP Load Single Analysis
% Choose your data! Participant,video,electrode (channel)
idParticipant = 4;
idVideo = 2;
% to change electrodes, change the member of the struct e.g. 'electrodes.F8'
idChannel = 1;

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

% Reconstruct coefficients in time domain
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
signalToTest = gamma;

samples = signalToTest.samples;
M = fix(numel(samples)/16);
nfft = 2^nextpow2(M);
freqBins = rate/2^(nextpow2(M))
overlap = 10;
display = 0;

tic
[bisp, waxis,zeroPos] = bispecd(samples,nfft,0,M,rate,overlap,display);
magnBisp = abs(bisp(zeroPos:end,zeroPos:end));
figure
mesh(waxis(zeroPos:end),waxis(zeroPos:end),magnBisp); grid on
figure
[bispFilt,peakInfo] = bispecdFilter(magnBisp,freqBins);
mesh(waxis(zeroPos:end),waxis(zeroPos:end),bispFilt)

offsetStart = fix(38/freqBins);
offsetEnd   = fix(46/freqBins);
start = zeroPos + offsetStart;
finish = zeroPos + offsetEnd;

figure
contour(waxis(start:finish),waxis(start:finish),bispFilt(offsetStart:offsetEnd,offsetStart:offsetEnd),8); grid on
hline = refline(1,0);
hline.Color = 'r';

figure
contour(waxis(start:finish),waxis(start:finish),magnBisp(offsetStart:offsetEnd,offsetStart:offsetEnd),8); grid on
hline = refline(1,0);
hline.Color = 'r';
toc

for i = 1:numel(peakInfo)
   fprintf('f1 = %0.2f, f2 = %0.2f, z = %0.2f\n',peakInfo(i).f1,peakInfo(i).f2,peakInfo(i).value);
end

%% PROBLEM: expected normalized bispectrum bounded 0-1 z axis but there is frequency incosistency with bispectrum 

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

%% Bulk Analysis Per Participant Bispectrum
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

idParticipant = 1;
fprintf('Bulk visualization per participant with id: %d ... \n',idParticipant)

% create folders to save plots
eval('mkdir plots')
participantDir = sprintf('participant_%d/gamma',idParticipant);
eval(['mkdir plots/' participantDir])
eval(['mkdir plots/' participantDir '/HVHA'])
eval(['mkdir plots/' participantDir '/LVHA'])
eval(['mkdir plots/' participantDir '/HVLA'])
eval(['mkdir plots/' participantDir '/LVLA'])

overlap = 10

numVideo = 40;    % max 40
numChannels = 32; % max 32

tic
for idVideo = 1:numVideo
    %f = figure;
    f = figure('visible','off');
    set(gcf, 'Position', get(0, 'Screensize')); % maximize figure
    countSubplot = 1;
    fprintf('Channels:\n')
    for idChannel = 1:numChannels
        fprintf('%d,',idChannel)
        [fullsignal,bands] = loadAndDecomposeDEAP(deapPath,idParticipant,idVideo,idChannel);
        gamma = bands{1};
        emotionLabel = fullsignal.label;
        clear bands
        
        signalToTest = gamma;
        samples = signalToTest.samples;
        M = fix(numel(samples)/16);
        nfft = 2^nextpow2(M);
        freqBins = rate/nfft;
        [bisp,waxis,zeroPos] = bispecd(signalToTest.samples,nfft,0,M,rate,overlap,0);
        clear fullsignal magnBisp
        
        magnBisp = abs(bisp(zeroPos:end,zeroPos:end));
        clear bisp
        
        % plot
        subplot(4,8,countSubplot)
        meshTarget = mesh(waxis(zeroPos:end),waxis(zeroPos:end),magnBisp); grid on
        
        % find peak
        [row,col,value] = maxMatrix(magnBisp);
        x = waxis(row+zeroPos-1);
        y = waxis(col+zeroPos-1);
        z = magnBisp(row,col);
        
        titleHandler = title(sprintf('%d,%s(%d)\n(%0.2f,%0.2f,%0.2f)', ...
            idVideo,electrodesString(idChannel),idChannel,x,y,z));
        titleHandler.FontSize = 8;
        
        countSubplot = countSubplot+1;
    end
    fprintf('\nVideo ID: %d | Emotion Label: %s\n',idVideo,emotionLabel)
    
    % save plots
    text = sprintf('Participant ID: %d | Video ID: %d | Emotion Label: %s\n',idParticipant,idVideo,emotionLabel);
    sgtitle(text);
    fullpath = sprintf('plots/%s/%s/bispecd_video%d',participantDir,emotionLabel,idVideo);
    print(f,fullpath,'-dpng','-r150')
    clf(f)
end
toc

%% Histogram coupling frequencies and Bispectrum plots Per Channel
% We want to extract results about peaks (concentrated or diffuse)
% There is no difference in the frequency pair observed in coupling (?)
% for the fullsignal

clc;

overlap = 10

numParticipants = 5;
numChannels = 1;  % max 32
numVideo = 5;     % max 40
numLabels = 4;

stringLabels = ["HVHA" "LVHA" "HVLA" "LVLA"];

% FP1 FP2 AF3 AF4  F7 T7 P7 O1 0z O2 P8 CP6 T8 F8 Cz
%channelRange = [1 2 4 8 12 14 15 17 18 21 24 26 27 30 32]
channelRange = [1];

% create folders to save plots
eval('mkdir plots')

% keep track number of occurences of coupling bispectrum frequencies (f1,f2)
% freqChannel: numChannels x numLabels x number of freq (f1,f2) x ?
% Last dimension size isn't fixed, because the number of peaks is unknown
freqPairPeaks = ones(numChannels,numLabels,2,2); % preallocate
countFreqPair = zeros(numChannels,numLabels);

% frequency of number of peaks
% Last dimension counts the number of peaks
maxPossible = numParticipants*numVideo;
numPeaks = zeros(numChannels,numLabels,maxPossible);
countNumPeaks = zeros(numChannels,numLabels);

% average bispectrum per channel
bispSizeX = 256; % M dependent
bispSizeY = 256;
bispAverage = zeros(numLabels,maxPossible,bispSizeX,bispSizeY); % preallocate
countBispAverage = zeros(numLabels);

for idChannel = channelRange    
    channelDir = sprintf('bispectrum/channel_%d/gamma',idChannel);
    eval(['mkdir plots/' channelDir])
    eval(['mkdir plots/' channelDir '/HVHA'])
    eval(['mkdir plots/' channelDir '/LVHA'])
    eval(['mkdir plots/' channelDir '/HVLA'])
    eval(['mkdir plots/' channelDir '/LVLA'])
    
    fprintf('Channel ID: %d\n',idChannel)
    for idVideo = 1:numVideo
        f = figure('visible','off');
        set(gcf, 'Position', get(0, 'Screensize')); % maximize figure
        countSubplot = 1; % keep track subplots
        for idParticipant = 1:numParticipants
            fprintf('%d,',idParticipant)
            [fullsignal,bands] = loadAndDecomposeDEAP(deapPath,idParticipant,idVideo,idChannel);
            gamma = bands{1};
            emotionLabel = fullsignal.label;
            clear bands
            
            signalToTest = gamma;
            samples = signalToTest.samples;
            M = fix(numel(samples)/16);
            nfft = 2^nextpow2(M);
            freqBins = rate/nfft;
            [bisp,waxis,zeroPos] = bispecd(signalToTest.samples,nfft,0,M,rate,overlap,0);
            clear fullsignal magnBisp
            
            magnBisp = abs(bisp(zeroPos:end,zeroPos:end));
            clear bisp

            % filtering bispectrum
            [bispFilt,peaksInfo] = bispecdFilter(magnBisp,freqBins);            
            
            % bispectrum plot
            subplot(4,8,countSubplot)
            %mesh(waxis(zeroPos:end),waxis(zeroPos:end),magnBisp); grid on
            offsetStart = fix(38/freqBins);
            offsetEnd   = fix(46/freqBins);
            start = zeroPos + offsetStart;
            finish = zeroPos + offsetEnd;
            %contour(waxis(start:finish),waxis(start:finish),magnBisp(offsetStart:offsetEnd,offsetStart:offsetEnd),12); grid on
            %contour(waxis(start:finish),waxis(start:finish),bispFilt(offsetStart:offsetEnd,offsetStart:offsetEnd),12); grid on
            hline = refline(1,0);
            hline.Color = 'r';
            
            % find peak
            %[row,col,value] = maxMatrix(magnBisp);
            %x = waxis(row+zeroPos-1);
            %y = waxis(col+zeroPos-1);
            %z = magnBisp(row,col);

            titleHandler = title(sprintf('%s,%d\n(%0.2f,%0.2f,%0.2f)', ...
                emotionLabel,idParticipant,peaksInfo(1).f1,peaksInfo(1).f2,peaksInfo(1).value));
            titleHandler.FontSize = 8;
        
            if emotionLabel == "HVHA"
                label = 1;
            elseif emotionLabel == "LVHA"
                label = 2;
            elseif emotionLabel == "HVLA"
                label = 3;
            elseif emotionLabel == "LVLA"
                label = 4;
            end

            % count number of peaks
            countNumPeaks(idChannel,label) = countNumPeaks(idChannel,label) + 1;
            numPeaks(idChannel,label,countNumPeaks(idChannel,label)) = numel(peaksInfo);
            
            % count occurences of coupling frequencies
            for i = 1:numel(peaksInfo)
                countFreqPair(idChannel,label) = countFreqPair(idChannel,label) + 1;
                iter = countFreqPair(idChannel,label);
                %fprintf('Iter %d, f1 = %0.2f, f2 = %0.2f, value = %0.2f\n',i,peaksInfo(i).f1,peaksInfo(i).f2,peaksInfo(i).value)
                freqPairPeaks(idChannel,label,1,iter) = peaksInfo(i).f1;
                freqPairPeaks(idChannel,label,2,iter) = peaksInfo(i).f2;
            end
            
            % bispectrum average
            countBispAverage(label) = countBispAverage(label) + 1;
            % normalize with max value to avoid losing peaks with small z-axis value
            bispAverage(label,countBispAverage(label),:,:) = magnBisp./max(magnBisp(:));

            countSubplot = countSubplot + 1;
        end
        fprintf('\nVideo ID: %d\n',idVideo)

        % save bispectrum plots
        text = sprintf('Channel: %d (%s) | Video ID: %d\n',...
            idChannel,electrodesString(idChannel),idVideo);
        sgtitle(text);
        fullpath = sprintf('plots/bispectrum/channel_%d/gamma/bispecd_contour_video%d_filt', ...
                idChannel,idVideo);
        %print(f,fullpath,'-dpng','-r150')
        clf(f)
    end
        
    % narrow the plot range for gamma band
    offsetStart = fix(38/freqBins);
    offsetEnd   = fix(46/freqBins);
    start = zeroPos + offsetStart;
    finish = zeroPos + offsetEnd;
    
    % export contours per emotion
%     for i = 1:4
%         for j = 1:countBispAverage(i)
%             f = figure('visible','off');
%             set(gcf, 'Position', get(0, 'Screensize')); % maximize figure
% 
%             bisp = squeeze(bispAverage(i,j,:,:));
%             contour(waxis(start:finish),waxis(start:finish),...
%             bisp(offsetStart:offsetEnd,offsetStart:offsetEnd),12); grid on
%             hline = refline(1,0); hline.Color = 'r';
% 
%             fullpath = sprintf('plots/bispectrum/channel_%d/gamma/%s/bispecd_contour%d', ...
%                 idChannel,stringLabels(i),j);
%             print(f,fullpath,'-dpng','-r150')
%             clf(f)
%         end
%     end
    
    countLabel = 1;
    figure
    for i = 1:4
        fprintf('Emotion %s: %d samples\n',stringLabels(i),countBispAverage(i))
        bispAvgPerLabel = squeeze(mean(squeeze(bispAverage(i,1:countBispAverage(i),:,:)),1));
        
        subplot(2,2,countLabel)
        contour(waxis(start:finish),waxis(start:finish),...
            bispAvgPerLabel(offsetStart:offsetEnd,offsetStart:offsetEnd),12); grid on
        hline = refline(1,0); hline.Color = 'r';
        title(stringLabels{i})
        
        countLabel = countLabel+1;
    end
    text = sprintf('Channel %d (%s)',idChannel,electrodesString(idChannel));
    sgtitle(text)
    clear bispAvgPerLabel bispAverage
end

%% Histogram frequency pairs

eval('mkdir plots')
eval('mkdir plots/histo')

% compute histogram and export plots
for idChannel = channelRange
    channelName = sprintf('channel_%d/gamma',idChannel);
    eval(['mkdir plots/histo/' channelName])
    for idLabel = 1:4
        %f = figure;
        f = figure('visible','off');
        set(gcf, 'Position', get(0, 'Screensize')); % maximize figure

        % remove zeros due to dynamic allocation
        x = nonzeros(squeeze(freqPairPeaks(idChannel,idLabel,1,:)));
        y = nonzeros(squeeze(freqPairPeaks(idChannel,idLabel,2,:)));
        %count = 1;
        %for i = 1:numel(nonzeros(numPeaks(idChannel,idLabel,:)))
            %for j = 1:nonzeros(numPeaks(idChannel,idLabel,i))
                %fprintf('f1 = %0.2f, f2 = %0.2f\n',x(count),y(count));
                %count = count + 1;
            %end
            %fprintf('\n')
        %end
        freqRange = 38:0.5:45;
        hist3([x, y],'CDataMode','auto','Ctrs',{freqRange freqRange});
        h = hist3([x, y],'CDataMode','auto','Ctrs',{freqRange freqRange});
        %histogram2(x,y,0:0.25:15,0:0.25:15,'DisplayStyle','tile');
        
        % get the most frequent element
        [row,col,value] = maxMatrix(h);
        f1 = (row-1)*freqBins;
        f2 = (col-1)*freqBins;
        
        title(stringLabels(idLabel))
        xticks(freqRange);yticks(freqRange)
        xlabel('f1');ylabel('f2')
        colorbar
        view(2)

        %print(fullpath,'-dpng','-r150')
        fullpath = sprintf('plots/histo/%s/%s.png',channelName,stringLabels(idLabel));
        exportgraphics(gcf,fullpath,'ContentType','image')
        clf(f)
    end
end

% Histogram number of peaks

for idChannel = channelRange
    channelName = sprintf('channel_%d/gamma',idChannel);
    for idLabel = 1:4
        %f = figure;
        f = figure('visible','off');
        set(gcf, 'Position', get(0, 'Screensize')); % maximize figure
        histogram(nonzeros(numPeaks(idChannel,idLabel,:)))
        fullpath = sprintf('plots/histo/%s/%s_peaks.png',channelName,stringLabels(idLabel));
        exportgraphics(gcf,fullpath,'ContentType','image')
        clf(f)
    end
end