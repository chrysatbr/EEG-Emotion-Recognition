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
idParticipant = 10;
idVideo = 4;
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


%% Single Analysis Bispectrum non redundant region (triangle region)
% Run first the "DEAP load" section

% testing filters
%a = ones(5);
%triangleFilter = flipud(tril(a)) .* triu(a);
%fillHelper = -1 .* a + triangleFilter;

M = 1024;
nfft = 1024;
freqBins = (rate/2)/(nfft/2-1);
overlap = 50;

fprintf('Bispectrum started...\n')
[bisp,waxis,zeroPos] = bispecd(fullsignal.samples,nfft,0,M,rate,overlap,0);

%maxLag = 512;
%[bisp,waxis,zeroPos] = bispeci(fullsignal.samples,maxLag,M,rate,overlap,'unbiased',nfft,1,0);

magnBisp = abs(bisp(zeroPos:end,zeroPos:end));
fprintf('Bispectrum ended...\n')

figure
mesh(waxis(zeroPos:end),waxis(zeroPos:end),magnBisp)
title('Original')

fprintf('Filtering non redundant region...\n')
triangleFilter = flipud(tril(ones(size(magnBisp)))) .* triu(ones(size(magnBisp)));
bispFilt = triangleFilter .* magnBisp;

figure
mesh(waxis(zeroPos:end),waxis(zeroPos:end),bispFilt); grid on
title('Non redundant region')

% find bands in the non redudant triangle region
fprintf('Finding bands...\n')

% Available options to pass: 
% 'fullsignal','gamma','beta','alpha','theta','delta' 
signalToTest = theta;

[x,y] = size(bispFilt);
startIndex  = max(round(signalToTest.start/freqBins),1);
finishIndex = min(round(signalToTest.finish/freqBins,y-1));
xRange = startIndex:finishIndex;
yRange = 1:round(32/freqBins);

startPlotIndex  = max(round(signalToTest.plotStart/freqBins),1);
finishPlotIndex = min(round(signalToTest.plotFinish/freqBins,y-1));
xPlotRange = startPlotIndex:finishPlotIndex; 
yPlotRange = 1:round(min(signalToTest.plotFinish,32)/freqBins);

bispBand = bispFilt(yRange,xRange);
figure
mesh(waxis(zeroPos-1+xRange),waxis(zeroPos-1+yRange),bispBand)
title('Band portion of non redundant triangle region')

[bispPeaks,peakInfo] = findPeaks(bispBand,startIndex,freqBins);
%size(bispBand)
%size(bispPeaks)
figure
mesh(waxis(zeroPos-1+xRange),waxis(zeroPos-1+yRange),bispPeaks); grid on
title('Peaks')

fprintf('Peaks coordinates:\n')
for i = 1:numel(peakInfo)
   fprintf('f1 = %0.2f, f2 = %0.2f, z = %0.2f\n',peakInfo(i).f1,peakInfo(i).f2,peakInfo(i).value);
end

% to detect triangle when dividing the triangle region to bands fill the
% non triangle region with -1

fillHelper = -1 .* ones(size(magnBisp)) + triangleFilter;
bispFiltFilled = bispFilt + fillHelper;

figure
mesh(waxis(zeroPos:end),waxis(zeroPos:end),bispFiltFilled); grid on
title('Non redundant region filled')

figure
step = max(bispFiltFilled(:))/10;
contourf(waxis(zeroPos:end),waxis(zeroPos:end),bispFiltFilled,0:step:max(bispFiltFilled(:)))
colorbar
title('Non redundant region filled')

plotBispBandFilled = bispFiltFilled(yPlotRange,xPlotRange);
figure
mesh(waxis(zeroPos-1+xPlotRange),waxis(zeroPos-1+yPlotRange),plotBispBandFilled)
title('Band portion of non redundant triangle region filled')

figure
step = max(plotBispBandFilled(:))/10;
contourf(waxis(zeroPos-1+xPlotRange),waxis(zeroPos-1+yPlotRange),plotBispBandFilled,0:step:max(plotBispBandFilled(:)))
colorbar
title('Band portion of non redundant triangle region filled')

bispBandFilled = bispFiltFilled(yRange,xRange);
bispFeatureVector = bispBandFilled(bispBandFilled~=-1);
figure
plot(bispFeatureVector)
title('Values of the feature vector')

fprintf('Finished\n')
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

M = 1024;
nfft = 1024;
freqBins = (rate/2)/(nfft/2-1);
overlap = 50;

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

%% Feature extraction

clc;

numParticipants = 32;
numChannels     = 32;  % max 32
numVideo        = 40;     % max 40
numLabels       = 4;

% FP1 FP2 AF3 AF4  F7 T7 P7 O1 0z O2 P8 CP6 T8 F8 Cz
%channelRange = [1 2 4 8 12 14 15 17 18 21 24 26 27 30 32]
channelRange = [1];

stringBands = ["fullsignal" "gamma" "beta" "alpha" "theta"];
numBands = 5;

% create folders to save plots
eval('mkdir plots')

%maxPossible   = numParticipants*numVideo;
maxPossible   = 439; % 439 HVHA, 269 HVLA, 298 LVHA, 274 LVLA
countPerLabel = zeros(numBands,numLabels);

% feature: max values
rawMax        = zeros(numBands,numLabels,maxPossible); % preallocate
normDivMax    = zeros(numBands,numLabels,maxPossible); % normalize by dividing baseline feature
normSignalMax = zeros(numBands,numLabels,maxPossible); % normalize by subtracting mean baseline from emotional state signal

% feature: std values
rawStd        = zeros(numBands,numLabels,maxPossible); % preallocate
normDivStd    = zeros(numBands,numLabels,maxPossible); % normalize by dividing baseline feature
normSignalStd = zeros(numBands,numLabels,maxPossible); % normalize by subtracting mean baseline from emotional state signal

% feature: mean values
rawMean        = zeros(numBands,numLabels,maxPossible); % preallocate
normDivMean    = zeros(numBands,numLabels,maxPossible); % normalize by dividing baseline feature
normSignalMean = zeros(numBands,numLabels,maxPossible); % normalize by subtracting mean baseline from emotional state signal

% feature: skewness
rawSkew         = zeros(numBands,numLabels,maxPossible);
normDivSkew     = zeros(numBands,numLabels,maxPossible);
normSignalSkew  = zeros(numBands,numLabels,maxPossible);

% feature: Kurtosis
rawKurt         = zeros(numBands,numLabels,maxPossible);
normDivKurt     = zeros(numBands,numLabels,maxPossible);
normSignalKurt  = zeros(numBands,numLabels,maxPossible);

% feature: bispectral entropy BE1
%rawBE1       = zeros(numBands,numLabels,maxPossible); 
%normDivBE1   = zeros(numBands,numLabels,maxPossible);
%normSignalE1 = zeros(numBands,numLabels,maxPossible);

% feature: bispectral entropy BE2
%rawBE2       = zeros(numBands,numLabels,maxPossible);
%normDivBE2   = zeros(numBands,numLabels,maxPossible);
%normSignalE2 = zeros(numBands,numLabels,maxPossible);

% feature: differential entropy DE
rawDE         = zeros(numBands,numLabels,maxPossible);
normDivDE     = zeros(numBands,numLabels,maxPossible);
normSignalDE  = zeros(numBands,numLabels,maxPossible);

% feature: number of peaks
numPeaks = zeros(numBands,numLabels,maxPossible);

% averaging bispectrums (magnitude) triangle region
avgBispFullsignal = zeros(numLabels,maxPossible,256,511);
avgBispGamma      = zeros(numLabels,maxPossible,256,256);
avgBispBeta       = zeros(numLabels,maxPossible,256,128);
avgBispAlpha      = zeros(numLabels,maxPossible,256,64);
avgBispTheta      = zeros(numLabels,maxPossible,256,32);

% keep track number of occurences of coupling bispectrum frequencies (f1,f2)
% freqChannel: numChannels x numLabels x number of freq (f1,f2) x ?
% Last dimension size isn't fixed, because the number of peaks is unknown
freqPairPeaks = zeros(numBands,numLabels,2,2);
countFreqPair = zeros(numBands,numLabels);

% emotional state low resolution to match baseline
MLowRes        = 384; % to match freq resolution with the baseline for the normalization (fair play)
nfftLowRes     = 384;
freqBinsLowRes = (rate/2)/(nfftLowRes/2-1);
%maxLagLowRes  = 256;           

% emotional state high resolution
MHighRes        = 1024;
nfftHighRes     = 1024;
freqBinsHighRes = (rate/2)/(nfftHighRes/2-1);
%maxLagHighRes  = 512;

overlap = 50;

for idChannel = channelRange    
    fprintf('Channel ID: %d\n',idChannel)
    for idVideo = 1:numVideo
        for idParticipant = 1:numParticipants
            fprintf('%d,',idParticipant)
            [fullsignal,bands] = loadAndDecomposeDEAP(deapPath,idParticipant,idVideo,idChannel);
            emotionLabel = fullsignal.label;
            
            % BASELINE
            %[bispBaseline,waxis,zeroPos] = bispeci(fullsignal.baseline,maxLagLowRes,MLowRes,overlap,'unbiased',nfftLowRes,1,0);
            [bispBaseline,waxis,zeroPos] = bispecd(fullsignal.baseline,nfftLowRes,0,MLowRes,rate,overlap,0);
            clear magnBispBaseline
            magnBispBaseline = abs(bispBaseline(zeroPos:end,zeroPos:end));
            clear bispBaseline
            
            % EMOTIONAL STATE LOW RES
            %[bisp,waxis,zeroPos] = bispeci(fullsignal.samples,maxLagLowRes,MLowRes,overlap,'unbiased',nfftLowRes,1,0);
            [bisp,waxis,zeroPos] = bispecd(fullsignal.samples,nfftLowRes,0,MLowRes,rate,overlap,0);
            clear magnBispLowRes
            magnBispLowRes = abs(bisp(zeroPos:end,zeroPos:end));
            clear bisp
            
            % EMOTIONAL STATE HIGH RES
            %[bisp,waxis,zeroPos] = bispeci(fullsignal.samples,maxLagHighRes,MHighRes,overlap,'unbiased',nfftHighRes,1,0); 
            [bisp,waxis,zeroPos] = bispecd(fullsignal.samples,nfftHighRes,0,MHighRes,rate,overlap,0);
            clear magnBispHighRes
            magnBispHighRes = abs(bisp(zeroPos:end,zeroPos:end));
            clear bisp
            
            % NORMALIZED SIGNAL
            [bisp,waxis,zeroPos] = bispecd(fullsignal.norm,nfftHighRes,0,MHighRes,rate,overlap,0);
            clear magnBispNorm
            magnBispNorm = abs(bisp(zeroPos:end,zeroPos:end));
            clear bisp
            
            if emotionLabel     == "HVHA"
                label = 1;
            % Low Valence—High Arousal (Surprise 4o)
            elseif emotionLabel == "LVHA"
                label = 2;
            %High Valence—Low Arousal (angry & disgust 2o)
            elseif emotionLabel == "HVLA"
                label = 3;
            %Low Valence—Low Arousal (Fear & Sad 3o)
            elseif emotionLabel == "LVLA"
                label = 4;
            end
            
            % filtering non redudant region
            triangleFilterLowRes = flipud(tril(ones(size(magnBispLowRes)))) .* triu(ones(size(magnBispLowRes)));
            bispFiltBaseline = triangleFilterLowRes .* magnBispBaseline;
            clear magnBispBaseline

            bispFiltLowRes = triangleFilterLowRes .* magnBispLowRes;
            clear magnBispLowRes;
            
            triangleFilterHighRes = flipud(tril(ones(size(magnBispHighRes)))) .* triu(ones(size(magnBispHighRes)));
            bispFiltHighRes = triangleFilterHighRes .* magnBispHighRes;
            clear magnBispHighRes;
            
            bispFiltNorm = triangleFilterHighRes .* magnBispNorm;
            clear magnBispNorm;
            
            % to detect triangle when dividing the triangle region to bands fill the non triangle region with -1
            fillHelper = -1 .* ones(size(bispFiltLowRes)) + triangleFilterLowRes;
            bispFiltBaseline = bispFiltBaseline + fillHelper;
            bispFiltLowRes = bispFiltLowRes + fillHelper;
            
            fillHelper = -1 .* ones(size(bispFiltHighRes)) + triangleFilterHighRes;
            bispFiltHighRes = bispFiltHighRes + fillHelper;
            bispFiltNorm = bispFiltNorm + fillHelper;
            
            % signals to test
            signalToTestMat = {fullsignal bands{1} bands{2} bands{3} bands{4}};
            
            for idBand = 1:numel(signalToTestMat)
                signalToTest = signalToTestMat{idBand};

                % ranges
                [yRangeLowRes,xRangeLowRes,yPlotRange,xPlotRange,startIndexLowRes] ...
                    = rangeBand(bispFiltLowRes,signalToTest,freqBinsLowRes);
                [yRangeHighRes,xRangeHighRes,yPlotRange,xPlotRange,startIndexHighRes] = ...
                    rangeBand(bispFiltHighRes,signalToTest,freqBinsHighRes);

                bispBandFilledBaseline = bispFiltBaseline(yRangeLowRes,xRangeLowRes);
                bispFeatureVectorBaseline = bispBandFilledBaseline(bispBandFilledBaseline~=-1);
                clear bispBandFilledBaseline

                bispBandFilledLowRes = bispFiltLowRes(yRangeLowRes,xRangeLowRes);
                bispFeatureVectorLowRes = bispBandFilledLowRes(bispBandFilledLowRes~=-1);
                clear bispBandFilledLowRes

                bispBandFilledHighRes = bispFiltHighRes(yRangeHighRes,xRangeHighRes);
                bispFeatureVectorHighRes = bispBandFilledHighRes(bispBandFilledHighRes~=-1);
                %clear bispBandFilledHighRes

                bispBandFilledNorm = bispFiltNorm(yRangeHighRes,xRangeHighRes);
                bispFeatureVectorNorm = bispBandFilledNorm(bispBandFilledNorm~=-1);
                clear bispBandFilledNorm

                countPerLabel(idBand,label) = countPerLabel(idBand,label) + 1;

                % feature: number of peaks
                [bispPeaks,peaksInfo] = findPeaks(bispBandFilledHighRes,startIndexHighRes,freqBinsHighRes);
                clear bispPeaks
                numPeaks(idBand,label,countPerLabel(idBand,label)) = numel(peaksInfo);

                % count occurences of coupling frequencies
                for i = 1:numel(peaksInfo)
                    countFreqPair(idBand,label) = countFreqPair(idBand,label) + 1;
                    iter = countFreqPair(idBand,label);
                    freqPairPeaks(idBand,label,1,iter) = peaksInfo(i).f1;
                    freqPairPeaks(idBand,label,2,iter) = peaksInfo(i).f2;
                end

                % averaging bispectrum
                if     idBand == 1
                    avgBispFullsignal(label,countPerLabel(1,label),:,:) = bispBandFilledHighRes;
                elseif idBand == 2
                    avgBispGamma(label,countPerLabel(2,label),:,:)      = bispBandFilledHighRes;
                elseif idBand == 3
                    avgBispBeta(label,countPerLabel(3,label),:,:)       = bispBandFilledHighRes;
                elseif idBand == 4
                    avgBispAlpha(label,countPerLabel(4,label),:,:)      = bispBandFilledHighRes;
                elseif idBand == 5
                    avgBispTheta(label,countPerLabel(5,label),:,:)      = bispBandFilledHighRes;
                end
                
                % feature: max value
                rawMax(idBand,label,countPerLabel(idBand,label)) = max(bispFeatureVectorHighRes);
                normDivMax(idBand,label,countPerLabel(idBand,label)) = max(bispFeatureVectorLowRes)/max(bispFeatureVectorBaseline);
                normSignalMax(idBand,label,countPerLabel(idBand,label)) = max(bispFeatureVectorNorm);

                % feature: mean value
                rawMean(idBand,label,countPerLabel(idBand,label)) = mean(bispFeatureVectorHighRes);
                normDivMean(idBand,label,countPerLabel(idBand,label)) = mean(bispFeatureVectorLowRes)/mean(bispFeatureVectorBaseline);
                normSignalMean(idBand,label,countPerLabel(idBand,label)) = mean(bispFeatureVectorNorm);

                % feature: std value
                rawStd(idBand,label,countPerLabel(idBand,label)) = std(bispFeatureVectorHighRes);
                normDivStd(idBand,label,countPerLabel(idBand,label)) = std(bispFeatureVectorLowRes)/std(bispFeatureVectorBaseline);
                normSignalStd(idBand,label,countPerLabel(idBand,label)) = std(bispFeatureVectorNorm);

                % feature: skewness
                rawSkew(idBand,label,countPerLabel(idBand,label)) = skewness(bispFeatureVectorHighRes,0,'all');
                normDivSkew(idBand,label,countPerLabel(idBand,label)) = skewness(bispFeatureVectorLowRes,0,'all')/skewness(bispFeatureVectorBaseline,0,'all');
                normSignalSkew(idBand,label,countPerLabel(idBand,label)) = skewness(bispFeatureVectorNorm,0,'all');

                % feature: Kurtosis
                rawKurt(idBand,label,countPerLabel(idBand,label)) = kurtosis(bispFeatureVectorHighRes,0,'all');
                normDivKurt(idBand,label,countPerLabel(idBand,label)) = kurtosis(bispFeatureVectorLowRes,0,'all')/kurtosis(bispFeatureVectorBaseline,0,'all');
                normSignalKurt(idBand,label,countPerLabel(idBand,label)) = kurtosis(bispFeatureVectorNorm,0,'all');

                % feature: DE
                DEfunc = @(bisp) 0.5*log(2*pi*exp(1)*var(bisp));
                rawDE(idBand,label,countPerLabel(idBand,label)) = DEfunc(bispFeatureVectorHighRes); 
                normDivDE(idBand,label,countPerLabel(idBand,label)) = DEfunc(bispFeatureVectorLowRes)/DEfunc(bispFeatureVectorBaseline); 
                normSignalDE(idBand,label,countPerLabel(idBand,label)) = DEfunc(bispFeatureVectorNorm); 
            end
        end
        fprintf('\nVideo ID: %d\n',idVideo)
    end
 
    % EXPORT PLOTS
    gamma = bands{1}; beta= bands{2}; alpha = bands{3}; theta = bands{4};
    for idBand = 1:numel(signalToTestMat)
        channelName = sprintf('channel_%d/%s',idChannel,stringBands(idBand));
        eval(['mkdir plots/histo/' channelName])
        eval(['mkdir plots/boxplot/' channelName])
        eval(['mkdir plots/average/' channelName])

        % average
        if     idBand == 1
            startPlotIndex  = max(round(fullsignal.plotStart/freqBinsHighRes),1);
            finishPlotIndex = min(round(fullsignal.plotFinish/freqBinsHighRes,511));
            xPlotRange = startPlotIndex:finishPlotIndex; 
            yPlotRange = 1:round(min(fullsignal.plotFinish,32)/freqBinsHighRes);
            exportBispAverage(avgBispFullsignal,waxis,squeeze(countPerLabel(idBand,:)),yPlotRange,xPlotRange,zeroPos,'avgFullsignal',channelName)
        elseif idBand == 2
            startPlotIndex  = max(round((gamma.start - gamma.plotStart)/freqBinsHighRes),1)
            finishPlotIndex = min(round((gamma.finish - gamma.plotFinish)/freqBinsHighRes),256)
            xPlotRange = startPlotIndex:finishPlotIndex; 
            yPlotRange = 1:round(min(gamma.plotFinish,32)/freqBinsHighRes);
            exportBispAverage(avgBispGamma,waxis,squeeze(countPerLabel(idBand,:)),yPlotRange,xPlotRange,zeroPos,'avgGammasignal',channelName)
        elseif idBand == 3
            startPlotIndex  = max(round((beta.start - beta.plotStart)/freqBinsHighRes),1)
            finishPlotIndex = 128;
            xPlotRange = startPlotIndex:finishPlotIndex; 
            yPlotRange = 1:round(min(beta.plotFinish,32)/freqBinsHighRes);
            exportBispAverage(avgBispBeta,waxis,squeeze(countPerLabel(idBand,:)),yPlotRange,xPlotRange,zeroPos,'avgBetasignal',channelName)
        elseif idBand == 4
            startPlotIndex  = max(round((alpha.start - alpha.plotStart)/freqBinsHighRes),1)
            finishPlotIndex = 64;
            xPlotRange = startPlotIndex:finishPlotIndex; 
            yPlotRange = 1:round(min(alpha.plotFinish,32)/freqBinsHighRes);
            exportBispAverage(avgBispAlpha,waxis,squeeze(countPerLabel(idBand,:)),yPlotRange,xPlotRange,zeroPos,'avgAlphasignal',channelName)
        elseif idBand == 5
            startPlotIndex  = max(round((theta.start - theta.plotStart)/freqBinsHighRes),1)
            finishPlotIndex = 32;
            xPlotRange = startPlotIndex:finishPlotIndex; 
            yPlotRange = 1:round(min(theta.plotFinish,32)/freqBinsHighRes);
            exportBispAverage(avgBispTheta,waxis,squeeze(countPerLabel(idBand,:)),yPlotRange,xPlotRange,zeroPos,'avgThetasignal',channelName)
        end
        
        % HISTO
        exportHisto(squeeze(normDivMax(idBand,:,:)),'normDivMax',channelName)
        exportHisto(squeeze(normSignalMax(idBand,:,:)),'normSignalMax',channelName)
        exportHisto(squeeze(rawMax(idBand,:,:)),'rawMax',channelName)

        exportHisto(squeeze(normDivMean(idBand,:,:)),'normDivMean',channelName)
        exportHisto(squeeze(normSignalMean(idBand,:,:)),'normSignalMean',channelName)
        exportHisto(squeeze(rawMean(idBand,:,:)),'rawMean',channelName)

        exportHisto(squeeze(normDivStd(idBand,:,:)),'normDivStd',channelName)
        exportHisto(squeeze(normSignalStd(idBand,:,:)),'normSignalStd',channelName)
        exportHisto(squeeze(rawStd(idBand,:,:)),'rawStd',channelName)

        exportHisto(squeeze(normDivSkew(idBand,:,:)),'normDivSkew',channelName)
        exportHisto(squeeze(normSignalSkew(idBand,:,:)),'normSignalSkew',channelName)
        exportHisto(squeeze(rawSkew(idBand,:,:)),'rawSkew',channelName);

        exportHisto(squeeze(normDivKurt(idBand,:,:)),'normDivKurt',channelName)
        exportHisto(squeeze(normSignalKurt(idBand,:,:)),'normSignalKurt',channelName)
        exportHisto(squeeze(rawKurt(idBand,:,:)),'rawKurt',channelName);
        
        exportHisto(squeeze(normDivDE(idBand,:,:)),'normDivDE',channelName)
        exportHisto(squeeze(normSignalDE(idBand,:,:)),'normSignalDE',channelName)
        exportHisto(squeeze(rawDE(idBand,:,:)),'rawDE',channelName);

        exportHisto(squeeze(numPeaks(idBand,:,:)),'numPeaks',channelName)

        exportHisto2D(squeeze(freqPairPeaks(idBand,:,:,:)),'frequencyPairPeaks',channelName,...
            signalToTestMat{idBand}.plotStart:0.5:signalToTestMat{idBand}.plotFinish,0:0.5:min(signalToTestMat{idBand}.plotFinish,32),freqBinsHighRes)

        % BOXPLOT
        exportBoxplotOutliers(squeeze(normDivMax(idBand,:,:)),'normDivMax_no_outliers',channelName);
        exportBoxplotOutliers(squeeze(normSignalMax(idBand,:,:)),'normSignalMax_no_outliers',channelName);
        exportBoxplotOutliers(squeeze(rawMax(idBand,:,:)),'rawMax_no_outliers',channelName);

        exportBoxplotOutliers(squeeze(normDivMean(idBand,:,:)),'normDivMean_no_outliers',channelName);
        exportBoxplotOutliers(squeeze(normSignalMean(idBand,:,:)),'normSignalMax_no_outliers',channelName);
        exportBoxplotOutliers(squeeze(rawMean(idBand,:,:)),'normMeanMax_no_outliers',channelName);

        exportBoxplotOutliers(squeeze(normDivStd(idBand,:,:)),'normDivStd_no_outliers',channelName);
        exportBoxplotOutliers(squeeze(normSignalStd(idBand,:,:)),'normSignalStd_no_outliers',channelName);
        exportBoxplotOutliers(squeeze(rawStd(idBand,:,:)),'rawStd_no_outliers',channelName);

        exportBoxplotOutliers(squeeze(normDivSkew(idBand,:,:)),'normDivSkew_no_outliers',channelName)
        exportBoxplotOutliers(squeeze(normSignalSkew(idBand,:,:)),'normSignalSkew_no_outliers',channelName)
        exportBoxplotOutliers(squeeze(rawSkew(idBand,:,:)),'rawSkew_no_outliers',channelName);

        exportBoxplotOutliers(squeeze(normDivKurt(idBand,:,:)),'normDivKurt_no_outliers',channelName)
        exportBoxplotOutliers(squeeze(normSignalKurt(idBand,:,:)),'normSignalKurt_no_outliers',channelName)
        exportBoxplotOutliers(squeeze(rawKurt(idBand,:,:)),'rawKurt_no_outliers',channelName)
        
        exportBoxplotOutliers(squeeze(normDivDE(idBand,:,:)),'normDivDE_no_outliers',channelName)
        exportBoxplotOutliers(squeeze(normSignalDE(idBand,:,:)),'normSignalDE_no_outliers',channelName)
        exportBoxplotOutliers(squeeze(rawDE(idBand,:,:)),'rawDE_no_outliers',channelName)

        exportBoxplotOutliers(squeeze(numPeaks(idBand,:,:)),'numPeaks_no_outliers',channelName)  
        
        % BOXPLOT WITH OUTLIERS
        exportBoxplot(squeeze(normDivMax(idBand,:,:)),'normDivMax',channelName);
        exportBoxplot(squeeze(normSignalMax(idBand,:,:)),'normSignalMax',channelName);
        exportBoxplot(squeeze(rawMax(idBand,:,:)),'rawMax',channelName);

        exportBoxplot(squeeze(normDivMean(idBand,:,:)),'normDivMean',channelName);
        exportBoxplot(squeeze(normSignalMean(idBand,:,:)),'normSignalMax',channelName);
        exportBoxplot(squeeze(rawMean(idBand,:,:)),'normMeanMax',channelName);

        exportBoxplot(squeeze(normDivStd(idBand,:,:)),'normDivStd',channelName);
        exportBoxplot(squeeze(normSignalStd(idBand,:,:)),'normSignalStd',channelName);
        exportBoxplot(squeeze(rawStd(idBand,:,:)),'rawStd',channelName);

        exportBoxplot(squeeze(normDivSkew(idBand,:,:)),'normDivSkew',channelName)
        exportBoxplot(squeeze(normSignalSkew(idBand,:,:)),'normSignalSkew',channelName)
        exportBoxplot(squeeze(rawSkew(idBand,:,:)),'rawSkew',channelName);

        exportBoxplot(squeeze(normDivKurt(idBand,:,:)),'normDivKurt',channelName)
        exportBoxplot(squeeze(normSignalKurt(idBand,:,:)),'normSignalKurt',channelName)
        exportBoxplot(squeeze(rawKurt(idBand,:,:)),'rawKurt',channelName)
        
        exportBoxplot(squeeze(normDivDE(idBand,:,:)),'normDivDE',channelName)
        exportBoxplot(squeeze(normSignalDE(idBand,:,:)),'normSignalDE',channelName)
        exportBoxplot(squeeze(rawDE(idBand,:,:)),'rawDE',channelName)

        exportBoxplotOutliers(squeeze(numPeaks(idBand,:,:)),'numPeaks',channelName) 
    end
end

%% qpcotr bispectrum calculation
%clc;
fprintf('Qpctor bispectrum calculation started ...\n')

% Available options to pass: 
% 'fullsignal','gamma','beta','alpha','theta','delta' 
signalToTest = fullsignal;

samples = signalToTest.samples;
M = fix(numel(samples)/8);
sp = (samples-mean(samples))/std(samples);
maxlag = fix(M/10);
ar_order = 29;
nfft = 512;
overlap = 0;
flag = 'unbiased'; %or 'biased'

tic
[ar_vec,bspec] = qpctor(sp,maxlag,ar_order,nfft,M,overlap,flag);
toc

%% Comparison of bispectrum between full and baseline signal
%clc;
fprintf('Bispectrum Direct started ...\n')

% Available options to pass: 
% 'fullsignal','gamma','beta','alpha','theta','delta' 
signalToTest = fullsignal;

samples = signalToTest.samples;
samples_b = signalToTest.baseline;
M = fix(numel(samples)/8); 
%M = 128;
M_b = fix(numel(samples_b)/3);
nfft = 2^nextpow2(M);
nfft_b = 2^nextpow2(M_b);
freqBins = 64/(nfft/2 - 1);
freqBins_b = 64/(nfft_b/2 - 1);
overlap = 50;
display = 0;

tic
[Bspec, waxis,zeroPos] = bispecd(samples,nfft,0,M,rate,overlap,display);
%[Bspec,waxis,zeroPos] = bispeci (samples,128,M, overlap,'unbiased', nfft, 1, display);
[Bspec_b, waxis_b,zeroPos_b] = bispecd(samples_b,nfft_b,0,M_b,rate,overlap,display);
toc

sdf = [0 8 0 8;0 16 8 16;8 16 8 16;0 32 32 64];
%choose region of bispectrum
q=3;

f1y = sdf(q,1);
f2y = sdf(q,2);
f1x = sdf(q,3);
f2x = sdf(q,4);

y=zeroPos+round(f1y/freqBins):zeroPos+round(f2y/freqBins);
y_b = zeroPos_b+round(f1y/freqBins_b):zeroPos_b+round(f2y/freqBins_b);

x=zeroPos+round(f1x/freqBins):zeroPos+round(f2x/freqBins);
x_b = zeroPos_b+round(f1x/freqBins_b):zeroPos_b+round(f2x/freqBins_b);

contourf(waxis(x),waxis(y),abs(Bspec(y,x)))
title('full signal')
figure
contourf(waxis_b(x_b),waxis_b(y_b),abs(Bspec_b(y_b,x_b)))
title('baseline signal')

%% peak dif
%clc;

numParticipants = 1;
numChannels = 2;  % max 32
numVideo = 40;     % max 40
numLabels = 4;
% FP1 FP2 AF3 AF4  F7 T7 P7 O1 0z O2 P8 CP6 T8 F8 Cz
%channelRange = [1 2 4 8 12 14 15 17 18 21 24 26 27 30 32]
channelRange = [1 2];
stringLabels = ["HVHA" "LVHA" "HVLA" "LVLA"];


M = 1024; 
%M = 128;
M_b = 384/3;
nfft = 2^nextpow2(M);
nfft_b = 2^nextpow2(M_b);
freqBins = 64/(nfft/2 - 1);
freqBins_b = 64/(nfft_b/2 - 1);
overlap = 50;
display = 0;


freqBand = [0 8 4 8;0 13 8 13;0 32 13 32;0 32 32 64];

euclideanD = zeros(4,4,numChannels*numVideo,numChannels,numParticipants);




for idParticipant = 1:numParticipants
    fprintf('Participant %d\n',idParticipant)
    countPerLabel = zeros(numLabels,1);
    for idVideo = 1:numVideo 
       for idChannel = channelRange
        [fullsignal,bands] = loadAndDecomposeDEAP(deapPath,idParticipant,idVideo,idChannel);
        emotionLabel = fullsignal.label;
        
         signalToTest = fullsignal;
         samples = signalToTest.samples;
         baseline = signalToTest.baseline;
         
         if emotionLabel == "HVHA"
                label = 1;
            elseif emotionLabel == "LVHA"
                label = 2;
            elseif emotionLabel == "HVLA"
                label = 3;
            elseif emotionLabel == "LVLA"
                label = 4;
         end
         countPerLabel(label) = countPerLabel(label) + 1;
         
         [Bspec, waxis,zeroPos] = bispecd(samples,nfft,0,M,rate,overlap,display);
         [Bspec_b, waxis_b,zeroPos_b] = bispecd(baseline,nfft_b,0,M_b,rate,overlap,display);
         clear triaBisp triaBisp_b partS partB 
         triaBisp=triu(flip(tril(flip(abs(Bspec(zeroPos:end,zeroPos:end))))));
         triaBisp_b=triu(flip(tril(flip(abs(Bspec_b(zeroPos_b:end,zeroPos_b:end))))));
          
         for ind =1:4
             f1y = round(freqBand(ind,1)/freqBins); f1y_b = round(freqBand(ind,1)/freqBins_b);
             f2y = round(freqBand(ind,2)/freqBins); f2y_b = round(freqBand(ind,2)/freqBins_b);
             f1x = round(freqBand(ind,3)/freqBins); f1x_b = round(freqBand(ind,3)/freqBins_b);
             f2x = round(freqBand(ind,4)/freqBins); f2x_b = round(freqBand(ind,4)/freqBins_b);
             
             partS = triaBisp(1+f1y:1+f2y,1+f1x:1+f2x);
             partB = triaBisp_b(1+f1y_b:1+f2y_b,1+f1x_b:1+f2x_b);
             smax = max(partS(:));
             [srow,scol] = find(triaBisp(1+f1y:1+f2y,1+f1x:1+f2x) == smax);
             bmax = max(partB(:));
             [brow,bcol] = find(triaBisp_b(1+f1y_b:1+f2y_b,1+f1x_b:1+f2x_b) == bmax);
             srow = freqBand(ind,1)+round(srow*freqBins);
             scol = freqBand(ind,3)+round(scol*freqBins);
             brow = freqBand(ind,1)+round(brow*freqBins_b);
             bcol = freqBand(ind,3)+round(bcol*freqBins_b);
             euclideanD(label,ind,countPerLabel(label),idChannel,idParticipant) = sqrt((srow-brow)^2 + (scol-bcol)^2);

             
         end
       end
    end
end
for idParticipant = 1:numParticipants
    for label = 1:4
        for idChannel = channelRange
        meandif.Theta(label,idChannel,idParticipant) = mean(nonzeros(euclideanD(label,1,:,idChannel,idParticipant)));
        stddif.Theta(label,idChannel,idParticipant) = std(nonzeros(euclideanD(label,1,:,idChannel,idParticipant)));
        
        meandif.Alpha(label,idChannel,idParticipant) = mean(nonzeros(euclideanD(label,2,:,idChannel,idParticipant)));
        stddif.Alpha(label,idChannel,idParticipant) = std(nonzeros(euclideanD(label,2,:,idChannel,idParticipant)));
        
        meandif.Beta(label,idChannel,idParticipant) = mean(nonzeros(euclideanD(label,3,:,idChannel,idParticipant)));
        stddif.Beta(label,idChannel,idParticipant) = std(nonzeros(euclideanD(label,3,:,idChannel,idParticipant)));
        
        meandif.Gamma(label,idChannel,idParticipant) = mean(nonzeros(euclideanD(label,4,:,idChannel,idParticipant)));
        stddif.Gamma(label,idChannel,idParticipant) = std(nonzeros(euclideanD(label,4,:,idChannel,idParticipant)));
        end
    end
end