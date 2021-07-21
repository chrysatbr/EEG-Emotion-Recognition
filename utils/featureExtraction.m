function featureExtraction(deapPath,idChannel)

numParticipants = 32;
numChannels     = 32;     % max 32
numVideo        = 40;     % max 40
numLabels       = 4;
rate            = 128;

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
        startPlotIndex  = max(round((gamma.start - gamma.plotStart)/freqBinsHighRes),1);
        finishPlotIndex = min(round((gamma.finish - gamma.plotFinish)/freqBinsHighRes),256);
        xPlotRange = startPlotIndex:finishPlotIndex; 
        yPlotRange = 1:round(min(gamma.plotFinish,32)/freqBinsHighRes);
        exportBispAverage(avgBispGamma,waxis,squeeze(countPerLabel(idBand,:)),yPlotRange,xPlotRange,zeroPos,'avgGammasignal',channelName)
    elseif idBand == 3
        startPlotIndex  = max(round((beta.start - beta.plotStart)/freqBinsHighRes),1);
        finishPlotIndex = 128;
        xPlotRange = startPlotIndex:finishPlotIndex; 
        yPlotRange = 1:round(min(beta.plotFinish,32)/freqBinsHighRes);
        exportBispAverage(avgBispBeta,waxis,squeeze(countPerLabel(idBand,:)),yPlotRange,xPlotRange,zeroPos,'avgBetasignal',channelName)
    elseif idBand == 4
        startPlotIndex  = max(round((alpha.start - alpha.plotStart)/freqBinsHighRes),1);
        finishPlotIndex = 64;
        xPlotRange = startPlotIndex:finishPlotIndex; 
        yPlotRange = 1:round(min(alpha.plotFinish,32)/freqBinsHighRes);
        exportBispAverage(avgBispAlpha,waxis,squeeze(countPerLabel(idBand,:)),yPlotRange,xPlotRange,zeroPos,'avgAlphasignal',channelName)
    elseif idBand == 5
        startPlotIndex  = max(round((theta.start - theta.plotStart)/freqBinsHighRes),1);
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
    exportBoxplotOutliers(squeeze(normSignalMean(idBand,:,:)),'normSignalMean_no_outliers',channelName);
    exportBoxplotOutliers(squeeze(rawMean(idBand,:,:)),'normMean_no_outliers',channelName);

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
    exportBoxplot(squeeze(normSignalMean(idBand,:,:)),'normSignalMean',channelName);
    exportBoxplot(squeeze(rawMean(idBand,:,:)),'normMean',channelName);

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
