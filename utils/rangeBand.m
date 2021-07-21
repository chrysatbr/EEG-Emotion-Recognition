function [yRange,xRange,yPlotRange,xPlotRange,startIndex] = rangeBand(bispFilt,signalToTest,freqBins)
    [x,y] = size(bispFilt);
    startIndex  = max(round(signalToTest.start/freqBins),1);
    finishIndex = min(round(signalToTest.finish/freqBins,y-1));
    xRange = startIndex:finishIndex;
    yRange = 1:round(32/freqBins);
    
    startPlotIndex  = max(round(signalToTest.plotStart/freqBins),1);
    finishPlotIndex = min(round(signalToTest.plotFinish/freqBins,y-1));
    xPlotRange = startPlotIndex:finishPlotIndex; 
    yPlotRange = 1:round(min(signalToTest.plotFinish,32)/freqBins);
end