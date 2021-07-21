function exportBispAverage(avgBisp,waxis,countPerLabel,yPlotRange,xPlotRange,zeroPos,fileName,channelName)
    f = figure;
    set(gcf, 'Position', get(0, 'Screensize'),'Visible','off'); % maximize figure

    stringLabels = ["HVHA" "LVHA" "HVLA" "LVLA"];
    countLabel = 1;
    for i = 1:4
        fprintf('Emotion %s: %d samples\n',stringLabels(i),countPerLabel(i))
        bispAvgPerLabel = squeeze(mean(squeeze(avgBisp(i,1:countPerLabel(i),:,:)),1));
        
        subplot(2,2,countLabel)
        plotBispBandFilled = bispAvgPerLabel(yPlotRange,xPlotRange);
        step = max(plotBispBandFilled(:))/10;
        contourf(waxis(zeroPos-1+xPlotRange),waxis(zeroPos-1+yPlotRange),plotBispBandFilled,0:step:max(plotBispBandFilled(:)))
        colorbar
        %hline = refline(1,0); hline.Color = 'r';
        title(stringLabels(i))
        
        countLabel = countLabel+1;
    end
    fullpath = sprintf('plots/average/%s/%s.png',channelName,fileName)
    exportgraphics(gcf,fullpath,'ContentType','image')
end