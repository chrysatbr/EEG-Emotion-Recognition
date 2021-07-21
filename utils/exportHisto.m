function exportHisto(data,fileName,channelName)
    stringLabels = ["HVHA" "LVHA" "HVLA" "LVLA"];
    f = figure;
    set(gcf, 'Position', get(0, 'Screensize'),'Visible','off'); % maximize figure
    count  = 1;
    for idLabel = 1:4
        subplot(2,2,count)
        dataPerLabel = nonzeros(data(idLabel,:));
        h = histogram(rmoutliers(dataPerLabel),'normalization','probability');
        morebins(h); morebins(h);
        area = numel(rmoutliers(dataPerLabel))/numel(dataPerLabel);
        text = sprintf('%s\nMean:%0.5f, Mean (no outliers):%0.5f, Area: %0.2f',...
            stringLabels(idLabel),mean(dataPerLabel),mean(rmoutliers(dataPerLabel)),area);
        title(text)
        count = count + 1;
    end
    %sgtitle('Histogram Standard Deviation Bispectrum')
    fullpath = sprintf('plots/histo/%s/%s.png',channelName,fileName)
    exportgraphics(gcf,fullpath,'ContentType','image')
end