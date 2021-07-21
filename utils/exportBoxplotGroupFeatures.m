function exportBoxplotGroupFeatures(data,featureNames,fileName,channelName)
    stringLabels = ["HVHA" "LVHA" "HVLA" "LVLA"];
    f = figure;
    set(gcf, 'Position', get(0, 'Screensize'),'Visible','off'); % maximize figure
    count  = 1;
    numFeatures = numel(data);
    for idLabel = 1:4
        subplot(2,2,count)
        
        dataClean = data{1};
        dataClean = nonzeros(dataClean(idLabel,:));
        dataPerLabel = zeros(numFeatures,numel(dataClean));
        
        for i = 1:numel(data)
            feature = data{i};
            feature = nonzeros(feature(idLabel,:));
            dataPerLabel(i,:) = feature;
        end
                
        boxplot(dataPerLabel',featureNames,'PlotStyle','compact','LabelOrientation','horizontal');
        title(stringLabels(idLabel))
        count = count + 1;
    end
    fullpath = sprintf('plots/boxplot/%s/%s.png',channelName,fileName)
    exportgraphics(gcf,fullpath,'ContentType','image')
end