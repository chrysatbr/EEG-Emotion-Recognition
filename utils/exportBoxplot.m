function exportBoxplot(data,fileName,channelName)
    f = figure;
    set(gcf, 'Position', get(0, 'Screensize'),'Visible','off'); % maximize figure
    
    vectorBoxplot = [nonzeros(data(1,:));nonzeros(data(2,:));nonzeros(data(3,:));nonzeros(data(4,:))];
    g1 = repmat({'HVHA'},numel(nonzeros(data(1,:))),1);
    g2 = repmat({'LVHA'},numel(nonzeros(data(2,:))),1);
    g3 = repmat({'HVLA'},numel(nonzeros(data(3,:))),1);
    g4 = repmat({'LVLA'},numel(nonzeros(data(4,:))),1);
    g = [g1; g2; g3; g4];
    boxplot(vectorBoxplot,g,'PlotStyle','compact','LabelOrientation','horizontal')
    
    fullpath = sprintf('plots/boxplot/%s/%s.png',channelName,fileName)
    exportgraphics(gcf,fullpath,'ContentType','image')
end