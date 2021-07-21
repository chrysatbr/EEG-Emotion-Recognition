function exportHisto2D(freqPairPeaks,fileName,channelName,freqRange,yRange,freqBins)
    stringLabels = ["HVHA" "LVHA" "HVLA" "LVLA"];    
    for idLabel = 1:4
        f = figure;
        set(gcf, 'Position', get(0, 'Screensize'),'Visible','off'); % maximize figure

        % remove zeros due to dynamic allocation
        x = nonzeros(squeeze(freqPairPeaks(idLabel,1,:)));
        y = nonzeros(squeeze(freqPairPeaks(idLabel,2,:)));

        hist3([x,y],'CDataMode','auto','Ctrs',{freqRange yRange});
        h = hist3([x,y],'CDataMode','auto','Ctrs',{freqRange yRange});
        
        % get the most frequent element
        [row,col,value] = maxMatrix(h);
        f1 = (row-1)*freqBins;
        f2 = (col-1)*freqBins;
        
        %text = sprintf("%s\nMost frequent: f1 = %0.2f, f2 = 0.2%f, value = %d",...
            %stringLabels(idLabel),f1,f2,value);
        %title(text)
        xticks(freqRange);yticks(yRange)
        xlabel('f1');ylabel('f2')
        colorbar
        view(2)

        fullpath = sprintf('plots/histo/%s/%s_%s.png',channelName,fileName,stringLabels(idLabel))
        exportgraphics(gcf,fullpath,'ContentType','image')
    end
end