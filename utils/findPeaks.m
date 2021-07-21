% find peaks
function [bispPeaks,peakInfo] = findPeaks(bispPeaks,xOffset,freqBins)
    [row,col] = size(bispPeaks);
    threshold = 0.8 * max(bispPeaks(:));
    for i = 1:row
        for j = 1:col
            if bispPeaks(i,j) <= threshold
                bispPeaks(i,j) = 0;
            end
        end
    end
    
    % max
    [yMax,xMax,value] = maxMatrix(bispPeaks);
    peakInfo(1).f2 = (yMax-1)*freqBins;
    peakInfo(1).f1 = (xMax-1)*freqBins + (xOffset-1)*freqBins;
    peakInfo(1).value = value;
            
    % clear nearby peaks widthHertz x widthHertz Hz patch
    widthHertz = 1;
    patchSize = fix(widthHertz/(2*freqBins));
    pks = unique(sort(nonzeros(bispPeaks),'descend'),'stable');
    for i = 1:numel(pks)
        if(numel(pks) == 1) break; end
        [x,y] = find(bispPeaks == pks(i));
        N = numel(x);
        for j = 1:N
            if N==1 posX = x; else posX = x(j); end
            if N==1 posY = y; else posY = y(j); end
            
            % avoid clearing out already cleared data
            if bispPeaks(posX,posY) == 0 continue; end
           
            % before
            %fprintf('before\n')
            %bispPeaks(max(posX-patchSize,1):min(posX+patchSize,row),max(posY-patchSize,1):min(posY+patchSize,col))
            
            bispPeaks(max(posX+1,1):min(posX+patchSize,row),max(posY-patchSize,1):min(posY+patchSize,col)) = 0;
            bispPeaks(max(posX-patchSize,1):min(posX-1,row),max(posY-patchSize,1):min(posY+patchSize,col)) = 0;
            bispPeaks(posX,max(posY-patchSize,1):max(posY-1,1)) = 0;
            bispPeaks(posX,min(posY+1,col):min(posY+patchSize,col)) = 0;
            
            % after
            %fprintf('after\n')
            %bispPeaks(max(posX-patchSize,1):min(posX+patchSize,row),max(posY-patchSize,1):min(posY+patchSize,col))
        end
    end
    
    count = 1;
    for i = 1:row
        for j = 1:col
            if bispPeaks(i,j) ~= 0
                peakInfo(count).f2    = (i-1)*freqBins;
                peakInfo(count).f1    = (j-1)*freqBins + (xOffset-1)*freqBins;
                peakInfo(count).value = bispPeaks(i,j); 
                count = count + 1;
            end
        end
    end
end