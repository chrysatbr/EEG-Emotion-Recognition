% filtering the main quadrant of bispectrum
function [bispFilt,peakInfo] = bispecdFilter(magnBisp,freqBins)
    [row,col] = size(magnBisp);
    filter = triu(ones(row,col));
    bispFilt = filter .* magnBisp;
    threshold = 0.7 * max(bispFilt(:));
    for i = 1:row
        for j = 1:col
            if bispFilt(i,j) <= threshold
                bispFilt(i,j) = 0;
            end
        end
    end
        
    % clear nearby peaks 1x1 Hz patch
    widthHertz = 1;
    patchSize = fix(widthHertz/(2*freqBins));
    pks = unique(sort(nonzeros(bispFilt),'descend'));
    for i = 1:numel(pks)
        [x,y] = find(bispFilt == pks(i));
        N = numel(x);
        for j = 1:N
            if N==1 posX = x; else posX = x(j); end
            if N==1 posY = y; else posY = y(j); end
            
            % avoid clearing out already cleared data
            if bispFilt(posX,posY) == 0 continue; end
            
            % before
            %fprintf('before\n')
            %bispFilt(posX-patchSize:posX+patchSize,posY-patchSize:posY+patchSize)
            
            bispFilt(posX+1:posX+patchSize,posY-patchSize:posY+patchSize) = 0;
            bispFilt(posX-patchSize:posX-1,posY-patchSize:posY+patchSize) = 0;
            bispFilt(posX,posY-patchSize:posY-1) = 0;
            bispFilt(posX,posY+1:posY+patchSize) = 0;
            
            % after
            %fprintf('after\n')
            %bispFilt(posX-patchSize:posX+patchSize,posY-patchSize:posY+patchSize)
        end
    end
    
    count = 1;
    for i = 1:row
        for j = 1:col
            if bispFilt(i,j) ~= 0
                peakInfo(count).f2    = (i-1)*freqBins;
                peakInfo(count).f1    = (j-1)*freqBins;
                peakInfo(count).value = bispFilt(i,j); 
                count = count + 1;
            end
        end
    end
end