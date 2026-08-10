function [] = displayhilberttransformparams(significantpixelfeatures, pixels, displayall)
    arguments
        significantpixelfeatures
        pixels
        displayall = false
    end
    
    if displayall == true
        topx = length(significantpixelfeatures(:,1));
        usedpixels = significantpixelfeatures;
    elseif ~isempty(pixels)
        topx = length(pixels);
        usedpixels = zeros(topx, size(significantpixelfeatures, 2));
        for i=1:length(pixels)
            idx = find(significantpixelfeatures(:,1) == pixels(i), 1);
            if ~isempty(idx)
                usedpixels(i, :) = significantpixelfeatures(idx, :);
            end
        end
    else 
        topx = 0;
        usedpixels = [];
    end

    figure;
    hold on;
    legend_labels = cell(1, topx);

    for k = 1:topx
        tempdata = usedpixels(k, 2:end);
        timeaxis = 1:length(tempdata);
        plot(timeaxis, tempdata);
        legend_labels{k} = sprintf('pixel %d', usedpixels(k, 1));
    end
    
    hold off;    
    legend(legend_labels);
    xlabel('frame (# of frames)');
    ylabel('instantaneous freq');
    title('freq over time');
end