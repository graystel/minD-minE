function [] = displayhilberttransformparams(significantpixelfeatures, displayall, pixels)
    if nargin < 3 || isempty(pixels)
        pixels = [];
    end
    if nargin < 2 || isempty(displayall)
        displayall = false;
    end
    
    if displayall == true
        topx = length(significantpixelfeatures(:,1));
        usedpixels = significantpixelfeatures;
    elseif ~isempty(pixels)
        topx = length(pixels);
        usedpixels = cell(length(pixels), length(significantpixelfeatures(:,1)));
        for i=1:length(pixels)
            idx = find(significantpixelfeatures(:,1) == pixels(i), 1);
            for j=1:length(significantpixelfeatures(1,:))
                usedpixels(i, j) = significantpixelfeatures(idx, j+1);
            end
        end
    else 
        topx = 0;
        usedpixels = [];
    end

    figure;
    legend_labels = cell(1, topx);

    for k = 1:topx
        hold on;
        tempdata = usedpixels(k, 2:end);
        plot(1:(N_ims-1), tempdata);
        legend_labels{k} = sprintf('pixel %d', usedpixels(k, 1));
    end    
        
    legend(legend_labels);
    xlabel('frame (# of frames)');
    ylabel('instantaneous freq');
    title('freq over time');
end