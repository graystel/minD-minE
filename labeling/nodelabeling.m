function [] = nodelabeling(im, pixeltoVertex)
    figure;
    %clf; %clear the figure of current window
    
    imagesc(squeeze(im(:,:,1))); %squeeze takes the 1st frame from the im array; imagesc displays an image and autoscales color brightness
    hold on; %keeps the current graphic active so later plots don't replace this graphic

    for i = 1:length(pixeltoVertex)
        [Y, X] = find(pixeltoVertex == i);
        
        plot(X,Y,'o','MarkerSize', 5, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); %colors all the x, y pixels red?
        hold on; %keep the graphic 
        
        if ~isempty(X) %as long as the vertex isn't empty, then....
            midIdx = max(1, round(length(X) / 2)); %find the middle index of the vertex to anchor the text
            labelX = X(midIdx); 
            labelY = Y(midIdx);
            text(labelX+3, labelY-3, num2str(i), ... %plot the text label
                'Color', 'black', ...
                'FontSize', 10, ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', ...
                'BackgroundColor', 'white', ...
                'Margin', 1);
        end
    end
end