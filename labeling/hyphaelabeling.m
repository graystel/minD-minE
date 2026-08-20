function [] = hyphaelabeling(im, branchData)
    figure;
    %clf; %clear the figure of current window
    
    imagesc(squeeze(im(:,:,1))); %squeeze takes the 1st frame from the im array; imagesc displays an image and autoscales color brightness
    hold on; %keeps the current graphic active so later plots don't replace this graphic

    for i = 1:length(branchData)
        X = branchData{i}(:, 1); %gets all the x-coordinates (the first column) of the ith branch
        Y = branchData{i}(:, 2); %gets all the y-coordinates (the second column) of the ith branch
        
        plot(X,Y,'r'); %plots the branch (skeleton) and draws a red line through it
        hold on; %keep the graphic 

        plot(X(1),Y(1),'ro'); %plots the first point of the branch as a red dot
        plot(X(end),Y(end),'o','color',[0 0.7 0]); %marks the end point of the branch as a green circle
        
        if ~isempty(X) %as long as the branch isn't empty, then....
            midIdx = max(1, round(length(X) / 2)); %find the middle index of the branch to anchor the text
            labelX = X(midIdx); 
            labelY = Y(midIdx);
            text(labelX, labelY, num2str(i), ... %plot the text label
                'Color', 'white', ...
                'FontSize', 10, ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', ...
                'BackgroundColor', 'black', ...
                'Margin', 1);
        end
    end
end