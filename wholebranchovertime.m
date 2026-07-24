setup

%% BRANCH INTENSITY

ibranch = 52;
thisbranchData = branchData{ibranch};

branchintensity = zeros(size(thisbranchData, 1), N_ims+4);
%makes a zero matrix of size (rows, columns); here, we have the same number of rows as the number of pixels along this branch
%originally, thisbranchData = #pixels x 4 (for the 4 cols of branchData{i}), but doing size(..., 1) picks the first dim
%the other dim, N_ims+4, includes the intensity info for all N_ims frames and the 4 extra are the 4 info cols

branchintensity(:,1:4) = thisbranchData; %chooses the first four columns of thisbranchData and pastes it into the branchintensity

for i_t = 1:N_ims
    im_smooth = squeeze(im(:,:,i_t)); %squeezes im at the i_t'th frame into a 2D matrix
    branchintensity(:,i_t+4) = im_smooth(branchintensity(:,3)); %here, branchintensity(:,3) are the linear position indices for pixels 
    %so we're kinda going pixel by pixel from the 2D frame from the previous line
    %and accessing the value at each pixel (the value = intensity)
    %and we add it to branchintensity's i_t+4th column
end