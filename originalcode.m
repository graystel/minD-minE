bfr = BioformatsImage('E:\MinD-MinE\20260701\old\OLD_WAVES_GFP_80.btf');
i_tstart = 0;
N_ims = 240;
if ~exist('im','var') || ~exist('imframe','var')
   for i_t = 1:N_ims
       imframe = getPlane(bfr, 1, 1, i_tstart+i_t); % pull out a specific image from the btf file
       if size(imframe,3)==3 % if it is in color, then remove the color data
           imframe = squeeze(imframe(:,:,2)); % remove color data
       end
      
       if i_t == 1
           im = zeros(size(imframe,1),size(imframe,2),N_ims);
       end
  
       im(:,:,i_t) = imframe; % make a 3D stack of images, containing all frames
      
   end
  
   % segment the hyphae, focusing on the ones that have a lot of intensity
   % change (meaning the reaction diffusion system is 'working'
  
   % [Ny,Nx] = size(imframe);
end
imactivity = sum(abs(diff(double(im),1,3)),3);
imactivity = (imactivity-min(imactivity(:)))/(max(imactivity(:))-min(imactivity(:)));
hyphaeThicknessRange = [10, 40]; % range of hyphal diameters for the Hessian filter
% Hessian filter to enhance tubular objects
% 'StructureSensitivity' controls how aggressively it rejects non-tubular shapes
enhancedImg = fibermetric(imactivity, hyphaeThicknessRange, 'StructureSensitivity', 0.01);
% use intensity thresholding to mask out the hyphae
baseThreshold = multithresh(enhancedImg,1); % allow for multiple gray thresholds
firstMask = enhancedImg>baseThreshold(1);
minArea = 2000; % reject anything smaller than this
finalMask = bwareaopen(firstMask, minArea);
cleanMask = imfill(finalMask, 'holes');
% Remove tiny jagged edge artifacts that cause spurs.
% We use morphological closing (fills tiny indentations) followed by opening (removes tiny bumps).
radius = 2; % Adjust this if your hyphae are particularly thick or thin
se = strel('disk', radius);
cleanMask = imclose(cleanMask, se);
cleanMask = imopen(cleanMask, se);
% Skeletonization and Pruning
pruningThreshold = 25;
% Generate the clean skeleton
skeleton = bwskel(cleanMask, 'MinBranchLength', pruningThreshold);
[Ny,Nx] = size(skeleton);
% 5a. Isolate individual branches by removing branch points
branchPoints = bwmorph(skeleton, 'branchpoints');
isolatedBranches = skeleton & ~branchPoints;
% 5b. Find all distinct, disconnected skeletal segments
cc = bwconncomp(isolatedBranches, 8);
% Preallocate a cell array to store the parameterization for each branch.
% Each cell will contain an Nx3 matrix: [X-coord, Y-coord, ArcLength_s]
branchData = cell(cc.NumObjects, 1);
for i = 1:cc.NumObjects
   % Create a temporary binary mask for just this single branch
   singleBranch = false(size(skeleton));
   singleBranch(cc.PixelIdxList{i}) = true;
   % Find the endpoints of this specific branch
   endPts = bwmorph(singleBranch, 'endpoints');
   [epY, epX] = find(endPts);
   if isempty(epX)
       % Edge case: The branch is a perfectly closed loop (no endpoints).
       % We pick an arbitrary starting point.
       [startY, startX] = find(singleBranch, 1, 'first');
   else
       % Normal case: Pick the first endpoint as the origin (s = 0)
       startY = epY(1);
       startX = epX(1);
   end
   % Calculate the geodesic distance from the start point to all other
   % pixels *along* the branch path (using quasi-euclidean metric for accuracy).
   geoDist = bwdistgeodesic(singleBranch, startX, startY, 'quasi-euclidean');
   % Extract the coordinates for this branch
   [y, x] = find(singleBranch);
   % Extract the corresponding arc length distances using the linear indices
   s = geoDist(cc.PixelIdxList{i});
   % Sort the pixels by their arc length so the coordinates are ordered
   % sequentially from the start of the branch to the end.
   [s_sorted, sortIdx] = sort(s);
   x_sorted = x(sortIdx);
   y_sorted = y(sortIdx);
   ind_sorted = sub2ind([Ny,Nx],y_sorted,x_sorted);
   % Store the parameterized data
   branchData{i} = [x_sorted, y_sorted, ind_sorted, s_sorted]; % contains 4 columns: the first two coliumns are the coordinates
   % the next column is the indexed ordinates of those points, and then
   % the arc lengths along the hypha
end
%% code snippet for isolating intensity traces for a given hypha
% pick a specific branch to get intensity data from, run the code snippet
% at the end to see what each branch is
ibranch = 52;
thisbranchData = branchData{ibranch};
branchintensity = zeros(size(thisbranchData,1),N_ims+4);
% output: first four columns are the x,y,index and s data
branchintensity(:,1:4) = thisbranchData;
for i_t = 1:N_ims
   % im_smooth = imgaussfilt(squeeze(im(:,:,i_t)),4); % cleans up the image a bit
   im_smooth = squeeze(im(:,:,i_t)); % option to add smoothing at this stage
   branchintensity(:,i_t+4) = im_smooth(branchintensity(:,3));
end
%% code snippet for plotting intensity against s
figure(1);
N_frames = 5; dt = 2;
cmap = colormap(parula(N_frames));
for i_frame = 1:N_frames
   plot(branchintensity(:,4),branchintensity(:,5+(i_frame-1)*dt),'color',cmap(i_frame,:)); % plot intensity against s
   hold on;
end
xlabel('s (pixels)'); ylabel('intensity');
% Finalize the plot with a title and legend
title('Branch Intensity Profiles');
legend(arrayfun(@(x) sprintf('Frame %d', x), 1:dt:1+(N_frames-1)*dt, 'UniformOutput', false));
hold off;
print('/Users/mroper/Library/CloudStorage/Dropbox/syncytium/intensityonhypha.png');
%% code snippet for velocimetry - incomplete!!!
vel = zeros(Ns,iframe);
windowsize = 31; % how many pixels are put in the frame measurement (I am being
% a bit careless that s is not a uniform grid)
hwindow = (windowsize-1)/2;
for i_frame = 1:N_frames-1
   kym1 = branchintensity(:,5+i_frame);
   kym2 = branchintensity(:,5+i_frame+1);
   Ns = size(kym1,1);
   for i_s = hwindow+1:Ns-hwindow       
        [r,lags] = xcorr(kym1(i_s-hwindow:i_s+hwindow),kym2,hwindow);
   end
end
%% code snippet for making kymographs
figure(4);
colormap parula(256);
imagesc(branchintensity(:,5:end).');
xlabel('arc length (s)'); ylabel('frame (t)');
%% code snippet for labeling the hyphae
figure(2); clf;
imagesc(squeeze(im(:,:,1))); hold on;
for i = 1:length(branchData)
   % Extract the X and Y coordinates for this specific branch
   X = branchData{i}(:, 1);
   Y = branchData{i}(:, 2);
   plot(X,Y,'r');
   hold on;
   plot(X(1),Y(1),'ro');
   plot(X(end),Y(end),'o','color',[0 0.7 0]);
   % Safeguard: Ensure the branch isn't empty
   if ~isempty(X)
       % Find the middle index of the branch to anchor the text
       midIdx = max(1, round(length(X) / 2));
       labelX = X(midIdx);
       labelY = Y(midIdx);
       % Plot the text label
       % Using cyan/yellow/magenta gives high contrast against black & white
       text(labelX, labelY, num2str(i), ...
           'Color', 'white', ...
           'FontSize', 10, ...
           'FontWeight', 'bold', ...
           'HorizontalAlignment', 'center', ...
           'BackgroundColor', 'black', ... % Optional: adds a background box for readability
           'Margin', 1);
   end
end

