%% IF NOTHING WORKS

%clear;
%clc;
%close all;

%% SET UP

oldwaves = BioformatsImage('E:\MinD-MinE\20260701\old\OLD_WAVES_GFP_80.btf'); %importing the image of interest
i_tstart = 0; %we expect to start from frame 0
N_ims = 240; %we expect to encounter/go through all the frames, i.e. all 240 frames

dt = 2;
N_frames = 5;

if ~exist('im', 'var') || ~exist('imframe', 'var') %i.e. if a variable named im or imframe don't exist, then we want to prepare them
    for i_t = 1:N_ims %for loop through all the frames
        imframe = getPlane(oldwaves, 1, 1, i_tstart+i_t); %gets current frame from the btf file
        if size(imframe, 3) == 3 %if it's in color, then we wanna remove this color data 
            imframe = squeeze(imframe(:,:,2)); %so we remove it here
        end

        if i_t == 1 %if we just started, we prepare a data structure to receive all the processed frames in later iterations
            im = zeros(size(imframe,1), size(imframe, 2), N_ims); %makes an array of zeros of size matching each frame, for all 240 frames
        end

        im(:, :, i_t) = imframe; %storing the frame into the array
    end

    %[Ny, Nx] = size(imframe); %we're going to keep track of the size of the frame
end

imactivity = sum(abs(diff(double(im), 1, 3)), 3); 
%converts the values in im to doubles for calculation
%then takes diff (1st difference, meaning btwn consecutive frames?) of each value in the array with regards to time (3rd dim, as def'd for im)
%since any direction of difference represents "activity" in some regard and thus we want them to sum together w/o regard to direction, we take abs
%sum sums all the consecutive differences between each frame for each value in the array
%as a result, imactivity is a 2D array
imactivity = (imactivity-min(imactivity(:)))/(max(imactivity(:))-min(imactivity(:))); %looks like normalization?

hyphaeThicknessRange = [10, 40]; %we anticipate that the hyphae we're interested in are going to have this range of diammeters 
enhancedImg = fibermetric(imactivity, hyphaeThicknessRange, 'StructureSensitivity', 0.01); 
%fibermetric is a Bioformats toolbox function that looks for tubular or line-like structures
%'StructureSensitivity' controls the aggressiveness for rejecting non-tubular shapes
%the resulting enhancedImg is a 2D array (modified imactivity) with the modifications being enhancement of the tubes

baseThreshold = multithresh(enhancedImg, 1); %multithresh is a Matlab func that tries to separate groups of pixels based on the 
%second argument's value (i.e. here, we want to find 1 threshold value) in enhancedImg
firstMask = enhancedImg > baseThreshold(1); %we're using the baseThreshold's first (only, in this case) threshold value to give a 
%matrix of binary outputs for each pixel and seeing if their threshold value is sufficiently high enough to represent actual activity / oscillation

minArea = 2000; %since we don't want hyphae that are too short, as they are uninteresting, we set a minimum area of the branch
finalMask = bwareaopen(firstMask, minArea); %removes any areas in the firstMask of the image that are too small

cleanMask = imfill(finalMask, "holes"); %imfill fills enclosed background regions (holes = background/0 regions that cannot reach the edge of the image)

radius = 2; %higher number = more smoothing along the length of the hyphae 
se = strel('disk', radius); %this structuring element is like a filter in convolutions
cleanMask = imclose(cleanMask, se); %this is the dilation into erosion which fills small gaps + connects nearby objects
cleanMask = imopen(cleanMask, se); %this is the erosion into dilation which removes small isolated objects and smooths edges

%result is cleaner


%% SKELETONIZATION

pruningThreshold = 25;
skeleton = bwskel(cleanMask, 'MinBranchLength', pruningThreshold); 
%bwskeleton makes a skeleton (1 pixel wide) of the filled/1 values of cleanMask with a minimum branch length of 25
[Ny, Nx] = size(skeleton);
branchPoints = bwmorph(skeleton, "branchpoints"); %asks to find the branch points (bwmorph performs morphological operations on binary images)
isolatedBranches = skeleton & ~branchPoints; %keeps the skeletons without the branchpoints

%% BRANCHDATA

cc = bwconncomp(isolatedBranches, 8); %bwconncomp returns groups of connected components (i.e. the branches)
branchData = cell(cc.NumObjects, 1); %creates an empty data structure (cells of length equal to the number of objects in cc, width 1)

for i = 1:cc.NumObjects
    singleBranch = false(size(skeleton)); %makes a binary image that's all blank/0/false thats the same size as the skeleton 
    singleBranch(cc.PixelIdxList{i}) = true; %fills the pixels of singlebranch with the current branch (i.e. just getting the branches iteratively)
    %PixelIdxList{i} is a cell array so curly brackets, and it returns the pixel list of the ith object in cc
    %setting them to true just fills in those pixel locations

    endPts = bwmorph(singleBranch, 'endpoints'); %finds the endpoints of this specific branch 
    [epY, epX] = find(endPts);

    if isempty(epX) %if there's a loop, manually sets the first pixel as a starting point
        [startX, startY]  = find(singleBranch, 1, 'first');
    else
        startX = epX(1); %column values
        startY = epY(1); %row values
    end

    geoDist = bwdistgeodesic(singleBranch, startX, startY, 'quasi-euclidean'); 
    %calc the distance from the start pt to all other pixels in the same branch along the skeleton...quasi euclidean matric for accuracy

    [y, x] = find(singleBranch); %gets the coordinates for this branch (i.e. every pixel along this branch)
    s = geoDist(cc.PixelIdxList{i}); %extract the arc length distances via linear indices (indices as if an image were just a long 1-dim line of pixels)
    %as a result, s is a list of distances for each pixel along a branch from the start point of the branch

    [s_sorted, sortIdx] = sort(s); %sort the pixels by their arc lengths s so the coordinates are ordered sequentially
    %sortIdx retains info about how s was sorted
    x_sorted = x(sortIdx); %and now we sort x values the same way 
    y_sorted = y(sortIdx); %and y values the same way
    ind_sorted = sub2ind([Ny, Nx], y_sorted, x_sorted); %creates a linear pixel index (i.e. effectively for each pixel on this branch, give it a position)

    %store this data
    branchData{i} = [x_sorted, y_sorted, ind_sorted, s_sorted]; %first two columns are the coords
    %the next is the "positions"
    %the last column are the arc lengths along the hypha
end
