addpath('../minD-minE/setup');
addpath('../minD-minE/labeling')
addpath('../minD-minE/temporalanalysis')
addpath('../minD-minE/tools')

% original file path 
% oldwaves = BioformatsImage('E:\MinD-MinE\20260701\old\OLD_WAVES_GFP_80.btf');


% %% VERIFIED TEMPORALANALYSIS FILES
% i_tstart = 0; 
% N_ims = 240; %we expect to encounter/go through all the frames, i.e. all 240 frames
% N_frames = 5;
% dt = 2;
% reset = true;
% hyphaeThicknessRange = [10, 40];
% minArea = 2000;
% 
% [im, branchData, ~, ~, ~] = setup('/Users/heatherli491/Downloads/OLD_WAVES_GFP_80.btf', N_ims, i_tstart, N_frames, dt, reset, hyphaeThicknessRange, minArea);
% 
% fprintf("checkpt 1");
% 
% %hyphaelabeling(im, branchData);
% 
% ibranch = 56;
% topx = 10;
% display = true;
% shiftin = 100;
% step = 10;
% bwindow = 11;
% a = 1;
% 
% [sorted_scores, sorted_pixel_indices] = singlepixelswithbestperiodicity(@singlepixelovertime, @singlepixelallfilters, @singlepixelxcorrshift, branchData, im, N_ims, ibranch, step, topx, display, shiftin, bwindow, a);
% 
% fprintf("checkpt 2");
% 
% significantpixelfeatures = gethilberttransformparams(@hilberttransform, sorted_pixel_indices, branchData, im, N_ims, ibranch);
% 
% fprintf("checkpt 3");
% 
% displayhilberttransformparams(significantpixelfeatures, sorted_pixel_indices(1:10));

%% VERIFIED EDGES 

i_tstart = 0; 
N_ims = 240; %we expect to encounter/go through all the frames, i.e. all 240 || 145 || 50 frames
N_frames = 5;
dt = 2;
reset = true;
hyphaeThicknessRange = [10, 40]; %[10, 40] for original 
minArea = 2000;

[im, branchData, pixeltoVertex, adj, edges] = setup('/Users/heatherli491/Downloads/OLD_WAVES_GFP_80.btf',@loopremover, N_ims, i_tstart, N_frames, dt, reset, hyphaeThicknessRange, minArea, true);
%WAVES_MinD_GFP_1_hyphae1.tif
%OLD_WAVES_GFP_80.btf // the original
fprintf("checkpt 1\n");

hyphaenodelabeling(im, branchData, pixeltoVertex);
nodelabeling(im, pixeltoVertex);

%% DOESN'T WORK YET
% i_tstart = 0; 
% N_ims = 240; %we expect to encounter/go through all the frames, i.e. all 240 frames
% N_frames = 5;
% dt = 2;
% reset = true;
% hyphaeThicknessRange = [10, 40];
% minArea = 2500;
% 
% [im, branchData, pixeltoVertex, adjacencyMatrix, edgeMatrix] = setup('/Users/heatherli491/Downloads/OLD_WAVES_GFP_80.btf', @loopremover, N_ims, i_tstart, N_frames, dt, reset, hyphaeThicknessRange, minArea, true);
% 
% fprintf("checkpt 1");
% 
% hyphaenodelabeling(im, branchData, pixeltoVertex);
% 
% ibranch = 47;%56
% topx = 10;
% display = true;
% shiftin = 100;
% step = 10;
% bwindow = 11;
% a = 1;
% 
% [sorted_scores, sorted_pixel_indices] = singlepixelswithbestperiodicity(@singlepixelovertime, @singlepixelallfilters, @singlepixelxcorrshift, branchData, im, N_ims, ibranch, step, topx, display, shiftin, bwindow, a);
% 
% fprintf("checkpt 2");
% 
% for i=1:length(sorted_pixel_indices)
%     pixel_idx = sorted_pixel_indices(i);
%     px = branchData{ibranch}(pixel_idx, 1);
%     py = branchData{ibranch}(pixel_idx, 2);
%     currentpixel = squeeze(double(im(py, px, :)));
%     [~, ~, currentfiltered] = singlepixelallfilters(currentpixel, pixel_idx, ibranch, bwindow, a, false);
%     currentparams = hilberttransform(currentfiltered);
%     singlepixelhilbertreconstruction(currentfiltered, currentparams, pixel_idx, true);
% end


