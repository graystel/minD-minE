addpath('../minD-minE/setup');
addpath('../minD-minE/preliminary')
addpath('../minD-minE/temporalanalysis')

% original file path 
% oldwaves = BioformatsImage('E:\MinD-MinE\20260701\old\OLD_WAVES_GFP_80.btf');

i_tstart = 0; 
N_ims = 240; %we expect to encounter/go through all the frames, i.e. all 240 frames
N_frames = 5;
dt = 2;
reset = true;
hyphaeThicknessRange = [10, 40];
minArea = 2000;

[im, branchData] = setup('/Users/heatherli491/Downloads/OLD_WAVES_GFP_80.btf', N_ims, i_tstart, N_frames, dt, reset, hyphaeThicknessRange, minArea);

fprintf("checkpt 1");

%hyphaelabeling(im, branchData);

ibranch = 56;
topx = 10;
display = true;
shiftin = 100;
step = 10;
bwindow = 11;
a = 1;

[sorted_scores, sorted_pixel_indices] = singlepixelswithbestperiodicity(@singlepixelovertime, @singlepixelallfilters, @singlepixelxcorrshift, branchData, im, N_ims, ibranch, step, topx, display, shiftin, bwindow, a);

fprintf("checkpt 2");

significantpixelfeatures = gethilberttransformparams(@hilberttransform, sorted_pixel_indices, branchData, im, N_ims, ibranch);

fprintf("checkpt 3");

displayhilberttransformparams(significantpixelfeatures, sorted_pixel_indices(1:10));

%singlepixelhilbertreconstruction(params, true);

