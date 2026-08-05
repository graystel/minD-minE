addpath('../minD-minE/setup');
addpath('../minD-minE/preliminary')
addpath('../minD-minE/temporalanalysis')

% original file path 
% oldwaves = BioformatsImage('E:\MinD-MinE\20260701\old\OLD_WAVES_GFP_80.btf');

i_tstart = 0; 
N_ims = 240; %we expect to encounter/go through all the frames, i.e. all 240 frames

[im, branchData] = setup('/Users/heatherli491/Downloads/OLD_WAVES_GFP_80.btf', N_ims, i_tstart);

ibranch = 52;
topx = 5;
display = true;
shiftin = 100;
step = 100;

[sorted_scores, sorted_pixel_indices] = singlepixelswithbestperiodicity(singlepixelovertime, singlepixelallfilters, singlepixelscorrshift, input, branchData, im, N_ims, ibranch, step, topx, display, shiftin);


