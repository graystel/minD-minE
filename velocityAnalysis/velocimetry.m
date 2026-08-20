%% set up


bfr = BioformatsImage('C:\mydata\School\RoperLab\trains\matlab\OLD_WAVES_GFP_80.btf');
i_tstart = 0;
N_ims = 240;

if ~exist('im','var') | ~exist('imframe','var')

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
    
    [Ny,Nx] = size(imframe);

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


sigma = 5;
filteredMask = imgaussfilt(double(cleanMask),sigma);

im_filt5 = im;

for i_t = 1:N_ims
    thisIm = squeeze(im(:,:,i_t));
    thisIm(~cleanMask) = 0; % remove any pixel outside of the mask
    filteredIm = imgaussfilt(thisIm,sigma)./filteredMask; % cleans up the image a bit   
    filteredIm(~cleanMask) = 0;
    im_filt5(:,:,i_t) = filteredIm;    
end







% pick a specific branch to get intensity data from, run the code snippet
% at the end to see what each branch is

ibranch = 56;
thisbranchData = branchData{ibranch};

branchintensity = zeros(size(thisbranchData,1),N_ims);
inds = thisbranchData(:,3);
% output: first four columns are the x,y,index and s data

for i_t = 1:N_ims
    filteredIm = squeeze(im_filt5(:,:,i_t));
    branchintensity(:,i_t) = filteredIm(inds);
end

%% clean
close all;
%% UI set up

winOne = figure('Name', 'Dataset A Diagnostics');
tgOne  = uitabgroup(winOne);

% Create a Second Window
winTwo = figure('Name', 'Dataset B Diagnostics');
tgTwo  = uitabgroup(winTwo);

% Create a Third Window
winThree = figure('Name', 'Parameter Comparisons');
tgThree  = uitabgroup(winThree);


%% Functions
function snr_map = calculateSNR(zncc_matrix, exclusion_radius)
    num_points = size(zncc_matrix, 1);
    num_shifts = size(zncc_matrix, 2);
    snr_map = zeros(1, num_points);
    for k = 1:num_points
        zncc = zncc_matrix(k, :);
    
        % Find the main peak
        [R_max, max_idx] = max(zncc);
        noise_mask = true(1, num_shifts);
        start_ex = max(1, max_idx - exclusion_radius);
        end_ex   = min(num_shifts, max_idx + exclusion_radius);
        noise_mask(start_ex : end_ex) = false;       
    
        if sum(noise_mask) > 2
            mu_noise    = mean(zncc(noise_mask));
            sigma_noise = std(zncc(noise_mask));
        else
            mu_noise    = mean(zncc);
            sigma_noise = std(zncc);
        end
    
        snr_map(k) = (R_max - mu_noise) / (sigma_noise + eps);
    end
end

function velTab = plotVelocimetryQC(sUniform, I1, I2, pos, v_pair, frameA, frameB, tabGroup)
    
    velTab = uitab(tabGroup);
    velLay = tiledlayout(velTab, 3,1);
    % 1. Interpolate velocity onto sUniform grid
    v_full = interp1(pos, v_pair, sUniform, 'linear', 'extrap');
    
    % 2. Sample Frame 2 at (s + v) using sUniform
    I2_aligned = interp1(sUniform, I2, sUniform + v_full, 'linear', NaN);


    % --- Panel 1: Original Unshifted Waves ---
    ax1 = nexttile(velLay);
    plot(ax1, sUniform, I1, 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('Frame %d', frameA));
    hold(ax1, 'on');
    plot(ax1, sUniform, I2, 'g--', 'LineWidth', 1.5, 'DisplayName', sprintf('Frame %d', frameB));
    hold(ax1, 'off'); grid(ax1, 'on');
    ylabel(ax1, 'Intensity');
    title(ax1, sprintf('1. Original Wave Propagation: Frame %d vs Frame %d', frameA, frameB));
    legend(ax1, 'Location', 'northeast');

    % --- Panel 2: Shift-Corrected Overlay ---
    ax2 = nexttile(velLay);
    plot(ax2, sUniform, I1, 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf( 'Frame %d', frameA));
    hold(ax2, 'on');
    plot(ax2, sUniform, I2_aligned, 'g--', 'LineWidth', 1.5, 'DisplayName', sprintf( 'Frame %d ( Shift-Corrected)', frameB));
    hold(ax2, 'off'); grid(ax2, 'on');
    ylabel(ax2, 'Intensity');
    ylim(ax2, [min( I1)*0.9, max( I1)*1.1]); 
    title(ax2, '2. Shift-Corrected Overlay');
    legend(ax2, 'Location', 'northeast');

    % --- Panel 3: Spatial Velocity Profile ---
    ax3 = nexttile(velLay);

    plot(ax3, pos, v_pair, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'MarkerSize', 4);
    yline(ax3, 0, 'k--', 'Zero Shift'); grid(ax3, 'on')
    xlabel(ax3, 'Arc Length s (pixels)'); 
    ylabel(ax3, 'Velocity (px/frame)');
    title(ax3, '3. Calculated Spatial Velocity Profile v (s)');


    linkaxes([ax1,ax2,ax3], 'x');
end

function dsTab = plotDSDT(sUniform, I1, I2, pos, frameA, frameB, tabGroup)

    dsTab = uitab(tabGroup);
    velLay = tiledlayout(dsTab, 3,1);
    
    half_win = 2; 
    N = length(I1);
    slopes = zeros(1, N); % Pre-allocate for performance

    
    for i = 1:N
        % Determine safe local range so index doesn't go < 1 or > N
        start_idx = max(1, i - half_win);
        end_idx   = min(N, i + half_win);
    
        % Extract local neighborhood AND center x around 0
        x_local = sUniform(start_idx:end_idx) - sUniform(i); 
        y_local = I1(start_idx:end_idx);
    
        % Fit parabola and extract the centered derivative
        p = polyfit(x_local, y_local, 2);
    
        % The derivative at the shifted center (0) is just the p(2) coefficient
        slopes(i) = p(2);
    end
    
    
    % --- Panel 1: Original Unshifted Waves ---
    ax1 = nexttile(velLay);
    plot(ax1, sUniform, I1, 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('Frame %d', frameA));
    grid(ax1, 'on');
    ylabel(ax1, 'Intensity');
    title(ax1, sprintf('1. Original Wave Propagation: Frame %d vs Frame %d', frameA, frameB));
    legend(ax1, 'Location', 'northeast');
    
    % --- Panel 2: Shift-Corrected Overlay ---
    ax2 = nexttile(velLay);
    plot(ax2, sUniform, I1, 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf( 'Frame %d', frameA));
    hold(ax2, 'on');
    plot(ax2, sUniform, I2, 'g--', 'LineWidth', 1.5, 'DisplayName', sprintf( 'Frame %d ( Shift-Corrected)', frameB));
    hold(ax2, 'off'); grid(ax2, 'on');
    ylabel(ax2, 'Intensity');
    ylim(ax2, [min( I1)*0.9, max( I1)*1.1]); 
    title(ax2, '2. Shift-Corrected Overlay');
    legend(ax2, 'Location', 'northeast');
    
    % --- Panel 3: Spatial Velocity Profile ---
    ax3 = nexttile(velLay);
    
    plot(ax3, sUniform, slopes, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'MarkerSize', 4);
    yline(ax3, 0, 'k--', 'Zero Shift'); grid(ax3, 'on')
    xlabel(ax3, 'Arc Length s (pixels)'); 
    ylabel(ax3, 'Velocity (px/frame)');
    title(ax3, '3. Calculated Spatial Velocity Profile v (s)');
    
    
    linkaxes([ax1,ax2,ax3], 'x');
end

%% section
sRaw = thisbranchData(:, 4); 
% 1. Strip duplicate arc lengths to prevent interp1 function from crashing
[sUnique, uniqueIdx] = unique(sRaw, 'stable');
branchintensityUnique = branchintensity(uniqueIdx, :);
ds = 1.0; 
sUniform = (min(sUnique) : ds : max(sUnique))'; 

% 2. Interpolate intensity profiles across all frames
branchintensityUniform = interp1(sUnique, branchintensityUnique, sUniform, 'pchip');

% 3. Set velocimetry parameters
window_size  = 101;   
search_range = 12; 
step_size    = 1;
dx           = ds;           
dt           = 1;            
N_ims        = size(branchintensityUniform, 2);
exclusion_radius = 5; % Needed for the SNR calculation

% 4. Compute velocity across all frame pairs
% Use the refactored function for the first pair to get exact array sizes
[v_first, s_pos, zncc_first] = estimateVelocity1Drefactor(...
    branchintensityUniform(:,1), branchintensityUniform(:,2), ...
    window_size, search_range, step_size, dx, dt);

% Calculate the first SNR map from the ZNCC matrix
snr_first = calculateSNR(zncc_first, exclusion_radius);

% Preallocate arrays
num_points    = length(v_first);
vel_map       = zeros(num_points, N_ims - 1);
snr_map       = zeros(num_points, N_ims - 1);

vel_map(:, 1) = v_first;
snr_map(:, 1) = snr_first;

% Loop through the rest of the frames
for i_frame = 2:(N_ims - 1)
    I1_loop = branchintensityUniform(:, i_frame);
    I2_loop = branchintensityUniform(:, i_frame + 1);

    % 1. Get Velocity and ZNCC Matrix
    [v_temp, ~, zncc_temp] = estimateVelocity1Drefactor(...
        I1_loop, I2_loop, window_size, search_range, 1, dx, dt);

    % 2. Get SNR from the ZNCC Matrix
    snr_temp = calculateSNR(zncc_temp, exclusion_radius);

    % 3. Store both safely
    vel_map(:, i_frame) = v_temp;
    snr_map(:, i_frame) = snr_temp; 
end








%% Code snippet for ZNCC Best-Fit & SNR Plot 
% function plotSNRDiagnostics(sUniform, I1, I2, positions, snr_array, target_pos, window_size, search_range, exclusion_radius)
% % PLOTSNRDIAGNOSTICS Visualizes the ZNCC curve for a single point and the SNR profile for the whole hypha.
% %
% % Inputs:
% %   sUniform   - Full spatial grid of the hypha
% %   I1, I2     - 1D intensity arrays for the two frames being compared
% %   positions  - The spatial coordinates where SNR was calculated
% %   snr_array  - The pre-calculated SNR values for this frame pair
% %   target_pos - The physical position (arc length) to deep-dive into
% %   window_size, search_range - Velocimetry parameters
% %   exclusion_radius - (Optional) Pixels to exclude around peak. Default: 5
% 
% if nargin < 9
%     exclusion_radius = 5;
% end
% 
% half_win = floor(window_size / 2);
% shifts   = -search_range : search_range;
% 
% % --- PART 1: Compute ZNCC on-the-fly for the single deep-dive point ---
% % Find the index in sUniform closest to the requested target_pos
% [~, target_idx] = min(abs(sUniform - target_pos)); 
% 
% % Re-run just this one point to get the curve for the plot
% zncc = computeZNCC(I1, I2, target_idx, half_win, shifts);
% [best_shift_px, R_max, max_idx] = findSubpixelPeak(zncc, shifts);
% 
% % Recreate the noise mask for plotting
% noise_mask = true(1, length(zncc));
% start_ex   = max(1, max_idx - exclusion_radius);
% end_ex     = min(length(zncc), max_idx + exclusion_radius);
% noise_mask(start_ex : end_ex) = false;
% 
% noise_points = zncc(noise_mask);
% mu_noise     = mean(noise_points);
% target_snr   = calculateSNR(zncc, R_max, max_idx);
% 
% % --- PART 2: Identify High Confidence Regions ---
% high_conf = snr_array > 3.0;
% 
% % --- PLOTTING ---
% fig = figure('Name', 'SNR Diagnostics', 'Position', [150, 150, 800, 650]);
% clf(fig);
% 
% % Panel 1: The ZNCC Curve for the specific target point
% subplot(2, 1, 1);
% plot(shifts, zncc, 'k-', 'LineWidth', 1.2, 'DisplayName', 'Full ZNCC Curve'); hold on;
% plot(shifts(noise_mask), noise_points, 'co', 'MarkerFaceColor', 'c', 'DisplayName', 'Noise Points');
% plot(best_shift_px, R_max, 'r*', 'MarkerSize', 12, 'LineWidth', 2, ...
%     'DisplayName', sprintf('Best Fit Peak (%.2f)', R_max));
% yline(mu_noise, 'm--', 'LineWidth', 1.5, ...
%     'DisplayName', sprintf('Mean Noise Floor (%.2f)', mu_noise));
% grid on;
% xlabel('Shift (pixels)');
% ylabel('ZNCC Score');
% title(sprintf('1. ZNCC Best-Fit at s = %.1f px (SNR = %.2f)', sUniform(target_idx), target_snr));
% legend('Location', 'northeast');
% hold off;
% 
% % Panel 2: The full SNR Profile along the hypha
% subplot(2, 1, 2);
% % Plotting directly against 'positions' since SNR isn't calculated at the margins
% plot(positions, snr_array, 'k-', 'LineWidth', 1.2, 'DisplayName', 'Local SNR'); hold on;
% yline(3.0, 'r--', 'Confidence Threshold (SNR = 3.0)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
% 
% % Highlight high confidence regions in green
% plot(positions(high_conf), snr_array(high_conf), 'g.', 'MarkerSize', 10, 'DisplayName', 'High Confidence Data');
% 
% grid on;
% xlim([min(sUniform), max(sUniform)]); % Keep X-axis scaled to the full hypha
% xlabel('Arc Length s (pixels)');
% ylabel('SNR Value');
% title('2. Full Hypha SNR Profile');
% legend('Location', 'northeast');
% hold off;
% end
% 
% % Choose which frame pair to look at
% frameA = 10;
% frameB = 11;
% 
% % Grab the specific data for that pair
% I1_test = branchintensityUniform(:, frameA);
% I2_test = branchintensityUniform(:, frameB);
% snr_test = snr_map(:, frameA); 
% 
% % Pick a physical location to deep-dive (e.g., the exact middle)
% middle_pos = max(sUniform) / 2;
% 
% % Call the function!
% plotSNRDiagnostics(sUniform, I1_test, I2_test, s_pos, snr_test, middle_pos, 81, 12);
% 


%% plotVelocimetry Section

frameA = 8;
frameB = frameA+1;
I1 = branchintensityUniform(:, frameA);
I2 = branchintensityUniform(:, frameB);
velOne = plotVelocimetryQC(sUniform, I1, I2, s_pos, vel_map(:, frameA), frameA, frameB, tgOne);
frameA = 20;
frameB = frameA+1;
I1 = branchintensityUniform(:, frameA);
I2 = branchintensityUniform(:, frameB);
velOne = plotVelocimetryQC(sUniform, I1, I2, s_pos, vel_map(:, frameA), frameA, frameB, tgOne);
velTwo = plotVelocimetryQC(sUniform, I1, I2, s_pos, vel_map(:, frameA), frameA, frameB, tgTwo);
plotDSDT(sUniform, I1, I2, s_pos, frameA, frameB, tgTwo);

%drawnow;


%added 




% % Setup parameters
% frameA = 10;
% frameB = 11;
% I1 = branchintensityUniform(:, frameA);
% I2 = branchintensityUniform(:, frameB);
% 
% window_size  = 81;   
% search_range = 12; 
% half_win = floor(window_size / 2);
% shifts   = -search_range : search_range;
% margin   = half_win + search_range;
% 
% % --- PART 1: Single-Point Deep Dive ---
% target_idx = round(length(sUniform) / 2); 
% 
% % 1. Compute ZNCC, find peak, and calculate SNR using helpers
% zncc = computeZNCC(I1, I2, target_idx, half_win, shifts);
% [best_shift_px, R_max, max_idx] = findSubpixelPeak(zncc, shifts);
% snr_val = calculateSNR(zncc, R_max, max_idx);
% 
% % 2. Recreate the mask just for the plot visuals (matching calculateSNR)
% exclusion_radius = 5; 
% noise_mask = true(1, length(zncc));
% start_ex   = max(1, max_idx - exclusion_radius);
% end_ex     = min(length(zncc), max_idx + exclusion_radius);
% noise_mask(start_ex : end_ex) = false;
% 
% noise_points = zncc(noise_mask);
% mu_noise     = mean(noise_points);
% 
% % 3. Plot ZNCC Best Fit and Noise Floor
% figure(11); clf;
% plot(shifts, zncc, 'k-', 'LineWidth', 1.2, 'DisplayName', 'Full ZNCC Curve'); hold on;
% plot(shifts(noise_mask), noise_points, 'co', 'MarkerFaceColor', 'c', 'DisplayName', 'Noise Points');
% plot(best_shift_px, R_max, 'r*', 'MarkerSize', 12, 'LineWidth', 2, ...
%     'DisplayName', sprintf('Best Fit Peak (%.2f)', R_max));
% yline(mu_noise, 'm--', 'LineWidth', 1.5, ...
%     'DisplayName', sprintf('Mean Noise Floor (%.2f)', mu_noise));
% grid on;
% xlabel('Shift (pixels)');
% ylabel('ZNCC Score');
% title(sprintf('ZNCC Best-Fit at s = %.1f px (SNR = %.2f)', sUniform(target_idx), snr_val));
% legend('Location', 'northeast');
% hold off;
% 
% % --- PART 2: Full Hypha SNR Profile ---
% snr_profile = zeros(1, length(sUniform));
% 
% % Loop through spatial locations along the hypha
% for idx = (margin + 1) : (length(sUniform) - margin)
% 
%     % Use helpers to compute everything in 3 lines
%     zncc_curve = computeZNCC(I1, I2, idx, half_win, shifts);
%     [~, R_max, max_idx] = findSubpixelPeak(zncc_curve, shifts);
% 
%     snr_profile(idx) = calculateSNR(zncc_curve, R_max, max_idx);
% end
% 
% % Chunk hypha into High SNR (> 3.0) and Low SNR (< 3.0) regions
% high_confidence_chunks = snr_profile > 3.0;