%% code snippet for calculating velocities


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
[v_first, s_pos, zncc_first] = estimateVelocity1D(...
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
    [v_temp, ~, zncc_temp] = estimateVelocity1D(...
        I1_loop, I2_loop, window_size, search_range, 1, dx, dt);

    % 2. Get SNR from the ZNCC Matrix
    snr_temp = calculateSNR(zncc_temp, exclusion_radius);

    % 3. Store both safely
    vel_map(:, i_frame) = v_temp;
    snr_map(:, i_frame) = snr_temp; 
end


% Set up time axes (X-axis)
frames_raw = 1:size(branchintensityUniform, 2); % For intensity (N frames)
frames_vel = 1:size(vel_map, 2);                % For velocity (N-1 frame pairs)

figure('Name', 'Kymograph Analysis', 'Position', [100, 100, 1000, 800]);

% --- Panel 1: Raw Intensity Kymograph ---
subplot(2, 1, 1);
imagesc(frames_vel, sUniform, branchintensityUniform);
colormap(gca, 'turbo'); 
cb1 = colorbar;
cb1.Label.String = 'Brightness';
title('Brightness Kymograph');
xlabel('Time (Frames)');
ylabel('Arc Length s (pixels)');
set(gca, 'YDir', 'normal'); 
%lines
y = (570) + 9.15*(t-30);
y2 = (570) + 1.99*(t-34);
t = linspace(20,44,2);
hold on;
plot(t, y, 'LineWidth', 2);
plot(t, y2, 'LineWidth', 2);
hold off;
%lines
subplot(2, 1, 2);
imagesc(frames_vel, s_pos, vel_map);
colormap(gca, 'jet'); 
cb1 = colorbar;
cb1.Label.String = 'Velocity(px)';
title('Vel Kymograph');
xlabel('Time (Frames)');
ylabel('Arc Length s (pixels)');
set(gca, 'YDir', 'normal'); 

%% section
% 1. Set your target index and local window size
       % The point where you want the slope
       I1 = branchintensityUniform(:, frameA);
       I2 = branchintensityUniform(:, frameB);
% 
% half_window = 2;   % Uses 5 points left & 5 points right (11 points total)
% for pt = I1
% 
%     % 2. Extract the local neighborhood around your point
%     x_local = x(idx - half_window : idx + half_window);
%     y_local = y(idx - half_window : idx + half_window);
%     % 3. Fit local parabola: y = p1*x^2 + p2*x + p3
%     p = polyfit(x_local, y_local, 2);
%     % 4. Exact tangent slope (derivative) at x(idx): d/dx = 2*p1*x + p2
%     slope_at_point = 2 * p(1) * x(idx) + p(2);
% end
% 






function plotVelocimetryQC(sUniform, I1, I2, pos, v_pair, frameA, frameB)
% 1. Interpolate velocity onto sUniform grid
v_full = interp1(pos, v_pair, sUniform, 'linear', 'extrap');
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

% 2. Sample Frame 2 at (s + v) using sUniform
I2_aligned = interp1(sUniform, I2, sUniform + 9.15, 'linear', NaN);

fig = figure('Name', 'Wave Alignment Diagnostics', 'Position', [100, 100, 900, 800]);
% --- Panel 1: Original Unshifted Waves ---
subplot(3, 1, 1);
plot(sUniform, slopes, 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('Frame %d', frameA));
hold on;
%plot(sUniform, I2, 'g--', 'LineWidth', 1.5, 'DisplayName', sprintf('Frame %d', frameB));
hold off; 
grid on;
ylabel('Intensity');
title(sprintf('1. Original Wave Propagation: Frame %d vs Frame %d', frameA, frameB));
legend('Location', 'northeast');

% --- Panel 2: Shift-Corrected Overlay ---
subplot(3, 1, 2);
plot(sUniform, I1, 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('Frame %d', frameA));
hold on;
plot(sUniform, I2_aligned, 'g--', 'LineWidth', 1.5, 'DisplayName', sprintf('Frame %d (Shift-Corrected)', frameB));
hold off; grid on;
ylabel('Intensity');
ylim([min(I1)*0.9, max(I1)*1.1]); 
title('2. Shift-Corrected Overlay');
legend('Location', 'northeast');

% --- Panel 3: Spatial Velocity Profile ---
subplot(3, 1, 3);
plot(pos, v_pair, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'MarkerSize', 4);
yline(0, 'k--', 'Zero Shift'); grid on;
xlabel('Arc Length s (pixels)'); 
ylabel('Velocity (px/frame)');
title('3. Calculated Spatial Velocity Profile v(s)');


linkaxes(findobj(fig, 'type', 'axes'), 'x');
end
frameA = 0;
frameB = 0;


my_list = [30];
for item = my_list
    frameA = item;
    frameB = frameA+1;
    I1 = branchintensityUniform(:, frameA);
    I2 = branchintensityUniform(:, frameB);
    snr_test = snr_map(:, frameA); 
    vel_temp = vel_map(:, frameA);
    plotVelocimetryQC(sUniform, I1, I2, s_pos, vel_temp, frameA, frameB);
end



% 
% % --- Panel 2: Velocity Kymograph (Filtered by SNR) ---
% subplot(2, 1, 2);
% % 1. Filter the data: Set any velocity with an SNR < 3.0 to NaN (Not a Number)
% snr_threshold = 3.0;
% vel_masked = vel_map;
% vel_masked(snr_map < snr_threshold) = NaN;
% 
% % 2. Plot the velocities
% h = imagesc(frames_vel, s_pos, vel_masked);
% 
% % 3. Apply the transparency mask
% set(h, 'AlphaData', ~isnan(vel_masked)); 
% 
% % 4. Set the background axis color to dark gray so the transparent (bad) areas stand out
% set(gca, 'Color', [0.2 0.2 0.2]); 
% colormap(gca, 'parula');
% cb2 = colorbar;
% cb2.Label.String = 'Velocity (px/frame)';
% 
% % Optional: Center color limits around zero
% % clim([-max(abs(vel_masked(:))), max(abs(vel_masked(:)))]); 
% 
% title(sprintf('2. Local Velocity Kymograph (High Confidence: SNR > %.1f)', snr_threshold));
% xlabel('Time (Frame Pairs)');
% ylabel('Arc Length s (pixels)');
% set(gca, 'YDir', 'normal');
% 
% % Link the X and Y axes
% linkaxes(findobj(gcf, 'type', 'axes'), 'xy');