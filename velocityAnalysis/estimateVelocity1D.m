function [velocities, positions, zncc_matrix] = estimateVelocity1D(I1, I2, window_size, search_range, stepSize, dx, dt)
% estimateVelocity1D Estimates velocity from two 1D intensity distributions.
%
% Inputs:
%   I1, I2       - 1D arrays of intensity distributions (must be same size)
%   window_size  - Length of the snippet/template to track (pixels/indices)
%   search_range - Max expected movement in either both direction (pixels)
%   stepSize     - How far the window moves each time
%   dx           - Physical distance per pixel/index (e.g., meters per pixel)
%   dt           - Time interval between I1 and I2 (e.g., seconds)
%
% Outputs:
%   velocities   - Array of estimated velocities at each position
%   positions    - Physical spatial coordinates where velocity was measured
%   zncc_matrix  - Zero-Mean Normalized Cross-Correlation

    
    arguments
        I1 double {mustBeVector}
        I2 double {mustBeVector}
        window_size (1,1) double  = 81
        search_range (1,1) double = 20
        stepSize (1,1) double = -1.0
        dx (1,1) double = 1.0
        dt (1,1) double  = 1.0
    end
    
    % Force arrays to both be rows
    I1 = I1(:).';
    I2 = I2(:).';
    % Force window size = odd so there's exact center pixel for tracking
    if mod(window_size, 2) == 0
        window_size = window_size + 1;
    end
    half_win = floor(window_size / 2);
   
    % Define the center points for where windows go to sample that area
    % Default step size = half window size so regions overlap
    if stepSize <= 0
        step_size = half_win;
    else
        step_size = stepSize;
    end

    N             = length(I1);
    center_points = (half_win + search_range + 1) : step_size : (N - half_win - search_range);
    if isempty(center_points)
        warning('estimateVelocity1D:StepSizeTooLarge', ...
            'stepSize (%g) exceeds the valid signal region. Returning empty arrays.', step_size);
        velocities = [];
        positions  = [];
        return;
    end
    num_points    = length(center_points);

    %outputs
    velocities  = zeros(1, num_points); %sci convention to keep velocity singular not plural (delete later)
    positions   = center_points * dx; % Convert indices to physical units (we don't care about this specific step)

    % Array of pixel offsets to test
    shifts = -search_range : search_range;
    zncc_matrix = zeros(num_points, length(shifts));

    % --- MAIN MAPPING LOOP ---
    for k = 1:num_points
        idx = center_points(k);
        
        zncc = computeZNCC(I1,I2,idx,half_win,shifts);

        [best_shift_pixels, ~,~] = findSubpixelPeak(zncc, shifts);

        % Calculate physical velocities ( v = dx/dt )
        velocities(k) = (best_shift_pixels * dx) / dt;
        zncc_matrix(k,:) = zncc;
    end
end





% ------------------ HELPER FUNCTIONS ------------------

%computes the zero mean normalized cross-correlation
function zncc = computeZNCC(I1, I2, idx, half_win, shifts) %change to integer shiftRange?

    % Extract template from I1
    shiftRange = length(shifts);
    template   = I1(idx - half_win : idx + half_win); % this will be compared with the second intensity trace
    t_zm = normalizeSignal(template);
    zncc = zeros(1, shiftRange);
    
    
    % Test every possible shift within the search range
    for s = 1:shiftRange
        shift = shifts(s);
    
        % Extract corresponding snippet from I2
        i2_start = idx - half_win + shift;
        i2_end   = idx + half_win + shift;
        snippet2 = I2(i2_start : i2_end); % what is compared to the template
    
        s2_zm     = normalizeSignal(snippet2); 
    
        % Calculate Normalized Cross-Correlation
        zncc(s) = sum(t_zm .* s2_zm) / length(template);
    end
end

%normalizes a signal to have zero mean and unit variance
function normSig = normalizeSignal(sig) 

    sig_mean = mean(sig);
    sig_std = std(sig);

    if sig_std == 0, sig_std = 1; end % Prevent division by zero
    normSig = (sig - sig_mean) / sig_std;
end

%Finds peak correlation and applies parabolic sub-pixel interpolation
function [best_shift_pixels, R_max, max_idx] = findSubpixelPeak(zncc, shifts)

    
    [R_max, max_idx]  = max(zncc);
    best_shift_pixels = shifts(max_idx);
    
    % Sub-pixel interpolation (Parabolic fit)
    if max_idx > 1 && max_idx < length(shifts)
        c_left   = zncc(max_idx - 1);
        c_center = zncc(max_idx);
        c_right  = zncc(max_idx + 1);

        % Shift correction using peak fitting
        sub_pixel_correction = (c_left - c_right) / (2 * (c_left - 2*c_center + c_right));
        best_shift_pixels = best_shift_pixels + sub_pixel_correction;
    end
end

% Calculates Peak-to-Sidelobe Ratio by masking out the main peak
