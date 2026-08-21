function [velocities, positions, zncc_matrix] = estimateVelocity1Drefactor(I1, I2, window_size, search_range, stepSize, dx, dt)
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
%   zncc_matrix  - Array of Zero-Mean Normalized Cross-Correlation values

    
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



