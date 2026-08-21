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

