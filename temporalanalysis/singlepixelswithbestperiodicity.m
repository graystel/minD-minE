function [sorted_scores, sorted_pixel_indices] = singlepixelswithbestperiodicity(singlepixelovertime, singlepixelallfilters, singlepixelscorrshift, input, branchData, im, N_ims, ibranch, step, topx, display, shiftin)
    
    %shiftin = amount i want to skip to start; i wanna start into the branch b/c
    %i suspect edges of branches will be kinda weird and not worth looking at (i skipped 100)
    if nargin < 12 || isempty(shiftin)
        shiftin = 100;
    end
    if nargin < 11 || isempty(display)
        display = true;
    end
    if nargin < 10 || isempty(topx)
        topx = 10;
    end
    if nargin < 9 || isempty(step)
        step = 10;
    end
    if nargin < 8 || isempty(ibranch)
        ibranch = 52;
    end
    if nargin < 7 || isempty(N_ims)
        N_ims = 240;
    end

    num_pixels = length(input);
    pixel_indices = shiftin:step:num_pixels; %skipping around by 10 px bc Gaussian will muddle across pixels
    peak_scores = zeros(2, length(pixel_indices));

    for i = 1:length(pixel_indices)
        %getting pixel info
        ipixel = pixel_indices(i);
        pixelintensity = singlepixelovertime(branchData, im, N_ims, ibranch, ipixel, false);
        [~, ~, filteredgauss] = singlepixelallfilters(pixelintensity, ipixel, ibranch, bwindow, a, false);
        [r, lags] = singlepixelscorrshift(filteredgauss, ibranch, ipixel, false); %coeff normalizes lag 0 to 1.0
        
        %store the max
        pos_r = r(lags > 10);
        pks = findpeaks(pos_r);
        if ~isempty(pks)
            peak_scores(1, i) = max(pos_r);
            peak_scores(2, i) = ipixel;
        else
            peak_scores(1, i) = 0;
            peak_scores(2, i) = 0;
        end
    end

    [~, sort_idx] = sort(peak_scores(1,:), "descend");
    sorted_peak_scores = sort(peak_scores, sort_idx);

    sorted_scores = sorted_peak_scores(1,:);
    sorted_pixel_indices = sorted_peak_scores(2,:);
    
    %display
    if display == true
        significantpixels = sorted_pixel_indices;
        for k = 1:topx
        fprintf('Rank %d Pixel: Index %d (Score = %.3f)\n', ...
            k, sorted_pixel_indices(k), sorted_scores(k));
        significantpixels(k) = sorted_pixel_indices(k);
        end
    end
end