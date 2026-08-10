function [sorted_scores, sorted_pixel_indices] = singlepixelswithbestperiodicity(singlepixelovertime, singlepixelallfilters, singlepixelxcorrshift, branchData, im, N_ims, ibranch, step, topx, display, shiftin, bwindow, a)
    arguments 
        singlepixelovertime
        singlepixelallfilters
        singlepixelxcorrshift
        branchData
        im
        N_ims = 240
        ibranch = 52
        step = 10
        topx = 20
        display = true
        shiftin = 100
        bwindow = 11
        a = 1
    end

    %shiftin = amount i want to skip to start; i wanna start into the branch b/c
    %i suspect edges of branches will be kinda weird and not worth looking at (i skipped 100)

    num_pixels = size(branchData{ibranch}, 1);
    disp(num_pixels);
    pixel_indices = shiftin:step:num_pixels; %skipping around by 10 px bc Gaussian will muddle across pixels
    peak_scores = zeros(2, length(pixel_indices));

    for i = 1:length(pixel_indices)
        %getting pixel info
        ipixel = pixel_indices(i);
        pixelintensity = singlepixelovertime(branchData, im, ipixel, N_ims, ibranch, false);
        [~, ~, filteredgauss] = singlepixelallfilters(pixelintensity, ipixel, ibranch, bwindow, a, false);
        [r, lags] = singlepixelxcorrshift(filteredgauss, ibranch, ipixel, false); %coeff normalizes lag 0 to 1.0
        
        %store the max
        pos_r = r(lags > 10);
        pks = findpeaks(pos_r);
        if ~isempty(pks)
            peak_scores(1, i) = max(pks);
            peak_scores(2, i) = ipixel;
        else
            peak_scores(1, i) = 0;
            peak_scores(2, i) = 0;
        end
    end

    [~, sort_idx] = sort(peak_scores(1,:), "descend");
    sorted_peak_scores = peak_scores(:, sort_idx);

    sorted_scores = sorted_peak_scores(1,:);
    sorted_pixel_indices = sorted_peak_scores(2,:);

    %CHECKPT 1.5
    fprintf("checkpt 1.5");
    disp(size(sorted_scores));

    
    %display
    if display == true
        significantpixels = sorted_pixel_indices;
        for k = 1:topx
            fprintf('Rank %d Pixel: Index %d (Score = %.3f)\n', k, sorted_pixel_indices(k), sorted_scores(k));
            significantpixels(k) = sorted_pixel_indices(k);
        end
    end
end