function significantpixelfeatures = gethilberttransformparams(hilberttransform, pixelsofinterest, branchData, im, N_ims, ibranch, bwindow, a)
    arguments 
        hilberttransform
        pixelsofinterest
        branchData
        im
        N_ims = 240
        ibranch = 56
        bwindow = 11
        a = 1
    end
    
    topx = length(pixelsofinterest);
    significantpixelfeatures = zeros(topx, N_ims); %2D matrix / dict with the pixel and its frequencies across all N_im frames
    
    for k = 1:topx
        ipixel = pixelsofinterest(k);
        thispixelData = branchData{ibranch}(ipixel, :);
        x_coord = thispixelData(1);
        y_coord = thispixelData(2);
        
        pixelintensity = zeros(1, N_ims);
        for i_t = 1:N_ims
            im_smooth = im(:,:,i_t);
            pixelintensity(i_t) = im_smooth(y_coord, x_coord); %row = Y, col = X
        end
        
        pixelsignalfilt = pixelintensity - mean(pixelintensity);
        if std(pixelsignalfilt) > 0
            pixelsignalfilt = pixelsignalfilt / std(pixelsignalfilt);
        end
        b = gausswin(bwindow);
        b = b / sum(b);
        filteredgauss = filtfilt(b, a, pixelsignalfilt);
        
        params = hilberttransform(filteredgauss); %using the hilbert transform
        
        %get the frequency of intensity experienced by each pixel
        significantpixelfeatures(k, 1) = ipixel;
        for i=1:length(params.Frequency_Hz)
            significantpixelfeatures(k, i+1) = params.Frequency_Hz(i);
        end
    end

end
