function pixelintensity = singlepixelovertime(branchData, im, ipixel, N_ims, ibranch,  display)
    arguments 
        branchData 
        im
        ipixel 
        N_ims = 240
        ibranch = 52
        display = true
    end

    thispixelData = branchData{ibranch}(ipixel,:);
    x_coord = thispixelData(1); 
    y_coord = thispixelData(2);

    pixelintensity = zeros(1, N_ims+4);
    for i_t = 1:N_ims
        im_smooth = squeeze(im(:,:,i_t));
        pixelintensity(1, i_t+4) = im_smooth(x_coord, y_coord);
    end
    if display == true
        figure;
        plot(1:N_ims, pixelintensity);
        xlabel("frames (# of frames)");
        ylabel("intensity")
        title("signal intensity over time of the " + ipixel + "th pixel of branch " + ibranch + "");
    end
end