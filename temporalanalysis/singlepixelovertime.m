function pixelintensity = singlepixelovertime(branchData, im, N_ims, ibranch, ipixel, display)
    if nargin < 6 || isempty(display)
        display = true;
    end
    if nargin < 5 || isempty(ipixel)
        ipixel = 561;
    end
    if nargin < 4 || isempty(ibranch)
        ibranch = 52;
    end 
    if nargin < 3 || isempty(N_ims)
        N_ims = 240;
    end 
    figure;
    
    thispixelData = branchData{ibranch}(ipixel,:);
    pixelintensity = zeros(1, N_ims+4);
    pixelintensity(:,1:4) = thispixelData;
    for i_t = 1:N_ims
        im_smooth = squeeze(im(:,:,i_t));
        pixelintensity(1, i_t+4) = im_smooth(pixelintensity(1, 3));
    end
    if display == true
        figure;
        plot(1:N_ims, pixelintensity);
        xlabel("frames (# of frames)");
        ylabel("intensity")
        title("signal intensity over time of the " + ipixel + "th pixel of branch " + ibranch + "");
    end
end