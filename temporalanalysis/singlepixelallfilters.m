function [filteredmed, filteredsgolay, filteredgauss] = singlepixelallfilters(pixelintensity, ipixel, ibranch, bwindow, a, display)
    arguments
        pixelintensity {mustBeVector}
        ipixel (1, 1) int16
        ibranch = 52
        bwindow = 11
        a = 1
        display = false
    end
    

    %set up
    pixelsignalfilt = pixelintensity(5:end);
    pixelsignalfilt = pixelsignalfilt - mean(pixelsignalfilt);
    pixelsignalfilt = pixelsignalfilt / std(pixelsignalfilt);
    pixelsignalfilt = double(pixelsignalfilt);

    % FILTERING
    filteredmed = medfilt1(pixelsignalfilt, 10);
    filteredsgolay = sgolayfilt(pixelsignalfilt, 3, 11);
    b = gausswin(bwindow); %gausswin is default gaussian filter, input of gausswin is length of window for smoothing 
    % (should be odd for center point to emerge) 
    b = b / sum(b); %normalization
    filteredgauss = filtfilt(b, a, pixelsignalfilt); %here, b = numerator coeffs (coeffs onto the input signal), 
    %a = denom coeffs (coeffs on the output signal), a = 1 for gaussian almost always
    %filtfilt is like filter, but doesn't shift peaks which filter will do unintentionally 

    if display == true
        figure;
        plot(1:N_ims, pixelsignalfilt, 1:N_ims, filteredmed, 1:N_ims, filteredsgolay, 1:N_ims, filteredgauss);
        xlabel("frame");
        ylabel("intensity");
        title("pixel intensity over time, branch " + ibranch + ", pixel " + ipixel + "");
        legend("fft w/ bp", "original", "median filter", "sgolay filter", "gaussian filter");
        legend("boxon");
    end
end

