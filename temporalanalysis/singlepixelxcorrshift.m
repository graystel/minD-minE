function [r, lags] = singlepixelxcorrshift(filtered, ibranch, ipixel, display)
    %XCORR SHIFTING
    filteredcenter = filtered - mean(filtered); %normalization
    [r, lags] = xcorr(filteredcenter, "coeff"); %coeff normalizes lag 0 to 1.0

    if display == true
        figure(6);
        plot(lags, r);
        xlabel("lag (# of frames)");
        ylabel ("autocorrelation");
        title("X-shift correlation, branch " + ibranch + ", pixel" + ipixel + "");
    end
end


