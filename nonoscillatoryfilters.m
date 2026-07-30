setup
hyphaelabeling
%% SINGLE PIXEL OVER TIME
figure(3);
ibranch = 52;
ipixel = 561;
thispixelData = branchData{ibranch}(ipixel,:);
pixelintensity = zeros(1, N_ims+4);
pixelintensity(:,1:4) = thispixelData;
for i_t = 1:N_ims
   im_smooth = squeeze(im(:,:,i_t));
   pixelintensity(1, i_t+4) = im_smooth(pixelintensity(1, 3));
end

%% NON-OSCILLATORY FILTERS
pixelsignalfilt = pixelintensity(5:end);
pixelsignalfilt = pixelsignalfilt - mean(pixelsignalfilt);
pixelsignalfilt = pixelsignalfilt / std(pixelsignalfilt);
filteredmed = medfilt1(pixelsignalfilt, 10);
filteredsgolay = sgolayfilt(pixelsignalfilt, 3, 11);
b = gausswin(11); %gausswin is default gaussian filter, input of gausswin is length of window for smoothing
% (should be odd for center point to emerge)
b = b / sum(b); %normalization
a = 1;
filteredgauss = filtfilt(b, a, pixelsignalfilt); %here, b = numerator coeffs (coeffs onto the input signal),
%a = denom coeffs (coeffs on the output signal), a = 1 for gaussian almost always
%filtfilt is like filter, but doesn't shift peaks which filter will do unintentionally
hold on;
figure(5);
plot(1:N_ims, filteredsignal, 1:N_ims, pixelsignalfilt, 1:N_ims, filteredmed, 1:N_ims, filteredsgolay, 1:N_ims, filteredgauss);
xlabel("frame");
ylabel("intensity");
title("pixel intensity over time, branch " + ibranch + ", pixel " + ipixel + "");
legend("fft w/ bp", "original", "median filter", "sgolay filter", "gaussian filter");
legend("boxon");





