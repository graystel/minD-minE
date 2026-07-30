setup;
hyphaelabeling;


%% SINGLE PIXEL OVER TIME
figure(3);
ibranch = 52;
ipixel = 815;
thispixelData = branchData{ibranch}(ipixel,:);
pixelintensity = zeros(1, N_ims+4);
pixelintensity(:,1:4) = thispixelData;
for i_t = 1:N_ims
   im_smooth = squeeze(im(:,:,i_t));
   pixelintensity(1, i_t+4) = im_smooth(pixelintensity(1, 3));
end

%% GAUSSIAN FILTERING

pixelsignalfilt = pixelintensity(5:end);

pixelsignalfilt = pixelsignalfilt - mean(pixelsignalfilt);
pixelsignalfilt = pixelsignalfilt / std(pixelsignalfilt);

b = gausswin(11); %gausswin is default gaussian filter, input of gausswin is length of window for smoothing 
% (should be odd for center point to emerge) 
b = b / sum(b); %normalization
a = 1;

filteredgauss = filtfilt(b, a, pixelsignalfilt); %here, b = numerator coeffs (coeffs onto the input signal), 
%a = denom coeffs (coeffs on the output signal), a = 1 for gaussian almost always
%filtfilt is like filter, but doesn't shift peaks which filter will do unintentionally 

%% XCORR SHIFTING

filteredgausszerocenter = filteredgauss - mean(filteredgauss); %normalization
[r, lags] = xcorr(filteredgausszerocenter, "coeff"); %coeff normalizes lag 0 to 1.0

figure(6);
plot(lags, r);
xlabel("lag (# of frames)");
ylabel ("autocorrelation");
title("X-shift correlation, branch " + ibranch + ", pixel" + ipixel + "");


