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
plot(1:N_ims, pixelintensity(1,5:end));
xlabel('Frame');
ylabel('Intensity');
title('Pixel intensity over time');

%% BANDPASS FILTER
pixelsignal = pixelintensity(5:end);
pixelsignal = pixelsignal - mean(pixelsignal);
pixelsignal = pixelsignal / std(pixelsignal);
fs = 1/dt;
y = fft(pixelsignal);
l = length(pixelsignal);
p2 = abs(y/l);
p1 = p2(1:floor(l/2)+1);
p1(2:end-1) = 2*p1(2:end-1);
f = fs*(0:floor(l/2))/l;
figure(4);
plot(f, p1);
xlabel('freq (hz)');
ylabel('amplitude');
title('fft of the pixel intensity');
filteredsignal = bandpass(pixelsignal, [0.02 0.04], 1/dt);
figure(5);
plot(1:N_ims, pixelsignal, "DisplayName", "norm pixel intensity");
xlabel('Frame');
ylabel('Intensity');
title('Pixel intensity over time');
hold on;
plot(1:N_ims, filteredsignal, "DisplayName", "filtered pixel intensity");
xlabel('Frame');
ylabel('Filtered intensity');
title('Filtered pixel intensity over time');








