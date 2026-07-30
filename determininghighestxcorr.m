setup;
hyphaelabeling;

extract_oscillation_params;

%% FINDING PIXEL WITH BEST PERIODICITY
ibranch = 52;
shift = 100; %i wanna start into the branch; i suspect edges of branches will be kinda weird and not worth looking at
num_pixels = 853;
pixel_indices = shift:10:num_pixels; %skipping around by 10 px bc Gaussian will muddle across pixels
peak_scores = zeros(1, length(pixel_indices));
%Gaussian params used later
b = gausswin(11);
b = b / sum(b);
a = 1;
for i = 1:length(pixel_indices)
  %getting pixel info
  ipixel = pixel_indices(i);
  thispixelData = branchData{ibranch}(ipixel, :);
  x_coord = thispixelData(1);
  y_coord = thispixelData(2);
  pixelintensity = zeros(1, N_ims+4);
  for i_t = 1:N_ims
      im_smooth = squeeze(im(:,:,i_t));
      pixelintensity(1, i_t+4) = im_smooth(y_coord, x_coord);
  end
  %smooth w/ Gaussian
  pixelsignalfilt = pixelintensity(5:end);
  pixelsignalfilt = pixelsignalfilt - mean(pixelsignalfilt);
  pixelsignalfilt = pixelsignalfilt / std(pixelsignalfilt);
 
  filteredgauss = filtfilt(b, a, pixelsignalfilt);
  filteredgausszerocenter = filteredgauss - mean(filteredgauss);
  [r, lags] = xcorr(filteredgausszerocenter, "coeff"); %coeff normalizes lag 0 to 1.0
  %store the max
  pos_r = r(lags > 10);
  pks = findpeaks(pos_r);
  if ~isempty(pks)
      peak_scores(i) = max(pos_r);
  else
      peak_scores(i) = 0;
  end
end
[sorted_scores, sorted_pixel_indices] = sort(peak_scores, "descend");
%top 10 pixels
topx = 20;
%display
significantpixels = sorted_pixel_indices;
for k = 1:topx
   fprintf('Rank %d Pixel: Index %d (Score = %.3f)\n', ...
       k, sorted_pixel_indices(k), sorted_scores(k));
   significantpixels(k) = sorted_pixel_indices(k);
end
%% GET PARAMS
significantpixelfeatures = zeros(topx, 2); %2D matrix / dict with the pixel and its frequency
for k = 1:topx
   ipixel = significantpixels(k);
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
   filteredgauss = filtfilt(b, a, pixelsignalfilt);
  
   params = extract_oscillation_params(filteredgauss); %using the hilbert transform
  
   %get the frequency of intensity experienced by each pixel
   significantpixelfeatures(k, 1) = ipixel;
   significantpixelfeatures(k, 2) = params.Frequency_Hz;
end

fprintf("frequencies: \n");
for k = 1:topx
   fprintf('Rank %d Pixel %d Frequency = %.3f\n', ...
       k, significantpixelfeatures(k, 1), significantpixelfeatures(k, 2));
end
