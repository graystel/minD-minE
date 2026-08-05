function [filteredsignal, freqs] = singlepixelfftbandpass(pixelintensity, dt, N_ims, showplot, showfreq)
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
    
    filteredsignal = bandpass(pixelsignal, [0.02 0.04], 1/dt);
    
    if showplot == true
        figure;
        plot(1:N_ims, pixelsignal, 1:N_ims, filteredsignal);
        xlabel('frame (# of frames)');
        ylabel('intensity');
        title('Filtered pixel intensity over time');
    end

    if showfreq == true
        freqs = zeros(length(p1), 2);
        for i = 1:length(freqs)
            freqs(i, 1) = f(i);
        end
        for i=1:length(freqs)
            freqs(i, 2) = p1(2);
        end   
    else 
        freqs = null;
    end
end






