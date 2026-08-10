function reconstructed = singlepixelhilbertreconstruction(params, display)
    reconstructed = params.Envelope .* cos(params.Phase);
    if display == true
        plot(1:length(reconstructed), reconstructed, 'o', 1:N_ims, filteredgauss, 1:N_ims, pixelsignal, 1:N_ims, params.Frequency);
        legend("reconstructed signal", "filtered signal", "original signal");
        xlabel('frame');
        ylabel('intensity');
        title('original (filtered) data and reconstruction from Hilbert transform');
    end
end
