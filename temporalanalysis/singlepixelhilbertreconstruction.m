function reconstructed = singlepixelhilbertreconstruction(filteredgauss, params, ipixel, display)
    reconstructed = params.Envelope .* cos(params.Phase);
    freq = params.Frequency_Hz;
    if display == true
        figure;
        yyaxis left;
        plot(1:length(reconstructed), reconstructed, 'o', 1:length(filteredgauss), filteredgauss);
        legend("reconstructed signal", "filtered signal");
        xlabel('frame');
        ylabel('intensity');

        yyaxis right;
        plot(1:length(freq), freq);
        ylabel('freq');
        legend("frequency");
        title("original (filtered) data and reconstruction of pixel " + ipixel + " from Hilbert transform");
    end
end
