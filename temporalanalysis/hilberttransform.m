function params = hilberttransform(signal, fps)
    arguments
        signal
        fps = 2
    end
    
    sig_centered = signal - mean(signal);
    
    z = hilbert(sig_centered); %hilbert applies the hilbert transform, which returns the amplitude and instantaneous phase?
    envelope = abs(z); %amplitude
    inst_phase = unwrap(angle(z)); %phase
    
    [r, lags] = xcorr(sig_centered, 'coeff'); %finds the autocorrelation for the period
    pos_r = r(lags > 0);
    pos_lags = lags(lags > 0);
    [pks, pk_lags] = findpeaks(pos_r, pos_lags, 'MinPeakDistance', 3);
    
    if ~isempty(pks)
        T_frames = pk_lags(1);
    else
        T_frames = NaN;
    end

    f = (fps / (2 * pi)) * diff(inst_phase);
    
    params.Period_frames = T_frames;
    params.Period_seconds = T_frames / fps;
    params.Frequency_Hz = f; %fps / T_frames;
    params.Mean_Amplitude = mean(envelope);
    params.Max_Amplitude = max(envelope);
    params.Envelope = envelope;
    params.Phase = inst_phase;
end