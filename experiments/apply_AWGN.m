% function noisy_sig = apply_AWGN(snr_dbs, sig)
%     rms_val = sqrt(mean(abs(sig).^2));
%     sig_W = rms_val^2;  
%     sig_dbW = 10 * log10(sig_W);
%     noise_dbW = sig_dbW - snr_dbs;
%     noise_var = 10.^(noise_dbW / 10);
%     noise_std = sqrt(noise_var);
%     complex_std = noise_std / sqrt(2);
%     noise_samples = complex_std * (randn(size(sig)) + 1j * randn(size(sig)));
%     noisy_sig = sig + noise_samples;
% end

function noisy_sig = apply_AWGN(snr_dbs,sig)
    % Compute the RMS of the signal
    rms = sqrt(mean(abs(sig).^2));
    sig_W = rms.^2;
    
    % Convert signal power to dB
    sig_dbW = 10 .* log10(sig_W / 1);
    
    % Compute noise power in dB
    noise_dbW = sig_dbW - snr_dbs;
    
    % Convert noise power to linear scale
    noise_var = 10.^(noise_dbW / 10);
    
    % Compute noise standard deviation
    noise_std = sqrt(noise_var);
    complex_std = noise_std ./ sqrt(2);
    
    % Generate complex Gaussian noise
    noise_samples = complex_std .* (randn(size(sig)) + 1i .* randn(size(sig)));
    
    % Add noise to the original signal
    noisy_sig = sig + noise_samples;
    
    % Save the noisy signal to a MAT-file
    save('noisy_signal.mat', 'noisy_sig');
end

