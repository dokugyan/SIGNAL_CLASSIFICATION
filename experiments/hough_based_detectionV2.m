PROTOCOLS = ["ax","g","n"];
for i = 1:length(PROTOCOLS)
    protocol = PROTOCOLS(i);
    vertical_edges_gt_no_multipath(protocol);  % Generate the ground truth data
    disp("data generated");
    
    load(sprintf('11_%s_shifted_signals_and_ground_truth_random_packets0.mat', protocol));

    % Parameters
    snr_levels = 0:5:30; % SNR levels in dB
    Fs = 80e6; % Sampling frequency in Hz
    Fs_downsampled = 20e6; % Downsampled frequency in Hz
    num_realizations = 1000; % Total number of realizations
    channel_centers = 0e6; % Channel center frequencies
    num_channels = 1; % Number of channels
    freq_band = 10e6; % Frequency band 10 MHz
    downsampling_factor = Fs / Fs_downsampled;
    
    % Preallocate storage
    detected_time_series = cell(num_realizations, length(snr_levels), num_channels);

    % Loop through each realization
    for realization = 1:num_realizations
        per_channel_info = ground_truth_coords(realization).per_channel;

        % Loop through each channel
        for channel_idx = 1:num_channels
            channel_center = channel_centers(channel_idx);
            channel_info = per_channel_info(channel_idx); % Access channel-specific info
            signal_shifted = channel_info.signal;

            % Demodulate the signal (you would apply your own demodulator here)
            shifted_signal = demodulator(signal_shifted, channel_center, Fs);
            
            % Skip if the signal is empty
            if isempty(shifted_signal)
                continue;
            end

            % Loop through each SNR level
            for snr_idx = 1:length(snr_levels)
                snr = snr_levels(snr_idx);
                
                % Apply AWGN noise
                noisy_signal = apply_AWGN(snr, shifted_signal);

                % Perform STFT and crop to the frequency range of the current channel
                [S, F, T] = stft(noisy_signal, Fs, 'Window', hann(256, 'periodic'), 'OverlapLength', 128);
                freq_range = (F >= (channel_center - freq_band)) & (F <= (channel_center + freq_band));
                S_cropped_with_zeros = zeros(size(S));
                S_cropped_with_zeros(freq_range, :) = S(freq_range, :);

                % Create spectrogram and detect edges
                spectrogram_cropped = 20 * log10(abs(S_cropped_with_zeros));
                edges = edge(spectrogram_cropped, 'canny', [0.35, 0.5]);

                % Apply Hough Transform
                [H, theta, rho] = hough(edges);
                peaks = houghpeaks(H, 90, 'threshold', ceil(0.3 * max(H(:))));

                % Extract lines and filter based on theta (horizontal lines)
                lines = houghlines(edges, theta, rho, peaks);
                filtered_lines = lines(abs([lines.theta]) <= 10);

                % Detect peaks in time scores
                logic_matrix = false(size(spectrogram_cropped));
                for k = 1:length(filtered_lines)
                    for x = filtered_lines(k).point1(1):filtered_lines(k).point2(1)
                        logic_matrix(filtered_lines(k).point1(2):filtered_lines(k).point2(2), x) = true;
                    end
                end
                time_scores = sum(logic_matrix, 1);
                detected_peaks = find(time_scores > 0);

                if isempty(detected_peaks)
                    fprintf('No peaks detected in Channel %d for Realization %d, SNR %.1f dB.\n', ...
                        channel_idx, realization, snr);
                    continue;
                end

                % Filter peaks with a minimum gap
                min_gap = 10;
                filtered_peaks = [];
                last_peak = -min_gap;
                for idx = 1:length(detected_peaks)
                    if detected_peaks(idx) >= last_peak + min_gap
                        filtered_peaks = [filtered_peaks, detected_peaks(idx)];
                        last_peak = detected_peaks(idx);
                    end
                end

                % Extract portions of the spectrogram between peaks
                signal_extracted_spectrogram = zeros(size(S));
                for i = 1:length(filtered_peaks) - 1
                    t_start = filtered_peaks(i);
                    t_end = filtered_peaks(i + 1);
                    signal_extracted_spectrogram(:, t_start:t_end) = S_cropped_with_zeros(:, t_start:t_end);
                end

                % ISTFT to retrieve time-domain signal
                signal_extracted = stftRecon(signal_extracted_spectrogram, Fs, hann(256, 'periodic'), 256, 128);

                % Downsample the signal to 20 MHz
                signal_downsampled = downsample(signal_extracted, downsampling_factor);

                % Chop zero-padding by finding non-zero region
                time_end_res = max(find(abs(signal_downsampled) > 0));
                time_start_res = min(find(abs(signal_downsampled) > 0));
                signal_chopped = signal_downsampled(time_start_res:time_end_res);

                % Store the chopped signal for model input
                detected_time_series{realization, snr_idx, channel_idx} = signal_chopped;
                waveform = signal_chopped;

                % Create the directory if it doesn't exist
                output_dir = sprintf("DATASET1_1_TEST/0/802_11%s", protocol);
                if ~isfolder(output_dir)
                    mkdir(output_dir);
                end

                % Save the waveform with a unique filename for each realization and SNR level
                save(sprintf("%s/802.11%s_IQ_frame_realization_%d_snr_%d.mat", output_dir, protocol, realization, snr), "waveform");

                % Free memory after saving
                clear signal_extracted spectrogram_cropped edges H theta rho time_scores filtered_peaks;
            end
        end
    end
end
