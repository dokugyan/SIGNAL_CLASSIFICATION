PROTOCOLS = ["ax","g","n"];
for i = 1:length(PROTOCOLS)
    CH = 0;
    protocol = PROTOCOLS(i);
    % vertical_edges_gt_multipathV2(protocol,CH);
    vertical_edges_gt_no_multipath(protocol)
    disp("data generated");
    % load(sprintf('11_%s_shifted_signals_mult ipath%d.mat', protocol,CH));
    load(sprintf('11_%s_shifted_signals_and_ground_truth_random_packets0.mat', protocol));
    % Parameters
    snr_levels = 0:5:30;
    Fs = 80e6;
    Fs_downsampled = 20e6;
    num_realizations = 1000;
    channel_centers = 0;
    num_channels = 1;
    freq_band = 10e6;
    downsampling_factor = Fs / Fs_downsampled;
    start_threshold = 0.015;
    end_threshold = -0.035;
    min_gap = 50;

    detected_time_series = cell(num_realizations, length(snr_levels), num_channels);

    for realization = 1:num_realizations
        per_channel_info = ground_truth_coords(realization).per_channel;
        for channel_idx = 1:num_channels
            channel_center = channel_centers(channel_idx);
            channel_info = per_channel_info(channel_idx);
            signal_shifted = channel_info.signal;
            shifted_signal = demodulator(signal_shifted, channel_center, Fs);

            if isempty(shifted_signal)
                continue;
            end

            for snr_idx = 1:length(snr_levels)
                snr = snr_levels(snr_idx);
                noisy_signal = apply_AWGN(snr, shifted_signal);
                [S, F, T] = stft(noisy_signal, Fs, 'Window', hann(256, 'periodic'), 'OverlapLength', 128);

                freq_range = F >= (channel_center - freq_band) & F <= (channel_center + freq_band);
                S_cropped_with_zeros = zeros(size(S));
                S_cropped_with_zeros(freq_range, :) = S(freq_range, :);

                subchannel_power = sum(abs(S_cropped_with_zeros).^2, 1);
                subchannel_power = subchannel_power / max(subchannel_power);
                smooth_power = smooth(subchannel_power, 30);
                grad_mag = diff(smooth_power');

                start_edges = grad_mag > start_threshold;
                end_edges = grad_mag < end_threshold;
                detected_peaks = find(start_edges | end_edges);
                filtered_peaks = [];
                last_peak = -min_gap;

                for idx = 1:length(detected_peaks)
                    if detected_peaks(idx) >= last_peak + min_gap
                        filtered_peaks = [filtered_peaks, detected_peaks(idx)];
                        last_peak = detected_peaks(idx);
                    end
                end

                signal_extracted_spectrogram = zeros(size(S));
                for j = 1:length(filtered_peaks) - 1
                    t_start = filtered_peaks(j);
                    t_end = filtered_peaks(j + 1);
                    signal_extracted_spectrogram(:, t_start:t_end) = S_cropped_with_zeros(:, t_start:t_end);
                end

                signal_extracted = stftRecon(signal_extracted_spectrogram, Fs, hann(256, 'periodic'), 256, 128);
                signal_downsampled = downsample(signal_extracted, downsampling_factor);
                time_end_res = max(find(abs(signal_downsampled) > 0));
                time_start_res = min(find(abs(signal_downsampled) > 0));
                signal_chopped = signal_downsampled(time_start_res:time_end_res);

                detected_time_series{realization, snr_idx, channel_idx} = signal_chopped;
                waveform = signal_chopped;
                base_dir = "DATASET1_1_TEST";
                subfolder = sprintf("%s/%d", base_dir, CH);
                if ~isfolder(subfolder)
                    mkdir(subfolder);
                end

                % Create the protocol-specific folder inside the numbered subfolder
                output_dir = sprintf("%s/802_11%s", subfolder, protocol);
                if ~isfolder(output_dir)
                    mkdir(output_dir);
                end

                % Save the waveform with a unique filename
                save(sprintf("%s/802.11%s_IQ_frame_realization_%d_snr_%d.mat", ...
                    output_dir, protocol, realization, snr), "waveform");

                % Clear variables to free memory
                clear signal_extracted spectrogram_cropped subchannel_power grad_mag start_edges end_edges;
            end
        end
    end

end
