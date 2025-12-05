function demodulated_signal = demodulator(Y_bandlimited, channel_centers, Fs)

    % Compute the discrete angular frequency shift
    wshiftDiscrete = 2 * pi * channel_centers / Fs;

    % Create the frequency shift exponential term
    n = 0:(length(Y_bandlimited) - 1); % Time indices
    freq_shift = exp(-1j * wshiftDiscrete * n);

    % Apply the frequency shift to demodulate the signal
    demodulated_signal = Y_bandlimited .* freq_shift;

    % Optional: Return only the real part if necessary
    % demodulated_signal = real(demodulated_signal);
end




% load('ssgt.mat');
% % Parameters
% snr_levels = 10; % SNR levels in dB
% 
% Fs = 80e6; % Sampling frequency in Hz
% num_realizations = 1; % Total number of realizations
% channel_centers = [-20e6, 0, 20e6]; % Channel center frequencies
% num_channels = 3; % Number of channels
% freq_band = 10e6; % Frequency band ±10 MHz
% snr_range = snr_levels - 10 * log10(4); % Adjust SNR range for the system
% disp(snr_range)
% % Preallocate error storage
% vertical_errors = nan(num_realizations, length(snr_range), num_channels);
% 
% % Missing data tracking
% error_log = [];
% missing_count = 0;
% 
% % Missing data tracking summary
% missing_data_summary = zeros(length(snr_range), num_channels); % Rows: SNR, Cols: Channels
% 
% % Loop through each realization
% for realization = 1:num_realizations
%     per_channel_info = ground_truth_coords(realization).per_channel;
% 
%     % Loop through each channel
%     for channel_idx = 1:num_channels
%         channel_center = channel_centers(channel_idx);
%         channel_info = per_channel_info(channel_idx); % Access channel-specific info
% 
% 
% 
%         % Extract the shifted signal for the current channel
%         signal_shifted = channel_info.signal;
% 
%         % Skip if the signal is empty
%         if isempty(signal_shifted)
%             fprintf('No signal for Channel %d in Realization %d. Skipping.\n', channel_idx, realization);
%             continue;
%         end
%         shifted_signal = demodulator(signal_shifted,channel_center(channel_idx),Fs);
% 
% 
%         % Loop through each SNR level
%         for snr_idx = 1:length(snr_range)
%             snr = snr_range(snr_idx);
% 
%             % Add AWGN to the shifted signal
%             noisy_signal = apply_AWGN(signal_shifted, snr);
% 
%             % Perform STFT and crop to the frequency range of the current channel
%             [S, F, T] = stft(noisy_signal, Fs, 'Window', hann(256, 'periodic'), 'OverlapLength', 128);
%             freq_range = F >= (channel_center - freq_band) & F <= (channel_center + freq_band);
%             S_cropped_with_zeros = zeros(size(S));
%             S_cropped_with_zeros(freq_range,:)= S(freq_range,:);
% 
%             spectrogram_cropped = 20 * log10(abs(S_cropped_with_zeros));
% 
%             % Detect edges using the Canny method
% 
% 
%             edges = edge(spectrogram_cropped, 'canny', [0.35, 0.5]);
% 
%             % Apply Hough Transform
%             [H, theta, rho] = hough(edges);
%             peaks = houghpeaks(H, 90, 'threshold', ceil(0.3 * max(H(:))));
%             lines = houghlines(edges, theta, rho, peaks);
%             % Apply Hough Transform
%             % [H, theta, rho] = hough(edges);
%             % 
%             % % Restrict theta to the range of interest (-5 to 5 degrees)
%             % theta_range = (theta >= -10 & theta <= 10);
%             % H_filtered = H(:, theta_range);  % Filter H for the restricted theta range
%             % theta_filtered = theta(theta_range);  % Filter theta values accordingly
%             % 
%             % % Detect peaks in the filtered Hough transform
%             % peaks = houghpeaks(H_filtered, 90, 'threshold', ceil(0.3 * max(H_filtered(:))));
%             % 
%             % % Map the peaks back to the original theta and rho
%             % lines = houghlines(edges, theta_filtered, rho, peaks);
%             % 
% 
%             % Filter lines based on theta value (horizontal lines)
%             filtered_lines = lines(abs([lines.theta]) <= 10); %10
% 
%             % Initialize logical matrix for time scoring
%             logic_matrix = zeros(size(spectrogram_cropped));
% 
%             % Create logical matrix for the detected lines
%             for k = 1:length(filtered_lines)
%                 for x = filtered_lines(k).point1(1):filtered_lines(k).point2(1)
%                     logic_matrix(filtered_lines(k).point1(2):filtered_lines(k).point2(2), x) = 1;
%                 end
%             end
% 
%             % Calculate time scores by summing across the frequency axis (columns)
%             time_scores = sum(logic_matrix, 1);
% 
%             % Identify detected peaks in time scores
%             detected_peaks = find(time_scores > 0);
%             if isempty(detected_peaks)
%                 fprintf('No peaks detected in Channel %d for Realization %d, SNR %.1f dB.\n', ...
%                     channel_idx, realization, snr);
%                 error_log = [error_log; realization, channel_idx, snr_idx];
%                 missing_count = missing_count + 1;
% 
%                 % Increment the missing data summary for this SNR and channel
%                 missing_data_summary(snr_idx, channel_idx) = missing_data_summary(snr_idx, channel_idx) + 1;
%                 continue;
%             end
% 
%             % Ensure there is a gap of at least 10 indices between detected peaks
%             min_gap = 10;
%             filtered_peaks = [];
%             last_peak = -min_gap;
% 
%             for idx = 1:length(detected_peaks)
%                 if detected_peaks(idx) >= last_peak + min_gap
%                     filtered_peaks = [filtered_peaks, detected_peaks(idx)];
%                     last_peak = detected_peaks(idx);
%                 end
%             end
% 
%             % Convert filtered peaks to time values
%             detected_edges = T(filtered_peaks)';
% 
%             %% Reconstructed signal 
%             % Get the start and end times
%             % Extract portions of the spectrogram between peaks
%             signal_extracted_spectrogram = zeros(size(S));
%             for i = 1:length(filtered_peaks) - 1
%                 t_start = filtered_peaks(i);
%                 t_end = filtered_peaks(i + 1);
%                 signal_extracted_spectrogram(:, t_start:t_end) = S_cropped_with_zeros(:, t_start:t_end);
%             end
% 
%             S_extracted = signal_extracted_spectrogram;
% 
%             nwin = 256;
%             win = hann(nwin, 'periodic');
%             waveform = stftRecon(S_extracted, Fs,win,nwin,nwin/2);
%             waveform = downsample(waveform,4);
%             %% Reverse Filtering 
% 
%             % Chop zero-padding by finding non-zero region
%             non_zero_indices = find(waveform ~= 0);
%             if ~isempty(non_zero_indices)
%                 start_idx = non_zero_indices(2);
%                 end_idx = non_zero_indices(end);
%                 waveform = waveform(start_idx:end_idx);
%             else
%                 waveform = []; % Handle all-zero case
%             end
% 
%             % Store the chopped signal for model input
%             % detected_time_series{realization, snr_idx, channel_idx} = signal_chopped;
% 
%             % Display results
%             fprintf('Realization %d, Channel %.1f MHz, SNR %.1f dB:\n', ...
%                 realization, channel_center / 1e6, snr);
%             fprintf('Chopped Signal Length: %d\n', length(waveform));
% 
% 
% 
%             figure
%             plot(abs(waveform),'r');
% 
% 
% 
%             save("DATASET1_1_TEST/802_11ax/802.11ax_IQ_frame_1.mat","waveform");
% 
% 
% 
%             % Clip ground truth edges that exceed T(end)
%             ground_truth_edges = channel_info.vertical_edges;
%             if length(detected_edges) < length(ground_truth_edges)
%                 detected_edges = [detected_edges, zeros(1, length(ground_truth_edges) - length(detected_edges))];
%             elseif length(detected_edges) > length(ground_truth_edges)
%                 detected_edges = detected_edges(1:length(ground_truth_edges));
%             end
% 
%             %gt_edges = ground_truth_edges(ground_truth_edges <= T(end));
%             edge_errors = detected_edges - ground_truth_edges;
%             % Compute the error between detected and ground truth edges
%             vertical_edges_errors(realization, snr_idx, 1:length(edge_errors)) = edge_errors;
% 
%             fprintf('Realization %d, Channel %.1f MHz, SNR %.1f dB:\n', ...
%                 realization, channel_centers(channel_idx) / 1e6, snr);
%             fprintf('Detected edges: %s\n', mat2str(detected_edges));
%             fprintf('Ground truth edges: %s\n', mat2str(ground_truth_edges));
%             fprintf('Edge errors: %s\n\n', mat2str(edge_errors));
%         end
%     end
% end