function vertical_edges_gt_multipathV2(protocol, CH)
    % Define parameters
    Fs = 80e6; % Sampling frequency in Hz
    Fcut = 9e6; % Cut-off frequency for low-pass filter
    num_realizations = 1000;
    channel_centers = 0;%[-20e6, 0, 20e6]; % Channel center frequencies
    nwin = 256;
    padding_time = 3000;
    overlap = 128;

    % Load waveform
    waveform = load('DATASET1_1_TEST/mat_ax.mat', sprintf('single_waveform_%s', protocol)); 
    waveform = waveform.(sprintf('single_waveform_%s', protocol));
    len = length(waveform);

    % Initialize
    shifted_signals = cell(num_realizations, 1);
    ground_truth_coords = struct('per_channel', {});
    vertical_edge_counts = zeros(3,1);

    % Channel model selection
    switch CH
        case 1
            base_channel = wlanTGacChannel('SampleRate', 20e6, 'ChannelBandwidth', 'CBW20', 'DelayProfile', 'Model-B');
        case 2
            base_channel = wlanTGnChannel('SampleRate', 20e6, 'DelayProfile', 'Model-D');
        case 3
            base_channel = comm.RayleighChannel('SampleRate', 20e6, ...
                                'PathDelays', [0 1e-6 2e-6], ...
                                'AveragePathGains', [0 -3 -6], ...
                                'DopplerSpectrum', doppler('Flat'), ...
                                'MaximumDopplerShift', 30);
        otherwise
            error('Unsupported CH value');
    end

    % Low-pass filter
    lowpassfilter = designfilt('lowpassfir', ...
        'PassbandFrequency', Fcut / (Fs/2), ...
        'StopbandFrequency', (Fcut + 1e6) / (Fs/2), ...
        'PassbandRipple', 1, ...
        'StopbandAttenuation', 45, ...
        'DesignMethod', 'kaiserwin');

    l = grpdelay(lowpassfilter);
    k = l(1);

    fprintf('Starting %d realizations...\n', num_realizations);
    global_start = tic;

    for realization = 1:num_realizations
        realization_start = tic;
        fprintf('Processing realization %d/%d...\n', realization, num_realizations);
        realization_signal = [];
        per_channel_info = struct('packet_starts', {}, 'packet_ends', {}, ...
                                  'T_start', {}, 'T_end', {}, ...
                                  'vertical_edges', {}, 'signal', {});

        last_packet_end = 0;

        for channel_idx = 1
            channel_signal = [];
            packet_start_indices = [];
            packet_end_indices = [];

            % Apply multipath
            reset(base_channel); % Reset channel between uses
            multipath_signal = base_channel(waveform(1:len));
            multipath_signal_resampled = resample(multipath_signal, 4, 1);
            filtered_signal = multipath_signal_resampled;% filter(lowpassfilter, multipath_signal_resampled);

            % Frequency shift
            wshiftDiscrete = 2 * pi * channel_centers(channel_idx) / Fs;
            freq_shifted_signal = multipath_signal_resampled' .* exp(1j * wshiftDiscrete * (0:length(filtered_signal)-1));

            % Apply packet padding
            time_shift = randi([padding_time, 2 * padding_time]);
            freq_shifted_signal = [zeros(1,time_shift), freq_shifted_signal, zeros(1,padding_time)];
            channel_signal = freq_shifted_signal;

            % Store signal and timing info
            start_idx = last_packet_end + time_shift + 1;
            end_idx = start_idx + length(filtered_signal) - 1;
            packet_start_indices = [packet_start_indices, start_idx];
            packet_end_indices = [packet_end_indices, end_idx];

            last_packet_end = last_packet_end + length(freq_shifted_signal);

            Time_start = ((floor((packet_start_indices - nwin) / 128) + 1 + 1.5) * 128) / Fs;
            Time_end = ((floor((packet_end_indices - nwin) / 128) + 1 + 1.5) * 128) / Fs;
            vertical_edges = sort([Time_start, Time_end]);

            vertical_edge_counts(channel_idx) = vertical_edge_counts(channel_idx) + length(vertical_edges);

            % Store ground truth
            per_channel_info(channel_idx).packet_starts = packet_start_indices;
            per_channel_info(channel_idx).packet_ends = packet_end_indices;
            per_channel_info(channel_idx).T_start = Time_start;
            per_channel_info(channel_idx).T_end = Time_end;
            per_channel_info(channel_idx).vertical_edges = vertical_edges;
            per_channel_info(channel_idx).signal = channel_signal;

            realization_signal = [realization_signal, channel_signal];
        end

        shifted_signals{realization} = realization_signal;
        ground_truth_coords(realization).per_channel = per_channel_info;

        elapsed = toc(realization_start);
        fprintf('Finished realization %d/%d in %.2f seconds.\n', realization, num_realizations, elapsed);
    end

    total_time = toc(global_start);
    fprintf('All realizations completed in %.2f seconds (%.2f minutes).\n', total_time, total_time/60);

    % Save output
    save(sprintf('11_%s_shifted_signals_multipath%d.mat', protocol, CH), ...
        'shifted_signals', 'ground_truth_coords', 'vertical_edge_counts', 'protocol');

    fprintf('Saved results to 11_%s_shifted_signals_multipath%d.mat\n', protocol, CH);
end