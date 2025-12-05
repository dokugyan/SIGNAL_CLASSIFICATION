function vertical_edges_gt_no_multipath(protocol)

% Define parameters
Fs = 80e6; % Sampling frequency in Hz
Fcut = 9e6; % Cut-off frequency for low-pass filter
%T = 5e-3; % Duration in seconds
num_realizations = 1000; % Number of realizations
channel_centers = 0e6; % Channel center frequencies
nwin = 256;
padding_time = 3000; % Minimum gap between packets in samples
overlap = 128;
% Initialize arrays to store shifted signals and ground truth coordinates
shifted_signals = cell(num_realizations, 1);
ground_truth_coords = struct('per_channel', {});
vertical_edge_counts = zeros(3,1);
waveform = load('DATASET1_1_TEST/mc3_6/mat_ax.mat', sprintf('single_waveform_%s', protocol));
waveform = waveform.(sprintf('single_waveform_%s', protocol));
len = length(waveform);
Y1 = waveform()';

%{
% Generate 11ax waveforms (mocked here, you would replace with actual waveform generation)
waveforms = cell(1, 3);
for i = 1:3
    waveforms{i} = randn(1, 6800); % Replace with actual 802.11ax waveform generation
end
%}
% Define multipath channel models for three transmitters
tgnChan = wlanTGnChannel('SampleRate', 20e6, 'DelayProfile', 'Model-D');
tgxChan = wlanTGacChannel('SampleRate', 20e6, 'ChannelBandwidth', 'CBW20', 'DelayProfile', 'Model-B');
rayleighChan = comm.RayleighChannel('SampleRate', 20e6, ...
                                    'PathDelays', [0 1e-6 2e-6], ...
                                    'AveragePathGains', [0 -3 -6], ...
                                    'DopplerSpectrum', doppler('Flat'), ...
                                    'MaximumDopplerShift', 30);

channels = {tgxChan,tgnChan,rayleighChan};

lowpassfilter = designfilt('lowpassfir', ...
    'PassbandFrequency', Fcut / (Fs/2), ...
    'StopbandFrequency', (Fcut + 0.2e6) / (Fs/2), ...
    'PassbandRipple', 1, ...
    'StopbandAttenuation', 45, ...
    'DesignMethod', 'kaiserwin');
l = grpdelay(lowpassfilter);
k = l(1);
% Process each realization
for realization = 1:num_realizations
    realization_signal = [];
    per_channel_info = struct('packet_starts', {}, 'packet_ends', {}, ...
                               'T_start', {}, 'T_end', {}, ...
                               'vertical_edges', {}, 'signal', {});
    
    % Loop over each channel/transmitter
    for channel_idx = 1
        channel_signal = [];
        packet_start_indices = [];
        packet_end_indices = [];
        Time_start = [];
        Time_end = [];
        last_packet_end = 0;

        % Random number of packets1 to 3) for this channel
        % num_packets = randi([1,3]);
        num_packets = 1;
        % Generate and process each packet
        for packet_idx = 1:num_packets
            % Apply multipath effects
            %multipath_signal = channels{channel_idx}(waveform(1:len));
            multipath_signal = waveform(1:len).';
            multipath_signal_resampled = resample(multipath_signal, 4, 1); % Upsample by 4
            
            filtered_signal = filter(lowpassfilter, multipath_signal_resampled);

            % Resample to 80 MHz
            %multipath_signal_resampled = multipath_signal_resampled.';
            % Apply frequency shift
            wshiftDiscrete = 2 * pi * channel_centers(channel_idx) / Fs;
            freq_shifted_signal = multipath_signal_resampled .* exp(1j * wshiftDiscrete * (0:length(multipath_signal_resampled) - 1));
            
            % Ensure a gap between packets
            time_shift = randi([padding_time, 2 * padding_time]);
            freq_shifted_signal = [zeros(1,time_shift), freq_shifted_signal, zeros(1,padding_time)];
            
            % Update channel signal
            channel_signal = [channel_signal, freq_shifted_signal];
            
            % Record packet start and end indices
            start_idx = last_packet_end + time_shift + 1;
            end_idx = start_idx + length(filtered_signal) - 1;
            packet_start_indices = [packet_start_indices, start_idx];
            packet_end_indices = [packet_end_indices, end_idx];
            
            % Update last packet end position
            last_packet_end = last_packet_end + length(freq_shifted_signal);
        end
        
        % Calculate time indices and vertical edges
        Time_start = ((floor((packet_start_indices - nwin) / 128) + 1 + 1.5) * 128) / Fs;
        Time_end = ((floor((packet_end_indices - nwin) / 128) + 1 + 1.5) * 128) / Fs;
        vertical_edges = sort([Time_start, Time_end]);
        
        vertical_edge_counts(channel_idx) = vertical_edge_counts(channel_idx) + length(vertical_edges);


        % Save per-channel ground truth
        per_channel_info(channel_idx).packet_starts = packet_start_indices;
        per_channel_info(channel_idx).packet_ends = packet_end_indices;
        per_channel_info(channel_idx).T_start = Time_start;
        per_channel_info(channel_idx).T_end = Time_end;
        per_channel_info(channel_idx).vertical_edges = vertical_edges;
        per_channel_info(channel_idx).signal = channel_signal;
        
        % Add channel signal to the realization signal
        if isempty(realization_signal)
            realization_signal = channel_signal;
        else
            realization_signal = [realization_signal, channel_signal];
        end
    end
    
    % Save the realization signal and per-channel ground truth
    shifted_signals{realization} = realization_signal;
    ground_truth_coords(realization).per_channel = per_channel_info;
end
% disp('Vertical edge counts per channel:');
% disp(vertical_edge_counts);
% % Plot STFT to verify signals in different frequency bands
% figure;
% [S, F, T] = stft(realization_signal, Fs, 'Window', hamming(nwin), 'OverlapLength', 128);
% spectrogram_dB = 20 * log10(abs(S));
% imagesc(T, F, spectrogram_dB);
% %stft(real(realization_signal), Fs, 'Window', hamming(nwin), 'OverlapLength', nwin / 2, 'FFTLength', nwin);
% title('STFT of Combined Signal');
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');

% Save the shifted signals and ground truth
save(sprintf('11_%s_shifted_signals_and_ground_truth_random_packets0.mat', protocol), 'shifted_signals', 'ground_truth_coords', 'vertical_edge_counts', 'protocol');

% clc;
% clear;
end
