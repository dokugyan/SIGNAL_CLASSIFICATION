import scipy.io as sio
import time
import numpy as np
import torch
import os
from glob import glob
import argparse
import re  # For filename parsing
import sys

sys.path.append('/home/jay/t-prime-STABLE-2')
from baseline_models.model_MCFormer import MCformer

#TRANS_PATH = '../TPrime_transformer/model_cp'
PROTOCOLS = ['802_11ax', '802_11b', '802_11n', '802_11g']
CHANNELS = ['None']


def chan2sequence(obs):
    seq = np.empty(obs.size)
    seq[0::2] = obs[0]
    seq[1::2] = obs[1]
    return seq


def process_all_packets(model, class_map, seq_len, sli_len, channel, protocol, test_data_path, output_base_path):
    input_path = os.path.join(test_data_path, protocol)
    mat_list = sorted(glob(os.path.join(input_path, '*.mat')))
    print(f"\n✅ Found {len(mat_list)} .mat files for protocol: {protocol}")

    output_dir = os.path.join(output_base_path or test_data_path, protocol)
    os.makedirs(output_dir, exist_ok=True)

    for idx, signal_path in enumerate(mat_list):
        print(f"\n🔍 Processing file {idx + 1}/{len(mat_list)}: {signal_path}")
        sig = sio.loadmat(signal_path)

        if channel == 'None':
            sig = sig['waveform']
        len_sig = sig.shape[0]
        if len_sig == 0:
            continue

        noisy_sig = sig  # already preprocessed (no noise added here)
        if noisy_sig.shape == 0:
            print(f"⚠️ Skipping file due to empty noisy signal: {signal_path}")
            continue

        obs = np.stack((noisy_sig.real, noisy_sig.imag))
        obs = np.squeeze(obs, axis=2)
                
        idxs = list(range(sli_len, len_sig, sli_len))
        obs = np.split(obs, idxs, axis=1)[:-1]
        
        X = torch.tensor(np.asarray(obs), dtype=torch.float32).to(device)

        if X.numel() == 0 or X.shape[0] == 0:
            print(f"⚠️ Skipping file due to empty input tensor: {signal_path}")
            continue

        y = torch.full((len(X),), class_map[protocol], dtype=torch.long).to(device)

        row_percentages_list = []

        with torch.no_grad():
            start_time = time.time()
            pred = model(X)
            end_time = time.time()

            probabilities = torch.softmax(pred, dim=1)
            percentages = probabilities * 100

            # Save row percentages as doubles (float32 or float64)
            row_percentages_list = percentages.cpu().numpy().astype(np.float64)  # Convert to float64 (double in MATLAB)

            print("📊 Class probabilities (%):")
            for i, row in enumerate(percentages):
                row_percentages = ', '.join([f"{p.item():.2f}%" for p in row])
                print(f"   - Segment {i + 1}: {row_percentages}")

        # 💾 Save output using real filename info
        input_basename = os.path.basename(signal_path)
        base_filename = os.path.splitext(input_basename)[0]

        match = re.search(r"realization_(\d+)_snr_(\d+)", base_filename)
        if match:
            realization = match.group(1)
            snr = match.group(2)
            output_filename = f"probabilities_realization_{realization}_snr_{snr}.mat"
        else:
            output_filename = f"{base_filename}_probabilities.mat"

        output_path = os.path.join(output_dir, output_filename)
        sio.savemat(output_path, {'row_percentages': row_percentages_list})
        print(f"💾 Saved class probabilities to: {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--protocols', nargs='+', default=['802_11ax'], choices=PROTOCOLS)
    parser.add_argument('--data_path', default=None, help='Path to test dataset folder')
    parser.add_argument('--output_path', default=None, help='Optional path to save output .mat files')
    parser.add_argument('--device', default='cpu', help='cpu or cuda')
    args = parser.parse_args()

    class_map = dict(zip(PROTOCOLS, range(len(PROTOCOLS))))
    device = torch.device(args.device)

    print(f"📦 Loading model on device: {device}")
    model = MCformer()
    model.load_state_dict(torch.load(f"/home/jay/t-prime-STABLE-2/baseline_models/mcf_weights/model.MCformer.random.range.pt", map_location=device)['model_state_dict'])
    model.to(device)
    model.eval()

    for protocol in args.protocols:
        for channel in CHANNELS:
            process_all_packets(model, class_map, seq_len=1, sli_len=128, channel=channel,
                                protocol=protocol, test_data_path=args.data_path,
                                output_base_path=args.output_path)

    print("\n🚀 All protocols processed.")

