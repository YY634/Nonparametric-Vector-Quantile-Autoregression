import os
import mne
import glob
import numpy as np
from pathlib import Path


def process_eeg_to_npy(input_dir=".", pattern="*.set", output_dir="./npy"):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    set_files = sorted(
        glob.glob(os.path.join(input_dir, "*", pattern))
    )
    if not set_files:
        print("Saved 0 files processed (no .set found).")
        return

    saved = 0
    for f in set_files:
        base = f.stem
        out_path = output_dir / f"{base}-ts.npy"
        print(f"Processing {f}...")
        raw = mne.io.read_raw_eeglab(str(f), preload=True, verbose=False)
        data, times = raw.get_data(return_times=True)
        np.save(out_path, data)
        saved += 1

    print(f"Saved {saved} files processed")

if __name__ == "__main__":
    process_eeg_to_npy(
        input_dir=r"C:/Users/Desktop/NVQAR/derivatives",
        pattern="*.set",
        output_dir=r"C:/Users/Desktop/NVQAR/derivatives/npy"
    )

