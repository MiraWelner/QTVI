#!/usr/bin/env python3
import os


def pad(val, length):
    """Ensure value is exactly 'length' bytes for the EDF header."""
    return str(val)[:length].ljust(length).encode("ascii")


def process_file(path):
    with open(path, "rb") as f:
        # 1. Main Header is exactly 256 bytes
        main_header = f.read(256)
        ns = int(main_header[252:256].decode().strip())  # Number of signals

        # 2. Read Signal Headers Column by Column
        labels = f.read(ns * 16)
        transducers = f.read(ns * 80)
        phys_dims = f.read(ns * 8)

        # We need the numbers from these columns to calculate new scaling
        phys_mins = [float(f.read(8).decode().strip()) for _ in range(ns)]
        phys_maxs = [float(f.read(8).decode().strip()) for _ in range(ns)]
        dig_mins = [int(f.read(8).decode().strip()) for _ in range(ns)]
        dig_maxs = [int(f.read(8).decode().strip()) for _ in range(ns)]

        prefilters = f.read(ns * 80)
        nr_samples = f.read(ns * 8)
        reserved = f.read(ns * 32)

        # 3. Read all signal data
        data = f.read()

    # 4. Calculate new Physical limits so the signal doesn't get "flat"
    # We maintain the original gain: (phys_range / dig_range)
    new_p_mins = []
    new_p_maxs = []
    for i in range(ns):
        gain = (phys_maxs[i] - phys_mins[i]) / (dig_maxs[i] - dig_mins[i])
        new_p_mins.append(phys_mins[i] + (-32768 - dig_mins[i]) * gain)
        new_p_maxs.append(phys_maxs[i] + (32767 - dig_maxs[i]) * gain)

    # 5. Write back the file with the expanded bounds
    with open(path, "wb") as f:
        f.write(main_header)
        f.write(labels)
        f.write(transducers)
        f.write(phys_dims)
        for val in new_p_mins:
            f.write(pad(val, 8))
        for val in new_p_maxs:
            f.write(pad(val, 8))
        for _ in range(ns):
            f.write(pad("-32768", 8))  # New Digital Min
        for _ in range(ns):
            f.write(pad("32767", 8))  # New Digital Max
        f.write(prefilters)
        f.write(nr_samples)
        f.write(reserved)
        f.write(data)


def main():
    import os
    base = r"D:\USERS\MiraWelner\QTVI\data\mesa_raw_files\bl_mi"
    ids = [
        folder.removesuffix("_EDF")
        for folder in os.listdir(base)
        if os.path.isdir(os.path.join(base, folder))
    ]
    for rid in ids:
        # e.g. ...\nondiabetic\3010970_20120104_EDF\3010970_20120104.edf
        edf_path = os.path.join(base, f"{rid}_EDF", f"{rid}.edf")
        if os.path.exists(edf_path):
            print(f"Unbounding: {edf_path}")
            try:
                process_file(edf_path)
            except Exception as e:
                print(f"  Error: {e}")
        else:
            print(f"  Not found: {edf_path}")

if __name__ == "__main__":
    main()
