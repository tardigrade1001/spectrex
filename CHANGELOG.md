# Changelog

All notable releases of spectrex are documented here. Newest entries appear first.

---

## v1.1.0 — Desktop workflow update

- Added a native Windows interface: drag in a folder or individual files, or choose them with the new buttons. The app shows per-file progress, a clear completion summary, and an expandable error log when anything fails.
- Corrected FDS timestamp parsing and excluded four footer values that could otherwise appear as false spectrum data beyond the requested end wavelength.
- Replaced the old icon with a spectrum-and-cuvette mark.

---

## v1.0.0 — First public release

First public release of spectrex.

### What it does

Converts Hitachi spectrophotometer binary files to plain CSV and quick-look PNG plots on any PC, without needing the original 1995-era software.

Supports:
- `.UDS` files from the Hitachi U-2900 UV-Vis Spectrophotometer (UV Solutions)
- `.FDS` files from the Hitachi F-4600 FL Spectrophotometer (FL Solutions)

### For non-Python users

Download `spectrex.exe` below (64 MB, single file, no installation).

Two ways to run it:
1. Double-click the exe, choose a folder, and press **Convert files**. The window displays live progress and a clear result summary.
2. Drag a folder *onto* `spectrex.exe` in File Explorer. The exe processes that folder and its subfolders.

A `spectrex.log` summary is written at the end of every run.

CSV and PNG files are written next to each original binary, and `spectrex.log` contains the complete run summary.

### For Python users

Clone the repo and run `python spectrex.py` against the same folder structure. The source has no compiled dependencies beyond `matplotlib`.

### What's in the CSV

Each converted file gets a CSV containing the full spectrum (wavelength, absorbance or fluorescence intensity) with the instrument metadata in the header comments:

- Sample name, timestamp, instrument model, serial number
- Sampling step, wavelength range
- UDS: slit width, scan speed, path length, lamp change wavelength, baseline correction, response setting (full parameter set)
- FDS: excitation wavelength

### Verification

Output was verified byte-for-byte against the original program's own TXT exports. Maximum absolute difference: 0.0005 for absorbance (matches the TXT's three-decimal display precision) and ~0.5 for fluorescence intensity (matches its four-significant-figure precision).

### Known limitation

A subset of FDS acquisition parameters (scan speed, EX/EM slit widths, PMT voltage, response time, delay) are stored as setting codes inside the binary, with the codes-to-values lookup compiled into the FL Solutions DLLs. These fields appear as a `# Note:` line in the FDS CSV header rather than as extracted values. The spectrum data itself is decoded correctly and is not affected. See the README for the full story.

### Antivirus note

PyInstaller-bundled exes occasionally trigger Windows Defender false positives. If your AV flags `spectrex.exe`, the source in this repo is what was built. Feel free to verify by building locally with `pyinstaller --onefile --windowed --icon=spectrex.ico spectrex.py`.
