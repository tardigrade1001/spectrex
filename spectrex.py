"""Convert Hitachi UV-Vis (.UDS) and Fluorescence (.FDS) binary files to CSV+PNG.

Usage:  python spectrex.py            # open the folder-picker window
        python spectrex.py <folder>   # walk the given folder recursively

UDS (U-2900 spectrophotometer):
  Magic IIHIITAG. Stores transmittance in scan order (high -> low wavelength).
  We convert T -> A = -log10(T) and emit ascending-wavelength CSV.

FDS (F-4600 fluorescence spectrometer):
  Magic IIHIDTAG. Each data record is 5 oversampled doubles (0.2 nm stride).
  We output the full 0.2 nm internal resolution.

Run output:
  OK   <relpath>: <kind> <N> pts <start>-<end> nm
  FAIL <relpath> [<stage>]: <error>     # diagnostic info follows on next lines

Errors are also appended to spectrex.log alongside the script for later review.
"""
import struct, os, sys, glob, math, traceback, datetime, threading, queue
import matplotlib.pyplot as plt

# ---------------- diagnostics ----------------

class ParseError(RuntimeError):
    """Raised by parsers with a human-readable explanation of what went wrong."""

LOG_FILE = None  # set in main()
LOG_SINK = print

def log(msg):
    LOG_SINK(msg)
    if LOG_FILE:
        try:
            with open(LOG_FILE, "a", encoding="utf-8") as f:
                f.write(msg.rstrip() + "\n")
        except OSError:
            pass

# ---------------- parsing ----------------

def _cstr(buf, off, field_name="string"):
    if off >= len(buf):
        raise ParseError(f"reached EOF while reading {field_name} (offset {off}, file size {len(buf)})")
    try:
        end = buf.index(b"\x00", off)
    except ValueError:
        raise ParseError(f"no null terminator found for {field_name} starting at offset {off}")
    return buf[off:end].decode("latin-1", errors="replace"), end + 1

def parse_uds(buf):
    if len(buf) < 64:
        raise ParseError(f"file too small to be a UDS ({len(buf)} bytes)")
    off = 16
    try:
        sample, off     = _cstr(buf, off, "sample name")
        timestamp, off  = _cstr(buf, off, "timestamp")
        instrument, off = _cstr(buf, off, "instrument model")
        serial, off     = _cstr(buf, off, "serial number")
        version, off    = _cstr(buf, off, "ROM version")
    except ParseError as e:
        raise ParseError(f"UDS header strings: {e}")

    valid_steps = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}
    search_end = min(len(buf), off + 200)
    data_off = None
    rejection_reason = "no [step, start_wavelength] pair satisfied all constraints"
    strings_end = off  # remember where the c-strings stopped
    for i in range(off, search_end - 16):
        step, start_wl = struct.unpack_from("<dd", buf, i)
        if step not in valid_steps:
            continue
        if not (180 <= start_wl <= 1500):
            continue
        if abs(start_wl / step - round(start_wl / step)) > 1e-6:
            continue
        v0 = struct.unpack_from("<d", buf, i + 16)[0]
        if not (0 < v0 < 1.5):
            rejection_reason = (f"found candidate step={step} start_wl={start_wl} at offset {i}, "
                                f"but first data double {v0!r} is not a plausible transmittance (0..1.5)")
            continue
        data_off = i + 16
        step_used, start_wl_used = step, start_wl
        break
    if data_off is None:
        raise ParseError(f"UDS data block not found in bytes {off}..{search_end}. "
                         f"Valid steps tried: {sorted(valid_steps)}. Last rejection: {rejection_reason}")
    step, start_wl = step_used, start_wl_used

    # ---- extra header params (best-effort) ----
    # Layout between end-of-strings and the [step, start_wl] doubles:
    #   slit_width (double), zeros, "None"/baseline_correction string,
    #   "Medium"/response string, lamp_change_wavelength (double).
    lamp_change_off = data_off - 24            # 8 bytes before step
    lamp_change_wl = struct.unpack_from("<d", buf, lamp_change_off)[0]
    # Slit width: first plausible (0.05..10) double in [strings_end .. lamp_change_off).
    slit_width = None
    for j in range(strings_end, lamp_change_off - 8, 1):
        v = struct.unpack_from("<d", buf, j)[0]
        if 0.05 <= v <= 10 and abs(v * 10 - round(v * 10)) < 1e-6:
            slit_width = v
            slit_end = j + 8
            break
    # baseline_correction, response: two c-strings sitting between slit and lamp_change,
    # with zero-padding before each.
    baseline_correction = response_setting = None
    if slit_width is not None:
        def _skip_zeros(p, limit):
            while p < limit and buf[p] == 0:
                p += 1
            return p
        try:
            p = _skip_zeros(slit_end, lamp_change_off)
            baseline_correction, p = _cstr(buf, p, "baseline correction")
            p = _skip_zeros(p, lamp_change_off)
            response_setting, _ = _cstr(buf, p, "response setting")
        except ParseError:
            pass

    n = 0
    while data_off + (n + 1) * 8 <= len(buf):
        v = struct.unpack_from("<d", buf, data_off + n * 8)[0]
        if v != v or abs(v) > 5:
            break
        n += 1
    if n == 0:
        raise ParseError(f"UDS data array is empty at offset {data_off}")
    if n < 5:
        raise ParseError(f"UDS data array suspiciously short ({n} points) at offset {data_off}")

    transmittance = struct.unpack_from("<" + "d" * n, buf, data_off)
    if any(t <= 0 for t in transmittance):
        zero_count = sum(1 for t in transmittance if t <= 0)
        log(f"     note: {zero_count}/{n} transmittance values are <=0 (very absorbing sample); "
            f"absorbance set to NaN at those points")
    absorbance = [-math.log10(t) if t > 0 else float("nan") for t in transmittance]
    absorbance.reverse()
    end_wl = start_wl - (n - 1) * step
    wavelengths = [end_wl + i * step for i in range(n)]

    # ---- footer params (best-effort): scan_speed at +16, end_wl at +32, path_length at +40 ----
    foot = data_off + n * 8
    scan_speed = path_length = end_wl_footer = None
    if foot + 48 <= len(buf):
        try:
            _, _, ss, _, ew, pl = struct.unpack_from("<6d", buf, foot)
            if 1 <= ss <= 5000:   scan_speed = ss
            if 180 <= ew <= 1500: end_wl_footer = ew
            if 0.1 <= pl <= 100:  path_length = pl
        except struct.error:
            pass

    return dict(kind="UDS", sample=sample, timestamp=timestamp,
                instrument=instrument, serial=serial, version=version,
                y_label="absorbance", wavelengths=wavelengths, values=absorbance,
                step_nm=step, slit_width_nm=slit_width,
                lamp_change_wl_nm=lamp_change_wl,
                baseline_correction=baseline_correction,
                response_setting=response_setting,
                scan_speed_nm_per_min=scan_speed,
                path_length_mm=path_length)

def parse_fds(buf):
    if len(buf) < 64:
        raise ParseError(f"file too small to be an FDS ({len(buf)} bytes)")
    try:
        anchor = buf.index(b"F-4600")
    except ValueError:
        raise ParseError("instrument anchor 'F-4600' not found in file "
                         "(spectrex only knows the F-4600 fluorescence format)")

    off = 16
    try:
        sample, off    = _cstr(buf, off, "sample name")
        operator, off  = _cstr(buf, off, "operator")
        comment, off   = _cstr(buf, off, "comment")
        timestamp, off = _cstr(buf, off, "timestamp")
    except ParseError as e:
        raise ParseError(f"FDS header strings: {e}")

    off = anchor
    try:
        instrument, off = _cstr(buf, off, "instrument model")
        rom, off        = _cstr(buf, off, "ROM version")
        serial, off     = _cstr(buf, off, "serial number")
    except ParseError as e:
        raise ParseError(f"FDS instrument-block strings: {e}")

    if off + 16 > len(buf):
        raise ParseError(f"FDS truncated before storage_step/start_wavelength doubles (offset {off}, size {len(buf)})")
    storage_step, start_nm = struct.unpack_from("<dd", buf, off)
    if not (0 < storage_step < 5):
        raise ParseError(f"FDS storage_step {storage_step!r} outside expected range (~0.2 nm). "
                         "Format may differ from what spectrex was reverse-engineered against.")
    if not (180 <= start_nm <= 1500):
        raise ParseError(f"FDS start_wavelength {start_nm!r} outside plausible range 180..1500 nm")
    data_off = off + 16

    try:
        reagent = buf.index(b"Reagent 1\x00", data_off)
    except ValueError:
        raise ParseError("footer marker 'Reagent 1\\0' not found after data block. "
                         "Software version may differ; spectrex needs that string to locate data end.")
    data_end = reagent - 32
    rec_size = 40
    if data_end <= data_off:
        raise ParseError(f"computed data_end ({data_end}) precedes data_off ({data_off})")
    if (data_end - data_off) % rec_size != 0:
        raise ParseError(f"data span {data_end - data_off} bytes is not a multiple of "
                         f"record size {rec_size}. Footer offset may have shifted.")
    n_records = (data_end - data_off) // rec_size

    values = []
    for i in range(n_records):
        values.extend(struct.unpack_from("<5d", buf, data_off + i*rec_size))
    step = storage_step

    # The last partial five-value record is padded with four footer doubles.
    # Its final double is the requested emission end wavelength, so retain
    # exactly the samples from start through that endpoint.
    end_nm = struct.unpack_from("<d", buf, data_end - 8)[0]
    expected_points = round((end_nm - start_nm) / step) + 1
    if (180 <= end_nm <= 1500 and end_nm >= start_nm and
            1 <= expected_points <= len(values)):
        values = values[:expected_points]

    wavelengths = [start_nm + i * step for i in range(len(values))]

    # Excitation wavelength: first double of the footer (at data_end)
    excitation_wl = None
    if data_end + 8 <= len(buf):
        v = struct.unpack_from("<d", buf, data_end)[0]
        if 180 <= v <= 1500:
            excitation_wl = v

    return dict(kind="FDS", sample=sample, timestamp=timestamp, operator=operator, comment=comment,
                instrument=instrument, serial=serial, version=rom,
                y_label="fluorescence", wavelengths=wavelengths, values=values,
                step_nm=step, excitation_wl_nm=excitation_wl)

def parse(path):
    try:
        buf = open(path, "rb").read()
    except OSError as e:
        raise ParseError(f"cannot read file: {e}")
    if len(buf) < 16:
        raise ParseError(f"file too small ({len(buf)} bytes) to contain a header")
    magic = buf[:8]
    if magic == b"IIHIITAG": return parse_uds(buf)
    if magic == b"IIHIDTAG": return parse_fds(buf)
    raise ParseError(f"unrecognized magic bytes {magic!r}. Expected b'IIHIITAG' (UDS) or b'IIHIDTAG' (FDS).")

# ---------------- output ----------------

def write_csv(p, csv_path):
    def fmt(v, suffix=""):
        return f"{v}{suffix}" if v is not None else "unknown"
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write(f"# Kind: {p['kind']}\n")
        f.write(f"# Sample: {p['sample']}\n")
        f.write(f"# Timestamp: {p['timestamp']}\n")
        if p.get("operator") is not None:
            f.write(f"# Operator: {p['operator']}\n")
        if p.get("comment"):
            f.write(f"# Comment: {p['comment']}\n")
        f.write(f"# Instrument: {p['instrument']}  SN {p['serial']}  v{p['version']}\n")
        f.write(f"# Points: {len(p['values'])}  Range: {p['wavelengths'][0]:.2f}-{p['wavelengths'][-1]:.2f} nm\n")
        f.write(f"# Sampling step: {fmt(p.get('step_nm'), ' nm')}\n")
        if p["kind"] == "UDS":
            f.write(f"# Slit width: {fmt(p.get('slit_width_nm'), ' nm')}\n")
            f.write(f"# Scan speed: {fmt(p.get('scan_speed_nm_per_min'), ' nm/min')}\n")
            f.write(f"# Path length: {fmt(p.get('path_length_mm'), ' mm')}\n")
            f.write(f"# Lamp change wavelength: {fmt(p.get('lamp_change_wl_nm'), ' nm')}\n")
            f.write(f"# Baseline correction: {fmt(p.get('baseline_correction'))}\n")
            f.write(f"# Response setting: {fmt(p.get('response_setting'))}\n")
        else:  # FDS
            f.write(f"# Excitation wavelength: {fmt(p.get('excitation_wl_nm'), ' nm')}\n")
            f.write(f"# Note: FDS scan speed, slit widths, PMT voltage, response, and delay\n")
            f.write(f"#       are not currently extracted (appear encoded in the binary; see README).\n")
        f.write(f"wavelength_nm,{p['y_label']}\n")
        for wl, v in zip(p["wavelengths"], p["values"]):
            f.write(f"{wl:.2f},{v:.4f}\n")

def write_plot(p, png_path, title):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(p["wavelengths"], p["values"], lw=1.5, color="#e91e63")
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel(p["y_label"].capitalize())
    ax.set_title(title)
    ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(png_path, dpi=130); plt.close(fig)

def write_overlay(group, png_path, ylabel):
    fig, ax = plt.subplots(figsize=(10, 6))
    for name, p in group:
        ax.plot(p["wavelengths"], p["values"], lw=1.4, label=name)
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{ylabel} spectra (overlay)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout(); fig.savefig(png_path, dpi=130); plt.close(fig)

# ---------------- driver ----------------

def find_spectrum_files(inputs):
    """Resolve folders and individual files into a de-duplicated input list."""
    paths = []
    for item in inputs:
        if os.path.isdir(item):
            for ext in ("*.UDS", "*.uds", "*.FDS", "*.fds"):
                paths.extend(glob.glob(os.path.join(item, "**", ext), recursive=True))
        elif os.path.isfile(item) and os.path.splitext(item)[1].lower() in {".uds", ".fds"}:
            paths.append(item)
    return sorted({os.path.normpath(path) for path in paths})


def convert_folder(root, report=print):
    """Convert all supported files below root and return (succeeded, failed)."""
    if not os.path.isdir(root):
        raise ValueError(f"'{root}' is not a directory")
    return convert_paths(find_spectrum_files([root]), root, report)


def convert_paths(paths, root=None, report=print):
    """Convert explicit UDS/FDS paths and return (succeeded, failed)."""
    global LOG_FILE, LOG_SINK
    LOG_SINK = report
    paths = sorted({os.path.normpath(path) for path in paths})
    if root is None:
        try:
            root = os.path.commonpath([os.path.dirname(path) for path in paths])
        except ValueError:
            root = os.path.dirname(paths[0]) if paths else os.getcwd()

    LOG_FILE = os.path.join(root, "spectrex.log")
    with open(LOG_FILE, "w", encoding="utf-8") as f:
        f.write(f"spectrex run at {datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"root: {root}\n\n")

    if not paths:
        log("No .UDS or .FDS files were selected or found in the chosen folders")
        return 0, 0

    n_ok = n_fail = 0
    uds_group, fds_group = [], []
    for path in paths:
        rel = os.path.relpath(path, root)
        # Parse
        try:
            p = parse(path)
        except ParseError as e:
            n_fail += 1
            log(f"FAIL {rel} [parse]: {e}")
            continue
        except Exception as e:
            n_fail += 1
            log(f"FAIL {rel} [parse, UNEXPECTED]: {type(e).__name__}: {e}")
            log("  traceback:\n    " + traceback.format_exc().replace("\n", "\n    ").rstrip())
            continue

        # Sanity checks on parsed result
        if not p.get("values"):
            n_fail += 1
            log(f"FAIL {rel} [parse]: empty data array after parsing")
            continue
        if len(p["values"]) != len(p["wavelengths"]):
            n_fail += 1
            log(f"FAIL {rel} [parse]: length mismatch "
                f"(values={len(p['values'])}, wavelengths={len(p['wavelengths'])})")
            continue

        base = os.path.splitext(path)[0]
        name = os.path.basename(base)

        # Write CSV
        try:
            write_csv(p, base + ".csv")
        except OSError as e:
            n_fail += 1
            log(f"FAIL {rel} [csv]: {e}")
            continue

        # Write PNG
        try:
            write_plot(p, base + ".png", name)
        except Exception as e:
            n_fail += 1
            log(f"FAIL {rel} [png]: {type(e).__name__}: {e}")
            continue

        n_ok += 1
        log(f"OK   {rel}: {p['kind']} {len(p['values'])} pts "
            f"{p['wavelengths'][0]:.0f}-{p['wavelengths'][-1]:.0f} nm")
        (fds_group if p["kind"] == "FDS" else uds_group).append((name, p))

    # Overlays (skipped when a group has fewer than 2 spectra)
    overlays = [(g, lbl, fn) for g, lbl, fn in
                [(uds_group, "Absorbance", "_overlay_uds.png"),
                 (fds_group, "Fluorescence", "_overlay_fds.png")]
                if len(g) >= 2]
    if overlays:
        plots_dir = os.path.join(root, "plots")
        try:
            os.makedirs(plots_dir, exist_ok=True)
        except OSError as e:
            log(f"FAIL [plots dir]: {e}")
        else:
            for group, label, fname in overlays:
                try:
                    write_overlay(group, os.path.join(plots_dir, fname), label)
                    log(f"  -> plots/{fname}")
                except Exception as e:
                    log(f"FAIL [overlay {fname}]: {type(e).__name__}: {e}")

    log(f"\nDone. {n_ok} succeeded, {n_fail} failed. Log: {os.path.relpath(LOG_FILE, root)}")
    return n_ok, n_fail


def launch_gui():
    """Polished native Windows interface used when the EXE is double-clicked."""
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk
    from tkinterdnd2 import DND_FILES, TkinterDnD

    app = TkinterDnD.Tk()
    app.title("SpectreX")
    app.geometry("760x590")
    app.minsize(640, 500)
    app.configure(bg="#f6f7fb")

    style = ttk.Style(app)
    style.theme_use("clam")
    style.configure("App.TFrame", background="#f6f7fb")
    style.configure("Header.TFrame", background="#312e81")
    style.configure("Card.TFrame", background="#ffffff")
    style.configure("Header.Title.TLabel", background="#312e81", foreground="#ffffff",
                    font=("Segoe UI", 23, "bold"))
    style.configure("Header.Subtitle.TLabel", background="#312e81", foreground="#d9d8ff",
                    font=("Segoe UI", 10))
    style.configure("Title.TLabel", background="#ffffff", foreground="#1e1b4b",
                    font=("Segoe UI", 13, "bold"))
    style.configure("Muted.TLabel", background="#ffffff", foreground="#6b7280", font=("Segoe UI", 10))
    style.configure("Status.TLabel", background="#ffffff", foreground="#374151", font=("Segoe UI", 10))
    style.configure("Spectre.Horizontal.TProgressbar", troughcolor="#ede9fe", background="#8b5cf6",
                    lightcolor="#8b5cf6", darkcolor="#8b5cf6", bordercolor="#ede9fe", thickness=8)
    style.configure("Accent.TButton", background="#db2777", foreground="#ffffff",
                    font=("Segoe UI", 10, "bold"), borderwidth=0, padding=(18, 10))
    style.map("Accent.TButton", background=[("active", "#be185d"), ("disabled", "#f2b5d1")])
    style.configure("Quiet.TButton", background="#ffffff", foreground="#4338ca",
                    font=("Segoe UI", 10, "bold"), borderwidth=0, padding=(10, 6))
    style.map("Quiet.TButton", background=[("active", "#eef2ff")])

    selected = tk.StringVar()
    status = tk.StringVar(value="Choose the folder containing your .UDS and/or .FDS files.")
    events = queue.Queue()
    last_folder = [None]
    details_visible = [False]
    selected_inputs = []
    processed = [0]
    selected_total = [0]

    header = ttk.Frame(app, style="Header.TFrame", padding=(32, 24))
    header.pack(fill="x")
    ttk.Label(header, text="SpectreX", style="Header.Title.TLabel").pack(anchor="w")
    ttk.Label(header, text="Turn locked Hitachi spectra into usable CSV files and plots.",
              style="Header.Subtitle.TLabel").pack(anchor="w", pady=(2, 0))

    frame = ttk.Frame(app, style="App.TFrame", padding=(32, 24))
    frame.pack(fill="both", expand=True)

    card = ttk.Frame(frame, style="Card.TFrame", padding=22)
    card.pack(fill="x")
    ttk.Label(card, text="Drop spectra here, or choose what to convert", style="Title.TLabel").pack(anchor="w")
    ttk.Label(card, text="Drop a folder to scan it and its subfolders, or drop individual .UDS and .FDS files.",
              style="Muted.TLabel", wraplength=650).pack(anchor="w", pady=(3, 15))

    chooser = ttk.Frame(card, style="Card.TFrame")
    chooser.pack(fill="x")
    entry = ttk.Entry(chooser, textvariable=selected, font=("Segoe UI", 10))
    entry.pack(side="left", fill="x", expand=True, ipady=7)

    def choose_folder():
        folder = filedialog.askdirectory(title="Choose a folder of spectrum files")
        if folder:
            set_selection([folder])

    def choose_files():
        paths = filedialog.askopenfilenames(
            title="Choose spectrum files",
            filetypes=[("Hitachi spectrum files", "*.uds *.UDS *.fds *.FDS"), ("All files", "*.*")])
        if paths:
            set_selection(list(paths))

    def set_selection(inputs):
        selected_inputs[:] = [os.path.normpath(path) for path in inputs]
        if len(inputs) == 1:
            selected.set(selected_inputs[0])
        else:
            selected.set(f"{len(inputs)} individual spectrum files selected")

    browse = ttk.Button(chooser, text="Browse...", style="Quiet.TButton", command=choose_folder)
    browse.pack(side="left", padx=(8, 0))
    files_button = ttk.Button(chooser, text="Choose files", style="Quiet.TButton", command=choose_files)
    files_button.pack(side="left", padx=(4, 0))

    controls = ttk.Frame(card, style="Card.TFrame")
    controls.pack(fill="x", pady=(18, 0))
    run_button = ttk.Button(controls, text="Convert spectra", style="Accent.TButton", state="disabled")
    run_button.pack(side="left")
    ttk.Label(controls, text="CSV and PNG files are saved beside each original file.",
              style="Muted.TLabel").pack(side="left", padx=(14, 0))

    result_card = ttk.Frame(frame, style="Card.TFrame", padding=22)
    result_card.pack(fill="x", pady=(18, 0))
    result_title = ttk.Label(result_card, text="Ready when you are", style="Title.TLabel")
    result_title.pack(anchor="w")
    ttk.Label(result_card, textvariable=status, style="Status.TLabel", wraplength=650).pack(anchor="w", pady=(4, 12))
    progress_text = tk.StringVar(value="No files selected yet")
    progress = ttk.Progressbar(result_card, style="Spectre.Horizontal.TProgressbar", mode="determinate", maximum=1)
    progress.pack(fill="x")
    ttk.Label(result_card, textvariable=progress_text, style="Muted.TLabel").pack(anchor="w", pady=(5, 12))

    result_actions = ttk.Frame(result_card, style="Card.TFrame")
    result_actions.pack(fill="x")
    open_button = ttk.Button(result_actions, text="Open output folder", style="Quiet.TButton", state="disabled")
    open_button.pack(side="left")
    details_button = ttk.Button(result_actions, text="Show details", style="Quiet.TButton")
    details_button.pack(side="left", padx=(8, 0))

    output = tk.Text(frame, height=10, wrap="word", state="disabled", font=("Cascadia Mono", 9),
                     bg="#15143f", fg="#e9e8ff", insertbackground="#ffffff", relief="flat", padx=12, pady=10)

    def append(text):
        output.configure(state="normal")
        output.insert("end", text + "\n")
        output.see("end")
        output.configure(state="disabled")

    def toggle_details():
        details_visible[0] = not details_visible[0]
        if details_visible[0]:
            output.pack(fill="both", expand=True, pady=(14, 0))
            details_button.configure(text="Hide details")
        else:
            output.pack_forget()
            details_button.configure(text="Show details")

    def open_output_folder():
        if last_folder[0]:
            os.startfile(last_folder[0])

    def receive_drop(event):
        inputs = list(app.tk.splitlist(event.data))
        usable = [path for path in inputs if os.path.isdir(path) or
                  (os.path.isfile(path) and os.path.splitext(path)[1].lower() in {".uds", ".fds"})]
        if not usable:
            messagebox.showwarning("SpectreX", "Drop a folder or one or more .UDS / .FDS files.")
            return
        set_selection(usable)
        status.set("Ready to convert the dropped selection.")

    def poll_events():
        try:
            while True:
                kind, payload = events.get_nowait()
                if kind == "log":
                    append(payload)
                    if payload.startswith("OK   ") or payload.startswith("FAIL "):
                        processed[0] += 1
                        progress.configure(value=processed[0])
                        progress_text.set(f"Processed {processed[0]} of {selected_total[0]} file(s)")
                else:
                    ok, failed, error = payload
                    browse.configure(state="normal"); files_button.configure(state="normal")
                    run_button.configure(state="normal" if selected.get().strip() else "disabled")
                    if error:
                        result_title.configure(text="Conversion could not start")
                        status.set("Conversion could not start.")
                        messagebox.showerror("SpectreX", error)
                    else:
                        result_title.configure(text="Conversion complete")
                        status.set(f"{ok} file(s) converted, {failed} failed. CSV and PNG files are beside the originals.")
                        open_button.configure(state="normal")
                        if failed:
                            details_button.configure(text=f"Show {failed} error detail{'s' if failed != 1 else ''}")
                            if not details_visible[0]:
                                toggle_details()
        except queue.Empty:
            pass
        app.after(100, poll_events)

    def start_conversion():
        if not selected_inputs:
            messagebox.showwarning("SpectreX", "Choose or drop a folder or spectrum file first.")
            return
        paths = find_spectrum_files(selected_inputs)
        if not paths:
            messagebox.showwarning("SpectreX", "No .UDS or .FDS files were found in that selection.")
            return
        try:
            root = os.path.commonpath([os.path.dirname(path) for path in paths])
        except ValueError:
            root = os.path.dirname(paths[0])
        output.configure(state="normal"); output.delete("1.0", "end"); output.configure(state="disabled")
        processed[0] = 0
        selected_total[0] = len(paths)
        progress.configure(maximum=len(paths), value=0)
        progress_text.set(f"Preparing {len(paths)} file(s)")
        browse.configure(state="disabled"); files_button.configure(state="disabled"); run_button.configure(state="disabled")
        open_button.configure(state="disabled")
        result_title.configure(text="Converting spectra")
        status.set("Working through the selected folder. You can expand details to follow each file.")
        last_folder[0] = root
        def worker():
            try:
                ok, failed = convert_paths(paths, root, report=lambda m: events.put(("log", m)))
                events.put(("done", (ok, failed, None)))
            except Exception as exc:
                events.put(("done", (0, 0, str(exc))))
        threading.Thread(target=worker, daemon=True).start()

    run_button.configure(command=start_conversion)
    open_button.configure(command=open_output_folder)
    details_button.configure(command=toggle_details)
    for widget in (card, entry, chooser):
        widget.drop_target_register(DND_FILES)
        widget.dnd_bind("<<Drop>>", receive_drop)
    selected.trace_add("write", lambda *_: run_button.configure(state="normal" if selected.get().strip() else "disabled"))
    app.after(100, poll_events)
    app.mainloop()


def main():
    if len(sys.argv) > 1:
        try:
            convert_folder(sys.argv[1])
        except ValueError as e:
            print(f"ERROR: {e}")
            sys.exit(1)
    else:
        launch_gui()

if __name__ == "__main__":
    main()
