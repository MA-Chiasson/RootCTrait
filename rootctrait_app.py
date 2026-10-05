"""rootctrait_app.py — RootCTrait desktop app.

Features: session saving (including the parent folder), interactive Treeview table with
striped rows and scrolling, green progress bar with percentage, ETA estimate, console
log saving, and a log size guard that prevents the interface from freezing.
"""
import os
import json
import queue
import time
import threading
import webbrowser
import subprocess
import platform
import tkinter as tk
import tkinter.font as tkfont
from tkinter import ttk, scrolledtext, filedialog, messagebox

import pipeline_api as api

CONFIG_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "app_session.json")

def _ensure_dnd():
    try:
        import tkinterdnd2
        return True
    except Exception:
        pass
    import tkinter as _tk
    splash = _tk.Tk()
    splash.title("RootCTrait")
    splash.geometry("360x90")
    _tk.Label(splash, text="First launch: installing drag-and-drop support...\n(one-time, a few seconds)", padx=12, pady=16).pack()
    splash.update()
    ok = False
    try:
        import subprocess, sys
        subprocess.run([sys.executable, "-m", "pip", "install", "tkinterdnd2"], check=True, capture_output=True, timeout=180)
        import importlib
        importlib.invalidate_caches()
        import tkinterdnd2
        ok = True
    except Exception:
        ok = False
    splash.destroy()
    return ok

_DND = _ensure_dnd()
if _DND:
    from tkinterdnd2 import TkinterDnD, DND_FILES


class Tooltip:
    """Small information window shown when hovering over a widget,
    hidden as soon as the cursor leaves it."""

    def __init__(self, widget, text, delay=450, wraplength=280):
        self.widget = widget
        self.text = text
        self.delay = delay
        self.wraplength = wraplength
        self.tipwindow = None
        self._after_id = None
        widget.bind("<Enter>", self._schedule, add="+")
        widget.bind("<Leave>", self._hide, add="+")
        widget.bind("<ButtonPress>", self._hide, add="+")

    def _schedule(self, _event=None):
        self._cancel()
        self._after_id = self.widget.after(self.delay, self._show)

    def _cancel(self):
        if self._after_id is not None:
            try:
                self.widget.after_cancel(self._after_id)
            except Exception:
                pass
            self._after_id = None

    def _show(self):
        if self.tipwindow or not self.text:
            return
        x = self.widget.winfo_rootx() + 20
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 6
        self.tipwindow = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(True)
        tw.wm_geometry(f"+{x}+{y}")
        tw.attributes("-topmost", True)
        frame = tk.Frame(tw, background="#2d3748", borderwidth=0)
        frame.pack()
        label = tk.Label(
            frame, text=self.text, justify="left",
            background="#2d3748", foreground="#f7fafc",
            font=("Segoe UI", 9), wraplength=self.wraplength,
            padx=10, pady=7,
        )
        label.pack()

    def _hide(self, _event=None):
        self._cancel()
        if self.tipwindow is not None:
            try:
                self.tipwindow.destroy()
            except Exception:
                pass
            self.tipwindow = None


class ImportWindow(tk.Toplevel):
    def __init__(self, master, current, initial_parent, on_ok):
        super().__init__(master)
        self.master_app = master
        self.title("Import folders")
        self.geometry("650x520")
        self.configure(bg="#f8f9fa")
        self.on_ok = on_ok
        self.folders = list(current)
        self.parent_folder = initial_parent

        # Bring the import window to the front
        self.attributes("-topmost", True)
        self.grab_set()

        note = ("Import folders for analysis. Click 'Browse Parent Folder...' to choose a root folder "
                "(like 'data'), then click on rows to check/uncheck the batch folders you want.")
        ttk.Label(self, text=note, wraplength=610, style="Info.TLabel").pack(fill="x", padx=20, pady=(20, 10))

        # --- Top area: parent folder selection ---
        top_frame = ttk.Frame(self, padding=(20, 5))
        top_frame.pack(fill="x")
        ttk.Button(top_frame, text="📁 Browse Parent Folder...", command=self._browse_parent, style="Accent.TButton").pack(side="left")
        
        display_name = os.path.basename(self.parent_folder) if self.parent_folder else "No parent folder selected"
        self.parent_lbl = ttk.Label(top_frame, text=display_name, font=("Segoe UI", 9, "italic"), foreground="#718096")
        self.parent_lbl.pack(side="left", padx=15)

        # --- Central area: Treeview table with graphical check marks ---
        mid = ttk.Frame(self, padding=20)
        mid.pack(fill="both", expand=True)
        
        self.tree = ttk.Treeview(mid, columns=("Batch", "Path"), show="tree headings", selectmode="none")
        self.tree.heading("#0", text="Import", anchor="center")
        self.tree.heading("Batch", text="Batch Name")
        self.tree.heading("Path", text="Full Folder Path")
        self.tree.column("#0", width=70, anchor="center")
        self.tree.column("Batch", width=140, anchor="w")
        self.tree.column("Path", width=360, anchor="w")
        self.tree.pack(side="left", fill="both", expand=True)
        
        sb = ttk.Scrollbar(mid, command=self.tree.yview)
        sb.pack(side="left", fill="y")
        self.tree.configure(yscrollcommand=sb.set)

        # Blue colour for checked folders
        self.tree.tag_configure("checked", background="#3182ce", foreground="white")

        # Click binding to check/uncheck
        self.tree.bind("<ButtonRelease-1>", self._on_tree_click)

        if _DND:
            self.tree.drop_target_register(DND_FILES)
            self.tree.dnd_bind("<<Drop>>", self._on_drop)

        # --- Bottom area: action buttons ---
        btns = ttk.Frame(self, padding=20)
        btns.pack(fill="x")
        
        ttk.Button(btns, text="Select All", command=self._select_all).pack(side="left")
        ttk.Button(btns, text="Deselect All", command=self._deselect_all).pack(side="left", padx=10)
        ttk.Button(btns, text="Remove Selected", command=self._remove).pack(side="left")
        
        ttk.Button(btns, text="Cancel", command=self.destroy).pack(side="right")
        ttk.Button(btns, text="OK", command=self._ok, style="Accent.TButton").pack(side="right", padx=10)

        hint = "  (Drag-and-drop folders directly into the list is also supported)" if _DND else ""
        self.hint_lbl = ttk.Label(self, text=f"{len(self.folders)} folder(s) selected{hint}", style="Hint.TLabel")
        self.hint_lbl.pack(anchor="w", padx=20, pady=(0, 15))
        
        # If a parent folder was saved, load it automatically at startup
        if self.parent_folder and os.path.exists(self.parent_folder):
            self._load_subfolders(self.parent_folder)
        else:
            self._refresh()

    def _load_subfolders(self, root_abs):
        try:
            for entry in sorted(os.listdir(root_abs)):
                full_path = os.path.join(root_abs, entry)
                if os.path.isdir(full_path):
                    if full_path not in self.folders:
                        self.tree.insert("", "end", text="   ", values=(entry, full_path))
        except Exception as e:
            messagebox.showerror("Error", f"Could not read folder content: {e}", parent=self)
        self._refresh()

    def _refresh(self):
        """Refresh the display while keeping the initial, stable alphabetical order."""
        current_items = {}
        for item in self.tree.get_children():
            values = self.tree.item(item, "values")
            if values:
                current_items[values[1]] = values[0]

        for f in self.folders:
            current_items[f] = os.path.basename(f)

        self.tree.delete(*self.tree.get_children())
        
        for path, batch in sorted(current_items.items()):
            if path in self.folders:
                self.tree.insert("", "end", text="✔", values=(batch, path), tags=("checked",))
            else:
                self.tree.insert("", "end", text="   ", values=(batch, path))

        if hasattr(self, 'hint_lbl'):
            self.hint_lbl.configure(text=f"{len(self.folders)} folder(s) selected")

    def _browse_parent(self):
        self.attributes("-topmost", False)
        from tkinter import filedialog
        root_dir = filedialog.askdirectory(title="Select the Parent Folder (e.g., 'data')", parent=self)
        
        if root_dir:
            root_abs = os.path.abspath(root_dir)
            self.parent_folder = root_abs
            self.master_app.parent_folder = root_abs
            self.parent_lbl.configure(text=os.path.basename(root_abs))
            
            self.tree.delete(*self.tree.get_children())
            self._load_subfolders(root_abs)
            
        self.attributes("-topmost", True)
        self.lift()
        self.focus_force()

    def _on_tree_click(self, event):
        item = self.tree.identify_row(event.y)
        if not item: return
        
        values = self.tree.item(item, "values")
        if not values or len(values) < 2: return
        
        batch, path = values
        status = self.tree.item(item, "text")
        
        if "✔" not in status:
            if path not in self.folders:
                self.folders.append(path)
        else:
            if path in self.folders:
                self.folders.remove(path)
                
        self._refresh()

    def _select_all(self):
        for item in self.tree.get_children():
            values = self.tree.item(item, "values")
            if values and len(values) > 1:
                path = values[1]
                if path not in self.folders:
                    self.folders.append(path)
        self._refresh()

    def _deselect_all(self):
        self.folders.clear()
        self._refresh()

    def _on_drop(self, event):
        for p in self.tk.splitlist(event.data):
            p_abs = os.path.abspath(p)
            if os.path.isdir(p_abs) and p_abs not in self.folders:
                self.folders.append(p_abs)
        self._refresh()

    def _remove(self):
        for item in self.tree.get_children():
            values = self.tree.item(item, "values")
            if values and len(values) > 1:
                path = values[1]
                if path in self.folders:
                    self.folders.remove(path)
        self.tree.delete(*self.tree.get_children())
        self._refresh()

    def _ok(self):
        if not self.folders:
            messagebox.showwarning("No folder", "Please select at least one folder by clicking on it.", parent=self)
            return
        self.on_ok(self.folders)
        self.destroy()

PARAM_FIELDS = [
    ("RESULTS_ROOT", "Output folder (results root)", "path",  True),
    ("VOXEL",        "Voxel size (depth,x,z) mm",    "text",  True),
    ("PATTERN",      "File name pattern",            "text",  True),
    ("PRUNE_VOX",    "Prune length (voxels)",        "int",   False),
    ("MIN_SEG_LEN_MM","Min segment length (mm)",     "float", False),
    ("BC_MIN",       "Sheet: parallel neighbours",   "int",   False),
    ("LIN_MAX",      "Sheet: max linearity",         "float", False),
    ("LEN_MAX",      "Sheet: max length (mm)",       "float", False),
    ("TIMEOUT",      "Timeout per sample (s)",       "int",   False),
    ("DENS_MAX",     "Dense rule: min density",      "float", False),
    ("DENS_LEN_MAX", "Dense rule: max length (mm)",  "float", False),
    ("RESCUE_MIN_MM","Rescue: min length (mm)",      "float", False),
    ("DROP_ORPHANS", "Drop orphan fragments",        "bool",  False),
    ("SAVE_FIGURES", "Save 3D figures",              "bool",  False),
    ("EXPORT_RSML",  "Export RSML files",            "bool",  False),
]
REQUIRED = {k for k, _, _, req in PARAM_FIELDS if req}
PRESET_EXCLUDE = {"RESULTS_ROOT"}
PRESETS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "params")

def list_presets():
    names = ["Default"]
    if os.path.isdir(PRESETS_DIR):
        for f in sorted(os.listdir(PRESETS_DIR)):
            if f.endswith(".json"):
                names.append(f[:-5])
    return names

def load_preset(name, base_defaults):
    vals = dict(base_defaults)
    if name != "Default":
        try:
            with open(os.path.join(PRESETS_DIR, name + ".json"), encoding="utf-8") as fh:
                vals.update(json.load(fh))
        except Exception:
            pass
    return vals

def save_preset(name, params):
    os.makedirs(PRESETS_DIR, exist_ok=True)
    data = {k: v for k, v in params.items() if k not in PRESET_EXCLUDE}
    with open(os.path.join(PRESETS_DIR, name + ".json"), "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)

def delete_preset(name):
    path = os.path.join(PRESETS_DIR, name + ".json")
    if os.path.exists(path):
        os.remove(path)

class ParamsWindow(tk.Toplevel):
    def __init__(self, master, current, on_ok):
        super().__init__(master)
        self.title("Analysis parameters")
        self.geometry("620x660")
        self.configure(bg="#f8f9fa")
        self.on_ok = on_ok
        # Keep the window above the main window and hold the focus.
        self.transient(master)
        self.lift()
        self.grab_set()
        self.base_defaults = api.default_params()
        self.vars = {}

        pf = ttk.Frame(self, padding=(20, 20, 20, 5))
        pf.pack(fill="x")
        ttk.Label(pf, text="Parameter set:").pack(side="left")
        self.preset_var = tk.StringVar(value="Default")
        self.combo = ttk.Combobox(pf, textvariable=self.preset_var, state="readonly", values=list_presets(), width=22)
        self.combo.pack(side="left", padx=10)
        self.combo.bind("<<ComboboxSelected>>", self._on_preset)
        ttk.Button(pf, text="Save as...", command=self._save_as).pack(side="left", padx=5)
        ttk.Button(pf, text="Delete", command=self._delete).pack(side="left")

        self.error_lbl = ttk.Label(self, text="", font=("Segoe UI", 10, "bold"), foreground="#e53e3e", background="#f8f9fa")
        self.error_lbl.pack(fill="x", padx=20, pady=(5, 0))

        ttk.Label(self, text="Fields marked * are required. The output folder descriptions describe a run, not an experiment type.", wraplength=580, style="Hint.TLabel").pack(fill="x", padx=20, pady=(5, 15))

        form_container = ttk.Frame(self, padding=(20, 0))
        form_container.pack(fill="both", expand=True)
        
        canvas = tk.Canvas(form_container, bg="#f8f9fa", highlightthickness=0)
        scrollbar = ttk.Scrollbar(form_container, orient="vertical", command=canvas.yview)
        form = ttk.Frame(canvas, padding=5)
        
        form.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=form, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        for i, (key, label, kind, req) in enumerate(PARAM_FIELDS):
            text = label + (" *" if req else "")
            ttk.Label(form, text=text, font=("Segoe UI", 10)).grid(row=i, column=0, sticky="w", pady=6, padx=(0, 15))
            val = current.get(key, "")
            if kind == "bool":
                v = tk.BooleanVar(value=bool(val))
                ttk.Checkbutton(form, variable=v, command=self._validate_fields).grid(row=i, column=1, sticky="w", pady=6)
            else:
                v = tk.StringVar(value=str(val))
                entry = ttk.Entry(form, textvariable=v, width=34)
                entry.grid(row=i, column=1, sticky="we", pady=6)
                entry.bind("<KeyRelease>", lambda e: self._validate_fields())
                if kind == "path":
                    ttk.Button(form, text="Browse...", command=lambda vv=v: self._browse(vv)).grid(row=i, column=2, padx=10, pady=6)
            self.vars[key] = (v, kind)
        form.columnconfigure(1, weight=1)

        self.btns_frame = ttk.Frame(self, padding=20)
        self.btns_frame.pack(fill="x")
        ttk.Button(self.btns_frame, text="Cancel", command=self.destroy).pack(side="right")
        self.ok_btn = ttk.Button(self.btns_frame, text="OK", command=self._ok, style="Accent.TButton")
        self.ok_btn.pack(side="right", padx=10)
        
        self._validate_fields()

    def _validate_fields(self):
        missing = []
        for key, label, _, req in PARAM_FIELDS:
            if req:
                v, _ = self.vars[key]
                if not str(v.get()).strip():
                    missing.append(label)
        if missing:
            self.error_lbl.configure(text=f"⚠️ Missing required: {', '.join(missing)}")
            self.ok_btn.configure(state="disabled")
        else:
            self.error_lbl.configure(text="")
            self.ok_btn.configure(state="normal")

    def _browse(self, var):
        d = filedialog.askdirectory(title="Select the output folder", parent=self)
        if d:
            var.set(d)
            self._validate_fields()
        # Bring the Parameters window back to the front after the dialog.
        self.lift()
        self.focus_force()

    def _collect(self):
        out = {}
        for key, (v, kind) in self.vars.items():
            raw = v.get()
            if kind == "int":
                out[key] = int(raw) if str(raw).strip() != "" else ""
            elif kind == "float":
                out[key] = float(raw) if str(raw).strip() != "" else ""
            elif kind == "bool":
                out[key] = bool(raw)
            else:
                out[key] = raw
        return out

    def _fill(self, params):
        for key, (v, kind) in self.vars.items():
            if key in PRESET_EXCLUDE: continue
            if key in params:
                v.set(bool(params[key]) if kind == "bool" else str(params[key]))
        self._validate_fields()

    def _on_preset(self, _evt=None):
        self._fill(load_preset(self.preset_var.get(), self.base_defaults))

    def _save_as(self):
        from tkinter import simpledialog
        name = simpledialog.askstring("Save parameter set", "Name for this set:", parent=self)
        if not name: return
        name = name.strip()
        if name.lower() == "default" or not name:
            messagebox.showerror("Invalid name", "'Default' is reserved; choose another name.")
            return
        try:
            save_preset(name, self._collect())
        except ValueError:
            messagebox.showerror("Invalid value", "Check numeric fields before saving.")
            return
        self.combo.configure(values=list_presets())
        self.preset_var.set(name)

    def _delete(self):
        name = self.preset_var.get()
        if name == "Default":
            messagebox.showinfo("Default", "The Default set cannot be deleted.")
            return
        if messagebox.askyesno("Delete?", f"Delete the parameter set '{name}'?"):
            delete_preset(name)
            self.combo.configure(values=list_presets())
            self.preset_var.set("Default")
            self._on_preset()

    def _ok(self):
        try:
            out = self._collect()
        except ValueError:
            self.error_lbl.configure(text="⚠️ Please check numeric fields for invalid values.")
            return
        self.on_ok(out)
        self.destroy()

class App(TkinterDnD.Tk if _DND else tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("RootCTrait")
        self.geometry("960x780")
        self.minsize(760, 560)
        self.configure(bg="#f1f3f5")
        
        self.msg_q = queue.Queue()
        self.running = False
        self.stop_requested = False
        self.folders = []
        self.parent_folder = ""
        self.params = api.default_params()
        self.params["RESULTS_ROOT"] = ""
        
        self.start_time = None
        self.total_samples = 0
        self.current_done = 0

        self.style = ttk.Style()
        self.style.theme_use("clam")

        # --- Named fonts: resized together with the window ---
        self.f_base       = tkfont.Font(family="Segoe UI", size=10)
        self.f_button     = tkfont.Font(family="Segoe UI", size=9,  weight="bold")
        self.f_tab        = tkfont.Font(family="Segoe UI", size=10)
        self.f_tab_sel    = tkfont.Font(family="Segoe UI", size=11, weight="bold")
        self.f_status     = tkfont.Font(family="Segoe UI", size=12, weight="bold")
        self.f_hint       = tkfont.Font(family="Segoe UI", size=9,  slant="italic")
        self.f_info       = tkfont.Font(family="Segoe UI", size=10)
        self.f_section    = tkfont.Font(family="Segoe UI", size=10, weight="bold")
        self.f_tree       = tkfont.Font(family="Segoe UI", size=10)
        self.f_tree_head  = tkfont.Font(family="Segoe UI", size=10, weight="bold")
        self.f_console    = tkfont.Font(family="Consolas", size=10)
        self.f_console_b  = tkfont.Font(family="Consolas", size=10, weight="bold")

        # (font, base size) for proportional resizing
        self._scalable_fonts = [
            (self.f_base, 10), (self.f_button, 9), (self.f_tab, 10), (self.f_tab_sel, 11),
            (self.f_status, 12), (self.f_hint, 9), (self.f_info, 10), (self.f_section, 10),
            (self.f_tree, 10), (self.f_tree_head, 10), (self.f_console, 10), (self.f_console_b, 10),
        ]
        self._base_rowheight = 26
        self._last_scale = None

        self.style.configure(".", background="#f1f3f5", foreground="#212529", font=self.f_base)
        self.style.configure("TFrame", background="#f1f3f5")

        self.style.configure("TNotebook", background="#f1f3f5", borderwidth=0)
        self.style.configure("TNotebook.Tab", background="#e9ecef", foreground="#495057", padding=(15, 4), font=self.f_tab)
        self.style.map("TNotebook.Tab",
                       background=[("selected", "#ffffff")],
                       foreground=[("selected", "#1a365d")],
                       font=[("selected", self.f_tab_sel)],
                       padding=[("selected", (18, 6))])

        self.style.configure("TButton", background="#e2e8f0", foreground="#2d3748", borderwidth=0, padding=(12, 6), font=self.f_button)
        self.style.map("TButton", background=[("active", "#cbd5e0"), ("disabled", "#edf2f7")], foreground=[("disabled", "#a0aec0")])
        self.style.configure("Accent.TButton", background="#3182ce", foreground="white")
        self.style.map("Accent.TButton", background=[("active", "#2b6cb0"), ("disabled", "#edf2f7")], foreground=[("disabled", "#a0aec0")])
        self.style.configure("Danger.TButton", background="#e53e3e", foreground="white")
        self.style.map("Danger.TButton", background=[("active", "#c53030"), ("disabled", "#edf2f7")])
        self.style.configure("Success.TButton", background="#48bb78", foreground="white")
        self.style.map("Success.TButton", background=[("active", "#38a169")])

        self.style.configure("Status.TLabel", font=self.f_status, foreground="#c05621")
        self.style.configure("Hint.TLabel", font=self.f_hint, foreground="#718096", background="#f8f9fa")
        self.style.configure("Info.TLabel", font=self.f_info, foreground="#2d3748", background="#f8f9fa")
        self.style.configure("TCombobox", fieldbackground="white", background="#e2e8f0")
        self.style.configure("TEntry", fieldbackground="white")

        self.style.configure("Treeview", font=self.f_tree, rowheight=26, background="white", fieldbackground="white", borderwidth=0)
        self.style.configure("Treeview.Heading", font=self.f_tree_head, background="#e2e8f0", relief="flat")
        self.style.map("Treeview", background=[("selected", "#bee3f8")], foreground=[("selected", "#2b6cb0")])

        # Top container spanning the full width
        header_frame = ttk.Frame(self, padding=(10, 10, 10, 0))
        header_frame.pack(fill="x")
        
        nb = ttk.Notebook(header_frame)
        nb.pack(fill="x", expand=True, pady=(5, 0))

        # --- Tab Analysis ---
        tab_an = ttk.Frame(nb, padding=15)
        nb.add(tab_an, text="Analysis")
        
        # --- LEFT BLOCK: configuration and preparation ---
        left_frame = ttk.Frame(tab_an)
        left_frame.pack(side="left", fill="y")
        
        ttk.Button(left_frame, text="🔄 New Session", command=self.clear_session).pack(side="left")
        
        sep = ttk.Separator(left_frame, orient="vertical")
        sep.pack(side="left", fill="y", padx=12, pady=2)
        
        ttk.Button(left_frame, text="📂 Import Folders", command=self.open_import).pack(side="left")
        ttk.Button(left_frame, text="⚙️ Parameters", command=self.open_params).pack(side="left", padx=(8, 0))
        
        # --- CENTRAL BLOCK: run controls (forced to the exact centre) ---
        # Invisible expandable frame left of the centre
        spacer_left = ttk.Frame(tab_an)
        spacer_left.pack(side="left", fill="x", expand=True)
        
        center_frame = ttk.Frame(tab_an)
        center_frame.pack(side="left", fill="y")
        
        self.btn = ttk.Button(center_frame, text="▶ Analyse", command=self.start_analysis, state="disabled", style="Accent.TButton")
        self.btn.pack(side="left", padx=4)
        
        self.btn_stop = ttk.Button(center_frame, text="⏹ Stop", command=self.stop_analysis, state="disabled", style="Danger.TButton")
        self.btn_stop.pack(side="left", padx=4)
        
        # Invisible expandable frame right of the centre
        spacer_right = ttk.Frame(tab_an)
        spacer_right.pack(side="left", fill="x", expand=True)
        
        # --- BLOC DROIT : Statut dynamique ---
        right_frame = ttk.Frame(tab_an)
        right_frame.pack(side="right", fill="y")
        
        self.status = ttk.Label(right_frame, text="Import folders to start", style="Status.TLabel")
        self.status.pack(side="right", padx=(0, 10))

        # --- Tab Tools ---
        tab_tools = ttk.Frame(nb, padding=15)
        nb.add(tab_tools, text="Tools")

        # Single compact row (same height as the Analysis tab). Full names, no
        # truncation. The format selector applies to the tools that write tables.
        ttk.Label(tab_tools, text="Table format:").pack(side="left", padx=(0, 5))
        self.fmt_var = tk.StringVar(value="both")
        ttk.Combobox(tab_tools, textvariable=self.fmt_var, state="readonly",
                     values=["both", "xlsx", "csv"], width=8).pack(side="left", padx=(0, 4))

        btn_merge = ttk.Button(tab_tools, text="🗜 Merge batches", command=self.tool_merge)
        btn_merge.pack(side="left", padx=4)
        Tooltip(btn_merge, "Merges the trait tables of all batches into a single combined file, "
                           "using the format selected above (both, xlsx or csv).")

        btn_ckpt = ttk.Button(tab_tools, text="🔍 Extract checkpoints", command=self.tool_checkpoint)
        btn_ckpt.pack(side="left", padx=4)
        Tooltip(btn_ckpt, "Extracts traits from the checkpoints (checkpoint_traits.jsonl) without rerunning the analysis, "
                          "using the format selected above (both, xlsx or csv).")

        ttk.Separator(tab_tools, orient="vertical").pack(side="left", fill="y", padx=12, pady=2)

        btn_report = ttk.Button(tab_tools, text="🖼 Figure report", command=self.tool_report)
        btn_report.pack(side="left", padx=4)
        Tooltip(btn_report, "Ranks the samples by their quality control values (qc_ranking.csv) and "
                            "generates the HTML report of interactive 3D Plotly figures, by batch, "
                            "most suspect samples first, for visual validation.")

        btn_del = ttk.Button(tab_tools, text="🗑 Delete checkpoints", command=self.tool_delete_checkpoints)
        btn_del.pack(side="left", padx=4)
        Tooltip(btn_del, "Deletes the checkpoints (checkpoint_traits.jsonl) of a results folder to start over: the "
                         "progress is reset to 0 and the analysis restarts from the beginning.")


        # Import table area with alternating striped rows and scrolling
        table_frame = ttk.Frame(self, padding=(10, 5, 10, 0))
        table_frame.pack(fill="both", expand=True)
        ttk.Label(table_frame, text="Current Import Queue & Status", font=self.f_section, foreground="#4a5568").pack(anchor="w", pady=(0, 3))

        tree_scroll_frame = ttk.Frame(table_frame)
        tree_scroll_frame.pack(fill="both", expand=True)

        self.queue_tree = ttk.Treeview(tree_scroll_frame, columns=("Batch", "Path", "Status"), show="headings", height=4)
        self.queue_tree.heading("Batch", text="Batch Name")
        self.queue_tree.heading("Path", text="Folder Path")
        self.queue_tree.heading("Status", text="Status")
        self.queue_tree.column("Batch", width=150, anchor="w")
        self.queue_tree.column("Path", width=600, anchor="w")
        self.queue_tree.column("Status", width=150, anchor="center")
        self.queue_tree.pack(side="left", fill="both", expand=True)
        
        # Alternating stripe colours
        self.queue_tree.tag_configure("even", background="#ffffff")
        self.queue_tree.tag_configure("odd", background="#f7f9fa")
        
        queue_sb = ttk.Scrollbar(tree_scroll_frame, orient="vertical", command=self.queue_tree.yview)
        queue_sb.pack(side="right", fill="y")
        self.queue_tree.configure(yscrollcommand=queue_sb.set)

        # --- FIX: restore the green progress bar used during the run (tool 3) ---
        progress_frame = ttk.Frame(self, padding=(10, 10, 10, 0))
        progress_frame.pack(fill="x")
        self.progress_bar = ttk.Progressbar(progress_frame, orient="horizontal", mode="determinate")
        self.progress_bar.pack(fill="x")

        # Bottom console area (safe layout)
        term_frame = ttk.Frame(self, padding=10)
        term_frame.pack(fill="both", expand=True)
        
        # 1. The console title stays alone at the top
        ttk.Label(term_frame, text="Live Console Output", font=self.f_section, foreground="#4a5568").pack(anchor="w", pady=(0, 5))

        # 2. The bottom action bar is packed SECOND with side="bottom"
        bottom_action_bar = ttk.Frame(term_frame, padding=(0, 6, 0, 0))
        bottom_action_bar.pack(side="bottom", fill="x")
        
        # Plain button without the dots, aligned to the right edge (padx=2)
        ttk.Button(bottom_action_bar, text="💾 Save Console Log", command=self.save_console_log).pack(side="right", padx=(0, 2))

        # 3. La console noire prend TOUT l'espace central restant
        self.term = scrolledtext.ScrolledText(term_frame, bg="#1a202c", fg="#e2e8f0", insertbackground="#e2e8f0", font=self.f_console, wrap="word", bd=0, highlightthickness=1, highlightbackground="#4a5568")
        self.term.pack(side="top", fill="both", expand=True)
        self.term.configure(state="disabled")
        
        self.term.tag_configure("error", foreground="#f56565")
        self.term.tag_configure("success", foreground="#48bb78")
        self.term.tag_configure("heading", foreground="#63b3ed", font=self.f_console_b)
        
        self._load_session()
        self.protocol("WM_DELETE_WINDOW", self._on_closing)
        self.after(80, self._drain_queue)
        self.bind("<Configure>", self._on_resize)

    def _on_resize(self, event):
        """Resize the fonts in proportion to the window size."""
        if event.widget is not self:
            return
        # Factor based on the 960x780 reference geometry, bounded to stay readable.
        scale = min(event.width / 960.0, event.height / 780.0)
        scale = max(0.85, min(scale, 1.35))
        if self._last_scale is not None and abs(scale - self._last_scale) < 0.02:
            return
        self._last_scale = scale
        for font_obj, base in self._scalable_fonts:
            font_obj.configure(size=max(7, int(round(base * scale))))
        try:
            self.style.configure("Treeview", rowheight=int(round(self._base_rowheight * scale)))
        except Exception:
            pass


    def _load_session(self):
        if os.path.exists(CONFIG_FILE):
            try:
                with open(CONFIG_FILE, "r", encoding="utf-8") as f:
                    data = json.load(f)
                    self.folders = data.get("folders", [])
                    self.parent_folder = data.get("parent_folder", "")
                    saved_params = data.get("params", {})
                    self.params.update(saved_params)
                    if self.folders:
                        self.btn.configure(state="normal")
                        self._update_queue_table()
                        if not self.params.get("RESULTS_ROOT", ""):
                            self.status.configure(text="⚙️ Set output folder in Parameters")
                        else:
                            self.status.configure(text=f"📊 {len(self.folders)} folder(s) ready (restored)")
                        self.log(f"Session restored: loaded {len(self.folders)} folder(s).")
            except Exception as e:
                self.log(f"Could not restore previous session: {e}")

    def _on_closing(self):
        try:
            with open(CONFIG_FILE, "w", encoding="utf-8") as f:
                json.dump({"folders": self.folders, "parent_folder": self.parent_folder, "params": self.params}, f, indent=2, ensure_ascii=False)
        except Exception:
            pass
        self.destroy()

    def _update_queue_table(self, active_batch="", finished_batches=None):
        finished_batches = finished_batches or []
        for item in self.queue_tree.get_children():
            self.queue_tree.delete(item)
        for idx, f in enumerate(self.folders):
            bname = os.path.basename(f)
            if bname == active_batch:
                status_text = "⏳ Processing..."
            elif bname in finished_batches:
                status_text = "✅ Done"
            elif self.running:
                status_text = "💤 Queued"
            else:
                status_text = "📋 Ready"
            
            # Alternate tags to obtain the grey/white striping
            row_tag = "even" if idx % 2 == 0 else "odd"
            self.queue_tree.insert("", "end", values=(bname, f, status_text), tags=(row_tag,))

    def open_results_folder(self):
        path = self.params.get("RESULTS_ROOT", "")
        if not path or not os.path.exists(path):
            messagebox.showerror("Error", "Results folder does not exist yet.")
            return
        current_os = platform.system()
        try:
            if current_os == "Windows":
                os.startfile(path)
            elif current_os == "Darwin":
                subprocess.run(["open", path], check=True)
            else:
                subprocess.run(["xdg-open", path], check=True)
        except Exception as e:
            messagebox.showerror("Error", f"Could not open folder: {e}")

    def log(self, line):
        self.term.configure(state="normal")
        
        # --- Guard against freezing / memory saturation ---
        # Above 2000 log lines, purge the oldest lines to relieve Tkinter and avoid freezing
        if float(self.term.index("end-1c")) > 2000.0:
            self.term.delete("1.0", "100.0")
            self.term.insert("1.0", "... [Logs cleared down for memory optimization] ...\n", "heading")

        if "ERROR" in line or "!" in line:
            self.term.insert("end", line + "\n", "error")
        elif "=== Starting" in line or "=== Merge" in line or "=== Extract" in line:
            self.term.insert("end", line + "\n", "heading")
        elif "Done" in line or "Merged" in line or "Success" in line:
            self.term.insert("end", line + "\n", "success")
        else:
            self.term.insert("end", line + "\n")
        self.term.see("end")
        self.term.configure(state="disabled")

    def open_import(self):
        if self.running: return
        ImportWindow(self, self.folders, self.parent_folder, self._set_folders)

    def open_params(self):
        if self.running: return
        ParamsWindow(self, self.params, self._set_params)

    def _set_params(self, params):
        self.params.update(params)
        self.log("Parameters updated. Output folder: " + str(self.params.get("RESULTS_ROOT")))
        if not self.running:
            self._refresh_ready_status()

    def _set_folders(self, folders):
        self.folders = folders
        self.log(f"Imported {len(folders)} folder(s): " + ", ".join(os.path.basename(f) for f in folders))
        self.btn.configure(state="normal")
        self._update_queue_table()
        self._refresh_ready_status()

    def _refresh_ready_status(self):
        """Update the status message after an import: remind the user to choose the
        output folder in Parameters if it is not set yet."""
        n = len(self.folders)
        if not n:
            self.status.configure(text="Import folders to start")
            return
        out = self.params.get("RESULTS_ROOT", "")
        if not out:
            self.status.configure(text="⚙️ Set output folder in Parameters")
        else:
            self.status.configure(text=f"📊 {n} folder(s) ready")

    def _drain_queue(self):
        try:
            while True:
                item = self.msg_q.get_nowait()
                if item is None:
                    self.running = False
                    self.btn.configure(state="normal" if self.folders else "disabled")
                    self.btn_stop.configure(state="disabled")
                    self.status.configure(text="✅ Done — 100%")
                    
                    all_batches = [os.path.basename(f) for f in self.folders]
                    self._update_queue_table(finished_batches=all_batches)
                    
                    self.progress_bar.stop()
                    self.progress_bar.configure(mode="determinate", value=100 if not self.stop_requested else 0)
                    if self.stop_requested:
                        self.status.configure(text="🛑 Stopped — 0%")
                        self._update_queue_table()
                    else:
                        self.open_results_folder()
                    self.stop_requested = False
                else:
                    self.log(item)
                    
                    if self.running:
                        import re
                        
                        if "batch" in item.lower() and "samples" in item.lower():
                            try:
                                match_name = re.search(r"BATCH\s+(\w+)", item, re.IGNORECASE)
                                if match_name:
                                    current_b = match_name.group(1)
                                    if not hasattr(self, 'finished_list'):
                                        self.finished_list = []
                                    if current_b not in self.finished_list:
                                        self._update_queue_table(active_batch=current_b, finished_batches=self.finished_list)
                                
                                match_total = re.search(r"(\d+)\s+samples", item, re.IGNORECASE)
                                match_done = re.search(r"(\d+)\s+already\s+done", item, re.IGNORECASE)
                                
                                if match_total:
                                    self.total_samples = int(match_total.group(1))
                                    self.current_done = int(match_done.group(1)) if match_done else 0
                                    self.start_time = time.time()
                                    
                                    if self.total_samples > 0:
                                        val = (self.current_done / self.total_samples) * 100
                                        self.progress_bar.configure(mode="determinate", value=val)
                                        self.status.configure(text=f"⏳ Running ({self.current_done}/{self.total_samples}) — {val:.1f}%")
                                        continue
                            except Exception:
                                pass

                        if "raw=" in item and re.search(r"\[\d+s\]", item):
                            try:
                                if hasattr(self, 'total_samples') and hasattr(self, 'current_done'):
                                    self.current_done += 1
                                    if self.current_done > self.total_samples:
                                        self.current_done = self.total_samples
                                        
                                    val = (self.current_done / self.total_samples) * 100
                                    self.progress_bar.configure(mode="determinate", value=val)
                                    
                                    eta_text = ""
                                    if self.start_time and self.current_done > 0:
                                        elapsed = time.time() - self.start_time
                                        samples_processed_now = self.current_done - (self.total_samples - int(re.search(r"(\d+)\s+to\s+do", item, re.IGNORECASE).group(1) if "to do" in item else 0))
                                        if samples_processed_now <= 0:
                                            samples_processed_now = 1
                                        
                                        avg_time = elapsed / samples_processed_now
                                        remaining_samples = self.total_samples - self.current_done
                                        remaining_seconds = remaining_samples * avg_time
                                        
                                        if remaining_seconds > 0:
                                            mins, secs = divmod(int(remaining_seconds), 60)
                                            eta_text = f" — Rem.: ~{mins}m {secs}s"
                                            
                                    self.status.configure(text=f"⏳ Running ({self.current_done}/{self.total_samples}) — {val:.1f}%{eta_text}")
                                    continue
                            except Exception:
                                pass
                                
                        if "done." in item.lower():
                            try:
                                for f in self.folders:
                                    bname = os.path.basename(f)
                                    if bname not in self.finished_list:
                                        self.finished_list.append(bname)
                                        break
                            except Exception:
                                pass

                        if str(self.progress_bar.cget("mode")) == "determinate" and not hasattr(self, 'total_samples'):
                            self.progress_bar.configure(mode="indeterminate")
                            self.progress_bar.start(10)
                            
        except queue.Empty:
            pass
        self.after(80, self._drain_queue)

    def _existing_checkpoints(self):
        root = self.params.get("RESULTS_ROOT", "")
        found = []
        for f in self.folders:
            bname = os.path.basename(f)
            cp = os.path.join(root, bname, "checkpoint_traits.jsonl")
            if os.path.exists(cp) and os.path.getsize(cp) > 0:
                found.append((bname, cp))
        return found

    def _handle_existing_checkpoints(self, existing):
        names = ", ".join(b for b, _ in existing)
        win = tk.Toplevel(self)
        win.title("Existing checkpoint found")
        win.geometry("500x280")
        win.configure(bg="#f8f9fa")
        win.grab_set()
        msg = f"A checkpoint already exists for: {names}.\n\nSkip: samples already done will be SKIPPED.\n\nDelete: remove checkpoint(s) and recompute from scratch.\n\nCancel: go back."
        ttk.Label(win, text=msg, wraplength=460, justify="left", style="Info.TLabel").pack(fill="both", expand=True, padx=20, pady=20)
        choice = {"action": "cancel"}
        def pick(a):
            choice["action"] = a
            win.destroy()
        btns = ttk.Frame(win, padding=(20, 0, 20, 20))
        btns.pack(fill="x")
        ttk.Button(btns, text="Cancel", command=lambda: pick("cancel")).pack(side="right")
        ttk.Button(btns, text="Delete", command=lambda: pick("delete"), style="Danger.TButton").pack(side="right", padx=10)
        ttk.Button(btns, text="Skip", command=lambda: pick("continue"), style="Accent.TButton").pack(side="right", padx=10)
        self.wait_window(win)
        
        if choice["action"] == "cancel": return False
        if choice["action"] == "continue": return True
        if len(existing) == 1: to_delete = existing
        else:
            to_delete = self._pick_checkpoints_to_delete(existing)
            if to_delete is None: return False
        for bname, cp in to_delete:
            try:
                os.remove(cp)
                self.log(f"Deleted checkpoint for {bname}.")
            except Exception as e:
                self.log(f"! could not delete checkpoint for {bname}: {e}")
        return True

    def _pick_checkpoints_to_delete(self, existing):
        win = tk.Toplevel(self)
        win.title("Select checkpoints to delete")
        win.geometry("440x360")
        win.configure(bg="#f8f9fa")
        win.grab_set()
        ttk.Label(win, text="Choose which ones to delete.", style="Info.TLabel").pack(fill="x", padx=20, pady=(20, 10))
        frame = ttk.Frame(win, padding=10)
        frame.pack(fill="both", expand=True, padx=20)
        vars_ = []
        for bname, cp in existing:
            v = tk.BooleanVar(value=True)
            ttk.Checkbutton(frame, text=bname, variable=v).pack(anchor="w", pady=4)
            vars_.append(v)
        result = {"folders": None}
        def confirm():
            result["folders"] = [e for e, v in zip(existing, vars_) if v.get()]
            win.destroy()
        btns = ttk.Frame(win, padding=20)
        btns.pack(fill="x")
        ttk.Button(btns, text="Cancel", command=lambda: win.destroy()).pack(side="right")
        ttk.Button(btns, text="Delete selected", command=confirm, style="Danger.TButton").pack(side="right", padx=10)
        self.wait_window(win)
        return result["folders"]

    def start_analysis(self):
        if self.running or not self.folders: return
        if not self.params.get("RESULTS_ROOT", ""):
            messagebox.showwarning(
                "Output folder required",
                "Choose an output folder (results root) in Parameters before starting the analysis.")
            self.open_params()
            return
        existing = self._existing_checkpoints()
        if existing and not self._handle_existing_checkpoints(existing): return
        
        self.running = True
        self.stop_requested = False
        self.btn.configure(state="disabled")
        self.btn_stop.configure(state="normal")
        
        self.finished_list = []
        self.start_time = time.time()
        if hasattr(self, 'total_samples'): delattr(self, 'total_samples')
        
        self._update_queue_table()
        self.status.configure(text="⏳ Preparing analysis...")
        self.progress_bar.configure(mode="determinate", value=0)
        self.log("=== Starting analysis ===")
        threading.Thread(target=self._worker, daemon=True).start()

    def stop_analysis(self):
        if not self.running: return
        self.stop_requested = True
        self.btn_stop.configure(state="disabled")
        self.status.configure(text="🛑 Stopping...")
        self.log("=== Stop requested by user ===")

    def save_console_log(self):
        content = self.term.get("1.0", "end-1c")
        if not content.strip() or "Live Console Output" in content and len(content.strip()) < 30:
            messagebox.showinfo("Information", "The console is currently empty. Nothing to save.")
            return
        self.attributes("-topmost", False)
        file_path = filedialog.asksaveasfilename(title="Save Console Log As...", defaultextension=".txt", filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")], parent=self)
        if file_path:
            try:
                with open(file_path, "w", encoding="utf-8") as f: f.write(content)
                self.log(f"✅ Console log successfully saved to: {file_path}")
            except Exception as e: messagebox.showerror("Error", f"Could not save file: {e}")
        self.attributes("-topmost", True)

    def clear_session(self):
        if self.running:
            messagebox.showwarning("Analysis Running", "Cannot clear session while an analysis is in progress.")
            return
        if messagebox.askyesno("New Session", "Are you sure you want to clear the current queue and terminal? This will reset the workspace."):
            self.folders.clear()
            self.parent_folder = ""
            self._update_queue_table()
            self.term.configure(state="normal")
            self.term.delete("1.0", "end")
            self.term.configure(state="disabled")
            self.progress_bar.configure(mode="determinate", value=0)
            self.status.configure(text="Import folders to start")
            self.btn.configure(state="disabled")
            if os.path.exists(CONFIG_FILE):
                try: os.remove(CONFIG_FILE)
                except Exception: pass
            self.log("✨ New session started. Workspace cleared.")

    def _run_bg(self, fn):
        if self.running: return
        self.running = True
        threading.Thread(target=self._tool_worker, args=(fn,), daemon=True).start()

    def _tool_worker(self, fn):
        import contextlib
        writer = api._LineWriter(lambda s: self.msg_q.put(s))
        try:
            with contextlib.redirect_stdout(writer): fn(self.fmt_var.get())
            writer.flush()
        except Exception as e: self.msg_q.put(f"ERROR: {e}")
        finally: self.msg_q.put(None)

    def tool_report(self):
        # No folder to pick: the report is always built on the project's results folder.
        from tools import figure_report as grf
        self.log(f"=== Generate figure report ({grf.RESULTS_ROOT}) ===")
        def do(fmt):
            # rank the samples first, so that the report lists the most suspect first
            from tools import qc_rank
            qc_rank.run(grf.RESULTS_ROOT)
            path = grf.generate()                    # uses the project's RESULTS_ROOT
            if path:
                try:
                    webbrowser.open("file://" + os.path.abspath(path))
                    print("Opened in your browser.")
                except Exception as e:
                    print(f"Report ready ({path}); could not auto-open: {e}")
        self._run_bg(do)

    def tool_merge(self):
        d = filedialog.askdirectory(title="Select the results folder to merge")
        if not d: return
        self.log(f"=== Merge batches in {d} ===")
        from tools import merge_batches as mb
        def do(fmt):
            frames, presence = [], {}
            paths = sorted(__import__("glob").glob(__import__("os").path.join(d, "*", "*.xlsx")))
            paths = [p for p in paths if "merged_traits" not in __import__("os").path.basename(p)]
            if not paths: print(f"No trait table found."); return
            import pandas as pd, os as _os
            for p in paths:
                batch = _os.path.basename(_os.path.dirname(p))
                df = mb.read_batch_table(p)
                if df is None or len(df) == 0: continue
                df.insert(0, "batch", batch); frames.append(df)
                for c in df.columns: presence[c] = presence.get(c, 0) + 1
            if not frames: print("Nothing to merge."); return
            merged = pd.concat(frames, ignore_index=True, sort=False)
            lead = [c for c in ["batch","ID","n_raw","n_removed","%removed","pivot_return","n_rescued",
                                "collar_raise","hypocotyl","time"] if c in merged.columns]
            merged = merged[lead + [c for c in merged.columns if c not in lead]]
            base = _os.path.join(d, "merged_traits")
            if fmt in ("xlsx","both"): merged.to_excel(base+".xlsx", index=False); print(f"  -> {base}.xlsx")
            if fmt in ("csv","both"):  merged.to_csv(base+".csv", index=False); print(f"  -> {base}.csv")
        self._run_bg(do)

    def tool_checkpoint(self):
        d = filedialog.askdirectory(title="Select the results folder (with checkpoints)")
        if not d: return
        self.log(f"=== Extract checkpoints in {d} ===")
        from tools import extract_checkpoint as ec
        import glob, os as _os
        def do(fmt):
            paths = sorted(glob.glob(os.path.join(d, "*", "*checkpoint*.jsonl"))) + sorted(glob.glob(os.path.join(d, "*checkpoint*.jsonl")))
            if not paths: print("No checkpoint file found."); return
            for pth in paths:
                df = ec.read_checkpoint(pth)
                if df is None: continue
                batch = _os.path.basename(_os.path.dirname(pth)) or "checkpoint"
                print(f"{batch}: {len(df)} samples done so far")
                base = _os.path.join(_os.path.dirname(pth), f"traits_{batch}_partial")
                ec.write(df, base, fmt)
        self._run_bg(do)

    def tool_delete_checkpoints(self):
        d = filedialog.askdirectory(title="Select the results folder")
        if not d: return
        import glob
        paths = sorted(glob.glob(os.path.join(d, "*", "*checkpoint*.jsonl"))) + sorted(glob.glob(os.path.join(d, "*checkpoint*.jsonl")))
        if not paths: self.log(f"No checkpoint found in {d}."); return
        listing = "\n".join(f"  - {os.path.relpath(p, d)}" for p in paths)
        if not messagebox.askyesno("Delete?", f"{len(paths)} checkpoints will be deleted:\n\n{listing}"): return
        for p in paths:
            try: os.remove(p); self.log(f"  deleted: {os.path.relpath(p, d)}")
            except Exception as e: self.log(f"  ! could not delete {p}: {e}")

    def _worker(self):
        try: api.run(self.folders, params=self.params, progress=lambda s: self.msg_q.put(s), should_stop=lambda: self.stop_requested)
        except Exception as e: self.msg_q.put(f"ERROR: {e}")
        finally: self.msg_q.put(None)

if __name__ == "__main__":
    import multiprocessing as mp
    mp.freeze_support()
    App().mainloop()