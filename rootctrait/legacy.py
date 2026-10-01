"""Former output column names, kept only to read files written by earlier versions.

Versions before 2.0.0 wrote some table columns and checkpoint keys with French
names. They are renamed to their English equivalents when an old checkpoint or
trait table is read, so that partial runs can be resumed and old tables merged.
"""

# former table column / trait name -> current name
LEGACY_COLUMNS = {
    "n_brut": "n_raw",
    "n_retire": "n_removed",
    "%retire": "%removed",
    "NRL_court_<5": "NRL_short_<5",
    "NRL_moyen_5_15": "NRL_medium_5_15",
}

# former checkpoint record key -> current key
LEGACY_RECORD_KEYS = {"n_brut": "n_raw", "n_ret": "n_removed"}


def upgrade_record(rec):
    """Return a checkpoint record with current key and trait names."""
    out = {LEGACY_RECORD_KEYS.get(k, k): v for k, v in rec.items()}
    T = out.get("T")
    if isinstance(T, dict):
        out["T"] = {LEGACY_COLUMNS.get(k, k): v for k, v in T.items()}
    return out


def upgrade_columns(df):
    """Rename former column names of a pandas DataFrame in place and return it."""
    df.rename(columns={k: v for k, v in LEGACY_COLUMNS.items() if k in df.columns}, inplace=True)
    return df
