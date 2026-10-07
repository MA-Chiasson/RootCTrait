"""Trait names. Version 2.5.0 renamed the traits whose abbreviation came from French
to English abbreviations. OLD_TO_NEW maps the names written by RootCTrait up to 2.4.x to
the current names; english_names() converts a trait table or a dict of traits written
by an older version, so that older results can be read with the current tools."""

OLD_TO_NEW = {
    'LRP': 'PRL',                     # primary root length
    'NRL': 'NLR',                     # number of lateral roots
    'NRL_short_<5': 'NLR_short_<5',
    'NRL_medium_5_15': 'NLR_medium_5_15',
    'NRL_long_>15': 'NLR_long_>15',
    'PM': 'MD',                       # maximum depth
    'LM': 'MW',                       # maximum width
    'RLP': 'WDR',                     # width to depth ratio
    'VRT': 'TRV',                     # total root volume
    'SRT': 'TRSA',                    # total root surface area
    'IC': 'CI',                       # compactness index
    'SRL': 'RLV',                     # root length per root volume
    'DR': 'LRD',                      # lateral root density
    'NTR': 'NCR',                     # number of roots emerging near the collar
    'DRP': 'PRD',                     # primary root diameter
    'DRS': 'MLD',                     # mean lateral root diameter
}


def english_names(obj):
    """Rename old trait names (keys of a dict, or columns of a pandas DataFrame)."""
    if hasattr(obj, 'rename') and hasattr(obj, 'columns'):
        return obj.rename(columns={k: v for k, v in OLD_TO_NEW.items() if k in obj.columns})
    if isinstance(obj, dict):
        return {OLD_TO_NEW.get(k, k): v for k, v in obj.items()}
    return obj
