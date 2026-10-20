"""
Code to interact with asimov projects.
"""

import os
import yaml


def get_results_file(project, event, analysis):
    """
    Find the PESummary metafile for a given event and analysis.
    """

    return os.path.join(project, "results", event, analysis, "posterior_samples.h5")


def _get_nested(obj, key):
    """Resolve a dot-notation key against a nested dict.  Returns the value, or
    ``None`` if any intermediate key is missing."""
    for part in key.split("."):
        if not isinstance(obj, dict) or part not in obj:
            return None
        obj = obj[part]
    return obj


def _production_matches(prod_data, conditions):
    """Return True if *prod_data* satisfies every condition in *conditions*.

    Each condition is a single-key dict.  Keys support dot notation
    (e.g. ``waveform.approximant``).  For the bare ``status`` key the
    production-level value is checked first; if that doesn't match,
    ``postprocessing.status`` is tried as a fallback.
    """
    for condition in conditions:
        for key, expected in condition.items():
            value = _get_nested(prod_data, key)
            # Special-case: bare 'status' may live under postprocessing
            if value != expected and key == "status":
                value = _get_nested(prod_data, "postprocessing.status")
            if value != expected:
                return False
    return True


def _is_approved(prod_data):
    """Return True if the most recent review of the production is APPROVED.

    Reviews are appended to, so a later REJECTED entry supersedes an earlier
    APPROVED one.
    """
    reviews = prod_data.get("review") or []
    return bool(reviews) and reviews[-1].get("status") == "APPROVED"


def get_approved_event_data(project, ledger_path, subjects=None, analyses_criteria=None):
    """Extract event data from an asimov ledger suitable for SmallMultiples.plot().

    For each event in *subjects* (or all events when ``None``), the last
    production that satisfies *analyses_criteria* **and** carries an APPROVED
    review is selected.  The metafile path is constructed from the standard
    asimov results layout.

    Parameters
    ----------
    project:
        Root directory of the asimov project.
    ledger_path:
        Path to the ledger YAML file.
    subjects:
        Optional list of event names to include.  All events are used when
        ``None``.
    analyses_criteria:
        List of AND-groups (each a list of ``{key: value}`` dicts).  A
        production matches if it satisfies any one AND-group.  ``None``
        means no filtering beyond the APPROVED requirement.

    Returns
    -------
    list of dict
        Each dict has ``name`` and ``metafile`` keys, ready for
        ``SmallMultiples.plot()``.
    """
    with open(ledger_path) as f:
        ledger = yaml.safe_load(f)

    events_out = []

    for event in ledger.get("events", []):
        name = event["name"]

        if subjects is not None and name not in subjects:
            continue

        best = None
        for prod_entry in event.get("productions", []):
            for prod_name, prod_data in prod_entry.items():
                if not _is_approved(prod_data):
                    continue
                if analyses_criteria is not None:
                    # OR of AND-groups
                    if not any(
                        _production_matches(prod_data, group)
                        for group in analyses_criteria
                    ):
                        continue
                best = prod_name  # keep last match — highest run number wins

        if best is None:
            continue

        events_out.append(
            {"name": name, "metafile": get_results_file(project, name, best)}
        )

    return events_out
    
