"""Tests for NanoComp --split_runs with input lacking run IDs (e.g. fastq)."""
import pandas as pd
import pytest

from nanocomp.utils import change_identifiers


def test_split_runs_without_runids_exits_cleanly():
    """Non-summary input has no runIDs column: must exit cleanly, not KeyError."""
    df = pd.DataFrame({"lengths": [100, 200], "dataset": ["a", "a"]})
    with pytest.raises(SystemExit):
        change_identifiers(df, {"run1": "runA"})


def test_split_runs_with_runids_applies_names():
    """Summary input keeps working: run IDs are renamed to the given names."""
    df = pd.DataFrame(
        {
            "lengths": [100, 200],
            "dataset": ["x", "x"],
            "runIDs": ["run1", "run1"],
        }
    )
    change_identifiers(df, {"run1": "runA"})
    assert list(df["dataset"]) == ["runA", "runA"]
