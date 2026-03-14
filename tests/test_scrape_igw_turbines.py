import json

import pandas as pd
import pytest

from scow.scrape_igw_turbines import (
    LEGACY_CACHE_NAME,
    download_and_process_igw_legacy_schema,
    validate_legacy_igw_frame,
)


class _FakeResponse:
    def __init__(self, text):
        self.text = text

    def raise_for_status(self):
        return None


def test_download_and_process_igw_legacy_schema_writes_canonical_cache(monkeypatch, tmp_path):
    wheels = {
        "1": [
            {
                "project": "Pretul",
                "x": "15.0",
                "y": "47.0",
                "facilityInfo": "2 Anlage(n), je 3.2 MW, errichtet: 2014",
                "typeInfo": "Type: ENERCON, E115/3.000<br>Nabenhöhe: 149, Rotordurchmesser: 115",
                "mweinzeln": "3.2",
            }
        ]
    }
    html = f"<div wheels='{json.dumps(wheels)}'></div>"

    monkeypatch.setattr(
        "scow.scrape_igw_turbines.requests.get",
        lambda *args, **kwargs: _FakeResponse(html),
    )

    result = download_and_process_igw_legacy_schema(tmp_path)
    cache_path = tmp_path / LEGACY_CACHE_NAME

    assert cache_path.exists()
    cached = pd.read_csv(cache_path, sep=";", decimal=",")
    assert list(cached.columns) == [
        "Name",
        "Betreiber1",
        "Betreiber2",
        "n_Anlagen",
        "kW",
        "Type",
        "Jahr",
        "x",
        "lat",
        "lon",
        "url",
        "Hersteller",
        "Nabenhöhe",
        "Rotordurchmesser",
    ]
    assert result.loc[0, "Name"] == "Pretul"
    assert result.loc[0, "Type"] == "E82 E4"
    assert result.loc[0, "Betreiber1"] == "Österreichische Bundesforste"


def test_validate_legacy_igw_frame_rejects_missing_coordinates():
    df = pd.DataFrame(
        {
            "Name": ["Broken"],
            "Betreiber1": [""],
            "Betreiber2": [""],
            "n_Anlagen": [1],
            "kW": [1000],
            "Type": ["Foo"],
            "Jahr": [2010],
            "x": [15.0],
            "lat": [47.0],
            "lon": [None],
            "url": [""],
            "Hersteller": ["Bar"],
            "Nabenhöhe": [80.0],
            "Rotordurchmesser": [90.0],
        }
    )

    with pytest.raises(ValueError, match="finite x, lon, and lat"):
        validate_legacy_igw_frame(df)
