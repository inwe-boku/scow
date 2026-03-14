#!/usr/bin/env python3
from __future__ import annotations

import json
import logging
import re
from html import unescape as html_unescape
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup


IGW_PAGE_URL = "https://www.igwindkraft.at/aktuelles/windrad-karte"
LEGACY_CACHE_NAME = "igwturbines.csv"


def _to_float_de(value: str | None) -> float:
    if value is None:
        return np.nan
    text = str(value).strip()
    if not text:
        return np.nan
    text = text.replace(" ", "")
    match = re.search(r"[-+]?\d+(?:[.,]\d+)?", text)
    if not match:
        return np.nan
    token = match.group(0).replace(",", ".")
    try:
        return float(token)
    except ValueError:
        return np.nan


def extract_wheels_json(html: str) -> dict:
    """Extract the IGW `wheels=...` JSON payload from the page HTML."""
    soup = BeautifulSoup(html, "html.parser")

    node = soup.select_one("[wheels]")
    if node is not None:
        wheels_raw = node.get("wheels")
        if wheels_raw:
            return json.loads(wheels_raw)

    match = re.search(r"wheels\s*=\s*'(.+?)'", html, flags=re.DOTALL)
    if not match:
        match = re.search(r'wheels\s*=\s*"(.+?)"', html, flags=re.DOTALL)
    if not match:
        raise RuntimeError("Could not locate wheels=... in HTML.")

    return json.loads(html_unescape(match.group(1)))


def parse_typeinfo(typeinfo_html: str | None) -> tuple[str, str, float, float]:
    """Parse the IGW turbine type HTML into manufacturer, model, hub height, and rotor diameter."""
    raw = typeinfo_html or ""
    norm = raw.replace("<br>", "\n").replace("<br/>", "\n").replace("<br />", "\n")
    norm = BeautifulSoup(norm, "html.parser").get_text(separator="\n", strip=True)

    manufacturer = ""
    model = ""
    hub_height = np.nan
    rotor_diameter = np.nan

    match = re.search(r"(?im)^\s*Type\s*:\s*(.+?)\s*(?:,\s*(.+?))?\s*$", norm)
    if match:
        manufacturer = (match.group(1) or "").strip()
        model = (match.group(2) or "").strip()
    else:
        first_line = norm.splitlines()[0].strip() if norm.strip() else ""
        fallback = re.search(r"(?i)type\s*:\s*(.+)", first_line)
        payload = fallback.group(1).strip() if fallback else first_line
        if "," in payload:
            manufacturer, model = [part.strip() for part in payload.split(",", 1)]
        else:
            manufacturer = payload.strip()

    hub_match = re.search(r"(?i)\bNabenhöhe\s*:\s*([0-9]+(?:[.,][0-9]+)?)", norm)
    if hub_match:
        hub_height = _to_float_de(hub_match.group(1))

    rotor_match = re.search(r"(?i)\bRotordurchmesser\s*:\s*([0-9]+(?:[.,][0-9]+)?)", norm)
    if rotor_match:
        rotor_diameter = _to_float_de(rotor_match.group(1))

    return manufacturer, model, hub_height, rotor_diameter


def parse_facilityinfo(facilityinfo: str | None) -> tuple[float, float, float]:
    """Parse the IGW facility string into turbine count, MW per unit, and build year."""
    if not facilityinfo:
        return np.nan, np.nan, np.nan

    text = str(facilityinfo)

    count_match = re.search(r"(?i)\b(\d+)\s*Anlage\(n\)\b", text)
    n_anlagen = float(int(count_match.group(1))) if count_match else np.nan

    mw_match = re.search(r"(?i)\bje\s*([0-9]+(?:[.,][0-9]+)?)\s*MW\b", text)
    mw_each = _to_float_de(mw_match.group(1)) if mw_match else np.nan

    year_match = re.search(r"(?i)\berrichtet\s*:\s*(\d{4})\b", text)
    built_year = float(int(year_match.group(1))) if year_match else np.nan

    return n_anlagen, mw_each, built_year


def normalize_legacy_igw_numeric_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize numeric columns in the canonical legacy IGW cache schema."""
    df = df.copy()
    numeric_cols = [
        "n_Anlagen",
        "kW",
        "Jahr",
        "x",
        "lat",
        "lon",
        "Nabenhöhe",
        "Rotordurchmesser",
    ]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col].astype(str).str.replace(",", ".", regex=False), errors="coerce")
    return df


def _request_headers() -> dict[str, str]:
    return {
        "User-Agent": (
            "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
            "AppleWebKit/605.1.15 (KHTML, like Gecko) Version/17.0 Safari/605.1.15"
        ),
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
        "Accept-Language": "de-AT,de;q=0.9,en;q=0.8",
    }


def _legacy_frame_from_wheels(wheels: dict) -> tuple[pd.DataFrame, list[dict]]:
    rows = []
    raw_rows = []
    for group_id, items in wheels.items():
        for idx, item in enumerate(items):
            data = {k: v for k, v in dict(item).items() if not str(k).isdigit()}
            raw_rows.append({"group_id": str(group_id), "idx_in_group": int(idx), **data})

            project = data.get("project", "")
            lon = data.get("x", np.nan)
            lat = data.get("y", np.nan)
            facility_info = data.get("facilityInfo", "")
            type_info = data.get("typeInfo", "")

            n_anlagen, mw_each, built_year = parse_facilityinfo(facility_info)
            manufacturer, model, hub_h, rotor_d = parse_typeinfo(type_info)

            mw_each = _to_float_de(data.get("mweinzeln", mw_each))
            kw_each = int(round(mw_each * 1000.0)) if np.isfinite(mw_each) else np.nan

            rows.append(
                {
                    "Name": project,
                    "Betreiber1": "",
                    "Betreiber2": "",
                    "n_Anlagen": n_anlagen,
                    "kW": kw_each,
                    "Type": model,
                    "Jahr": built_year,
                    "x": lon,
                    "lat": lat,
                    "lon": lon,
                    "url": "",
                    "Hersteller": manufacturer,
                    "Nabenhöhe": hub_h,
                    "Rotordurchmesser": rotor_d,
                }
            )

    igw = pd.DataFrame(
        rows,
        columns=[
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
        ],
    )

    type_replacements = {
        "E-40": "E40",
        "E40/5.40": "E40 5.40",
        "E40 5.4": "E40 5.40",
        "E66 18.7": "E66 18.70",
        "E66/18.70": "E66 18.70",
        "E66.18": "E66 18.70",
        "E66 20.7": "E66 20.70",
        "E70/E4": "E70 E4",
        "E70/20.71": "E70 E4",
        "E70": "E70 E4",
        "E-101": "E101",
        "E 101": "E101",
        "E115/3.000": "E115",
        "3.XM": "3XM",
        "V126/3450": "V126",
    }
    manufacturer_replacements = {"ENERCON": "Enercon", "DeWind": "Dewind"}

    igw["Type"] = igw["Type"].replace(type_replacements)
    igw["Hersteller"] = igw["Hersteller"].replace(manufacturer_replacements)

    igw.loc[(igw["Type"] == "E40") & (igw["kW"] == 500), "Type"] = "E40 5.40"
    igw.loc[(igw["Type"] == "E40") & (igw["kW"] == 600), "Type"] = "E40 6.44"
    igw.loc[(igw["Type"] == "E66") & (igw["kW"] == 1800), "Type"] = "E66 18.70"
    igw.loc[(igw["Type"] == "E82") & (igw["kW"] == 2300), "Type"] = "E82 E2"
    igw.loc[(igw["Type"] == "E115") & (igw["kW"] == 3200), "Type"] = "E115 E2"
    igw.loc[(igw["Type"] == "M114") & (igw["kW"] == 3170), "Type"] = "3.2M114"

    oberwaltersdorf = igw["Name"].astype(str).str.contains("Oberwaltersdorf", na=False)
    igw.loc[oberwaltersdorf, ["Type", "Nabenhöhe", "Rotordurchmesser"]] = ["V112", 140, 112]

    pretul = igw["Name"].astype(str).str.contains("Pretul", na=False)
    igw.loc[pretul, ["Type", "Betreiber1"]] = ["E82 E4", "Österreichische Bundesforste"]

    igw[["Nabenhöhe", "Rotordurchmesser"]] = (
        igw[["Nabenhöhe", "Rotordurchmesser"]].replace("", np.nan).astype(float)
    )
    return normalize_legacy_igw_numeric_columns(igw), raw_rows


def download_and_process_igw_legacy_schema(turbines_dir: Path) -> pd.DataFrame:
    """Download the current IGW map and write the canonical legacy cache schema."""
    turbines_dir.mkdir(parents=True, exist_ok=True)

    response = requests.get(IGW_PAGE_URL, headers=_request_headers(), timeout=60)
    response.raise_for_status()
    wheels = extract_wheels_json(response.text)

    igw, raw_rows = _legacy_frame_from_wheels(wheels)
    (turbines_dir / "turbines_raw.json").write_text(
        json.dumps(raw_rows, ensure_ascii=False),
        encoding="utf-8",
    )
    cache_path = turbines_dir / LEGACY_CACHE_NAME
    igw.to_csv(cache_path, sep=";", decimal=",", encoding="utf-8", index=False)
    logging.info("Wrote %s rows to %s", len(igw), cache_path)
    return igw


if __name__ == "__main__":
    from config import repo

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    out_dir = Path(repo) / "data" / "site" / "AUT_turbines"
    df = download_and_process_igw_legacy_schema(out_dir)
    print(f"Rows: {len(df)}")
    print(f"Wrote: {out_dir / LEGACY_CACHE_NAME}")
