#!/usr/bin/env python3

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal

import pandas as pd


META_STOP_COL = "height"  # sample columns start after this in filtered_for_mzmine.csv


def _parse_float_series(series: pd.Series) -> pd.Series:
    """Parse numeric columns that may contain 'NA', empty, or decimal commas."""
    s = series.astype(str).str.strip()
    s = s.replace({"": pd.NA, "NA": pd.NA, "NaN": pd.NA, "nan": pd.NA})
    s = s.str.replace(",", ".", regex=False)
    return pd.to_numeric(s, errors="coerce")


def _read_filtered(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)

    required = ["mz_range.min", "mz_range.max", "id"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in {path}: {missing}")

    df["mz_range.min"] = _parse_float_series(df["mz_range.min"])
    df["mz_range.max"] = _parse_float_series(df["mz_range.max"])
    df["id"] = df["id"].astype(str)

    if META_STOP_COL not in df.columns:
        raise ValueError(
            f"Expected column '{META_STOP_COL}' in {path} to detect sample columns. "
            f"Columns are: {list(df.columns)[:15]} ..."
        )

    stop_idx = df.columns.get_loc(META_STOP_COL)
    sample_cols = list(df.columns[stop_idx + 1 :])
    if not sample_cols:
        raise ValueError(f"No sample columns found after '{META_STOP_COL}' in {path}")

    for col in sample_cols:
        df[col] = _parse_float_series(df[col])

    return df


def _read_pools(path: Path) -> pd.DataFrame:
    # pools.csv is semicolon-separated and uses decimal commas
    df = pd.read_csv(path, sep=";", engine="python")

    required = ["Name", "M+H"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in {path}: {missing}")

    df = df[["Name", "M+H"]].copy()
    df["Name"] = df["Name"].astype(str).str.strip()
    df["M+H"] = _parse_float_series(df["M+H"])
    df = df.dropna(subset=["M+H"]).reset_index(drop=True)

    if df.empty:
        raise ValueError(f"No numeric M+H values parsed from {path}")

    return df


@dataclass(frozen=True)
class Match:
    feature_idx: int
    substrate: str
    pool_mh: float


def _match_features_to_pools(
    filtered: pd.DataFrame,
    pools: pd.DataFrame,
    tolerance: float,
) -> list[Match]:
    matches: list[Match] = []

    mh_values = pools["M+H"].to_numpy()
    names = pools["Name"].to_list()

    mz_min = filtered["mz_range.min"].to_numpy()
    mz_max = filtered["mz_range.max"].to_numpy()

    for i in range(len(filtered)):
        lo = mz_min[i]
        hi = mz_max[i]
        if pd.isna(lo) or pd.isna(hi):
            continue

        lo2 = lo - tolerance
        hi2 = hi + tolerance

        for substrate, mh in zip(names, mh_values):
            if lo2 <= mh <= hi2:
                matches.append(Match(feature_idx=i, substrate=substrate, pool_mh=float(mh)))

    return matches


def _aggregate(values: pd.Series, how: Literal["max", "sum", "mean"]) -> float:
    values = values.dropna()
    if values.empty:
        return float("nan")
    if how == "max":
        return float(values.max())
    if how == "sum":
        return float(values.sum())
    if how == "mean":
        return float(values.mean())
    raise ValueError(f"Unknown aggregation: {how}")


def build_area_matrix(
    filtered: pd.DataFrame,
    pools: pd.DataFrame,
    tolerance: float,
    aggregation: Literal["max", "sum", "mean"] = "max",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    stop_idx = filtered.columns.get_loc(META_STOP_COL)
    sample_cols = list(filtered.columns[stop_idx + 1 :])

    matches = _match_features_to_pools(filtered=filtered, pools=pools, tolerance=tolerance)
    if not matches:
        empty = pd.DataFrame(index=sample_cols)
        return empty, pd.DataFrame()

    # Long form: one row per (sample, substrate, matched feature)
    records: list[dict] = []
    for m in matches:
        row = filtered.iloc[m.feature_idx]
        feature_id = row.get("id")
        lo = row.get("mz_range.min")
        hi = row.get("mz_range.max")

        for sample in sample_cols:
            records.append(
                {
                    "sample": sample,
                    "substrate": m.substrate,
                    "area": row[sample],
                    "feature_id": feature_id,
                    "pool_mh": m.pool_mh,
                    "mz_range_min": lo,
                    "mz_range_max": hi,
                }
            )

    long_df = pd.DataFrame.from_records(records)

    # Aggregate (in case multiple features match a substrate for a sample)
    grouped = (
        long_df.groupby(["sample", "substrate"], as_index=False)["area"]
        .apply(lambda s: _aggregate(s, aggregation))
        .rename(columns={"area": "area"})
    )

    matrix = grouped.pivot(index="sample", columns="substrate", values="area").sort_index()
    matrix = matrix.sort_index(axis=1)

    return matrix, long_df


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Match filtered mzMine feature table peaks to pool substrates by M+H. "
            "A pool entry matches a feature if pool M+H is within [mz_range.min, mz_range.max] ± tolerance. "
            "Outputs a sample×substrate area matrix."
        )
    )
    parser.add_argument(
        "--filtered",
        type=Path,
        default=Path("filtered_for_mzmine.csv"),
        help="Path to filtered_for_mzmine.csv (comma-separated).",
    )
    parser.add_argument(
        "--pools",
        type=Path,
        default=Path("pools.csv"),
        help="Path to pools.csv (semicolon-separated, decimal commas).",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=0.0125,
        help="Allowed absolute deviation (in m/z units). Default: 0.0125",
    )
    parser.add_argument(
        "--aggregation",
        choices=["max", "sum", "mean"],
        default="max",
        help="How to combine multiple matching features per (sample, substrate).",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("matched_pool_areas.csv"),
        help="Output CSV path for the sample×substrate matrix.",
    )
    parser.add_argument(
        "--out-long",
        type=Path,
        default=Path("matched_pool_areas_long.csv"),
        help="Output CSV path for the long match table (traceability).",
    )

    args = parser.parse_args()

    filtered = _read_filtered(args.filtered)
    pools = _read_pools(args.pools)

    matrix, long_df = build_area_matrix(
        filtered=filtered,
        pools=pools,
        tolerance=args.tolerance,
        aggregation=args.aggregation,
    )

    matrix.to_csv(args.out, index=True)
    if not long_df.empty:
        long_df.to_csv(args.out_long, index=False)

    print(f"Wrote matrix: {args.out} (shape={matrix.shape})")
    if not long_df.empty:
        print(f"Wrote matches: {args.out_long} (rows={len(long_df)})")
    else:
        print("No matches found; wrote empty matrix.")


if __name__ == "__main__":
    main()
