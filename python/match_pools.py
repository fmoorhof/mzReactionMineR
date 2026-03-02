#!/usr/bin/env python3

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import pandas as pd


META_STOP_COL = "height"  # sample columns start after this in filtered_for_mzmine.csv
DEFAULT_RT_COL = "rt"
DEFAULT_POOLS_RT_COL = "Retention time"
DEFAULT_POOLS_POOL_COL = "Pool"


def _parse_float_series(series: pd.Series) -> pd.Series:
    """Parse numeric columns that may contain 'NA', empty, or decimal commas."""
    s = series.astype(str).str.strip()
    s = s.replace({"": pd.NA, "NA": pd.NA, "NaN": pd.NA, "nan": pd.NA})
    s = s.str.replace(",", ".", regex=False)
    return pd.to_numeric(s, errors="coerce")


def _unique_preserve_order(values: list[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for v in values:
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out


def _read_filtered(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)

    required = ["mz_range.min", "mz_range.max", "id", DEFAULT_RT_COL]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in {path}: {missing}")

    df["mz_range.min"] = _parse_float_series(df["mz_range.min"])
    df["mz_range.max"] = _parse_float_series(df["mz_range.max"])
    df[DEFAULT_RT_COL] = _parse_float_series(df[DEFAULT_RT_COL])
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

    keep_cols = ["Name", "M+H"]
    if DEFAULT_POOLS_RT_COL in df.columns:
        keep_cols.append(DEFAULT_POOLS_RT_COL)
    if DEFAULT_POOLS_POOL_COL in df.columns:
        keep_cols.append(DEFAULT_POOLS_POOL_COL)

    df = df[keep_cols].copy()
    df["Name"] = df["Name"].astype(str).str.strip()
    df["M+H"] = _parse_float_series(df["M+H"])
    if DEFAULT_POOLS_RT_COL in df.columns:
        df[DEFAULT_POOLS_RT_COL] = _parse_float_series(df[DEFAULT_POOLS_RT_COL])
    if DEFAULT_POOLS_POOL_COL in df.columns:
        df[DEFAULT_POOLS_POOL_COL] = df[DEFAULT_POOLS_POOL_COL].astype(str).str.strip()
    df = df.dropna(subset=["M+H"]).reset_index(drop=True)

    if df.empty:
        raise ValueError(f"No numeric M+H values parsed from {path}")

    return df


@dataclass(frozen=True)
class Match:
    feature_idx: int
    substrate: str
    pool_mh: float
    pool_rt: float
    feature_rt: float
    delta_rt: float
    pool: str | None


def _match_features_to_pools(
    filtered: pd.DataFrame,
    pools: pd.DataFrame,
    tolerance: float,
    rt_tol: float,
    rt_col: str = DEFAULT_RT_COL,
    pools_rt_col: str = DEFAULT_POOLS_RT_COL,
) -> list[Match]:
    if rt_tol <= 0:
        raise ValueError("rt_tol must be > 0")
    if rt_col not in filtered.columns:
        raise ValueError(f"Filtered table is missing RT column '{rt_col}'")
    if pools_rt_col not in pools.columns:
        raise ValueError(
            f"Pools table is missing RT column '{pools_rt_col}'. "
            f"Either add it to pools.csv or switch RT matching off (not supported in this mode)."
        )

    matches: list[Match] = []

    mh_values = pools["M+H"].to_numpy()
    names = pools["Name"].to_list()
    rt_values = pools[pools_rt_col].to_numpy()
    pools_values = pools[DEFAULT_POOLS_POOL_COL].to_list() if DEFAULT_POOLS_POOL_COL in pools.columns else [None] * len(pools)

    mz_min = filtered["mz_range.min"].to_numpy()
    mz_max = filtered["mz_range.max"].to_numpy()
    feature_rt = filtered[rt_col].to_numpy()

    for i in range(len(filtered)):
        lo = mz_min[i]
        hi = mz_max[i]
        rt_f = feature_rt[i]
        if pd.isna(lo) or pd.isna(hi) or pd.isna(rt_f):
            continue

        lo2 = lo - tolerance
        hi2 = hi + tolerance

        best_idx: int | None = None
        best_delta_rt: float | None = None

        for j, (substrate, mh, rt_p, pool_code) in enumerate(zip(names, mh_values, rt_values, pools_values)):
            if pd.isna(rt_p):
                continue  # hard RT filter: cannot match without pool RT
            if not (lo2 <= mh <= hi2):
                continue
            delta = abs(float(rt_f) - float(rt_p))
            if delta > rt_tol:
                continue

            if best_idx is None or delta < (best_delta_rt or float("inf")):
                best_idx = j
                best_delta_rt = delta

        if best_idx is not None:
            matches.append(
                Match(
                    feature_idx=i,
                    substrate=names[best_idx],
                    pool_mh=float(mh_values[best_idx]),
                    pool_rt=float(rt_values[best_idx]),
                    feature_rt=float(rt_f),
                    delta_rt=float(best_delta_rt),
                    pool=pools_values[best_idx],
                )
            )

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
    rt_tol: float,
    aggregation: Literal["max", "sum", "mean"] = "max",
    include_all_substrates: bool = True,
    substrate_order: Literal["pools-file", "alphabetical", "custom"] = "pools-file",
    substrate_order_file: Path | None = None,
    rt_col: str = DEFAULT_RT_COL,
    pools_rt_col: str = DEFAULT_POOLS_RT_COL,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    stop_idx = filtered.columns.get_loc(META_STOP_COL)
    sample_cols = list(filtered.columns[stop_idx + 1 :])

    if substrate_order == "custom":
        if substrate_order_file is None:
            raise ValueError("substrate_order='custom' requires substrate_order_file")
        custom = [
            line.strip()
            for line in substrate_order_file.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        expected_substrates = custom
    elif substrate_order == "alphabetical":
        expected_substrates = sorted(_unique_preserve_order(pools["Name"].astype(str).tolist()))
    else:
        expected_substrates = _unique_preserve_order(pools["Name"].astype(str).tolist())

    matches = _match_features_to_pools(
        filtered=filtered,
        pools=pools,
        tolerance=tolerance,
        rt_tol=rt_tol,
        rt_col=rt_col,
        pools_rt_col=pools_rt_col,
    )
    if not matches:
        empty = pd.DataFrame(index=sample_cols)
        if include_all_substrates and expected_substrates:
            empty = empty.reindex(columns=expected_substrates)
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
                    "rt_feature": m.feature_rt,
                    "rt_pool": m.pool_rt,
                    "delta_rt": m.delta_rt,
                    "pool": m.pool,
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

    if include_all_substrates and expected_substrates:
        matrix = matrix.reindex(columns=expected_substrates)
    else:
        matrix = matrix.sort_index(axis=1)

    return matrix, long_df


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Match filtered mzMine feature table peaks to pool substrates by M+H. "
            "A pool entry matches a feature if pool M+H is within [mz_range.min, mz_range.max] ± tolerance and RT is within ± rt_tol. "
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
        "--rt-tol",
        type=float,
        default=0.1,
        help="Allowed absolute deviation in retention time (same units as filtered 'rt'). Default: 0.1",
    )
    parser.add_argument(
        "--pool",
        type=str,
        default=None,
        help="Optional pool code to filter pools.csv (matches the 'Pool' column).",
    )
    parser.add_argument(
        "--include-all-substrates",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Ensure output matrix has a fixed substrate column set (fills missing with NaN). Default: enabled.",
    )
    parser.add_argument(
        "--substrate-order",
        choices=["pools-file", "alphabetical", "custom"],
        default="pools-file",
        help="Column order for substrates in the output matrix.",
    )
    parser.add_argument(
        "--substrate-order-file",
        type=Path,
        default=None,
        help="If --substrate-order custom: path to a text file with one substrate name per line.",
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

    if args.pool is not None:
        if DEFAULT_POOLS_POOL_COL not in pools.columns:
            raise ValueError(
                f"--pool was provided but pools.csv has no column '{DEFAULT_POOLS_POOL_COL}'."
            )
        pools = pools[pools[DEFAULT_POOLS_POOL_COL] == args.pool].reset_index(drop=True)
        if pools.empty:
            raise ValueError(f"No pools rows matched --pool {args.pool!r}")

    matrix, long_df = build_area_matrix(
        filtered=filtered,
        pools=pools,
        tolerance=args.tolerance,
        rt_tol=args.rt_tol,
        aggregation=args.aggregation,
        include_all_substrates=args.include_all_substrates,
        substrate_order=args.substrate_order,
        substrate_order_file=args.substrate_order_file,
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
