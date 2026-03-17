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
DEFAULT_XLSX_OUT = "matched_pool_areas.xlsx"


def _normalize_sample_key(value: str) -> str:
    s = str(value).strip()
    if s.startswith("datafile."):
        s = s[len("datafile.") :]
    # R sometimes prefixes column names that start with digits with 'X'
    if len(s) > 1 and s[0] == "X" and s[1].isdigit():
        s = s[1:]
    # keep only the filename part if paths slipped in
    s = Path(s).name
    return s


def _read_plate_subpools(
    path: Path,
    filename_col: str = "filename",
    subpool_col: str = "Subpool",
) -> dict[str, str]:
    df = pd.read_csv(path)
    missing = [c for c in [filename_col, subpool_col] if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in {path}: {missing}")

    out: dict[str, str] = {}
    for _, row in df.iterrows():
        fn = _normalize_sample_key(row[filename_col])
        sp = str(row[subpool_col]).strip()
        if fn and sp and sp.lower() != "nan":
            out[fn] = sp
            out[Path(fn).stem] = sp  # also allow matching without extension
    return out


def _parse_float_series(series: pd.Series) -> pd.Series:
    """Parse numeric columns that may contain 'NA', empty, or decimal commas."""
    s = series.astype(str).str.strip()
    s = s.replace(["", "NA", "NaN", "nan"], pd.NA)
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
    if not isinstance(stop_idx, int):
        raise ValueError(
            f"Column '{META_STOP_COL}' appears multiple times in {path}; cannot infer sample columns reliably."
        )
    sample_cols = list(df.columns)[stop_idx + 1 :]
    if not sample_cols:
        raise ValueError(f"No sample columns found after '{META_STOP_COL}' in {path}")

    for col in sample_cols:
        df[col] = _parse_float_series(df[col])

    # Optional columns (only used for reporting in extra sheets)
    if "mz" in df.columns:
        df["mz"] = _parse_float_series(df["mz"])

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


def _safe_sheet_name(name: str) -> str:
    # Excel sheet names: max 31 chars, cannot contain: : \ / ? * [ ]
    bad = [":", "\\", "/", "?", "*", "[", "]"]
    for ch in bad:
        name = name.replace(ch, "-")
    name = name.strip() or "Sheet"
    return name[:31]


def _match_features_to_pools(
    filtered: pd.DataFrame,
    pools: pd.DataFrame,
    tolerance: float,
    rt_tol: float,
    rt_col: str = DEFAULT_RT_COL,
    pools_rt_col: str = DEFAULT_POOLS_RT_COL,
) -> list[Match]:
    # IMPORTANT: We do *not* require an RT hit. RT is used only to rank matches.
    matches: list[Match] = []

    mh_values = pools["M+H"].to_numpy()
    names = pools["Name"].to_list()
    rt_values = (
        pools[pools_rt_col].to_numpy()
        if pools_rt_col in pools.columns
        else [float("nan")] * len(pools)
    )
    pools_values = (
        pools[DEFAULT_POOLS_POOL_COL].to_list()
        if DEFAULT_POOLS_POOL_COL in pools.columns
        else [None] * len(pools)
    )

    mz_min = filtered["mz_range.min"].to_numpy()
    mz_max = filtered["mz_range.max"].to_numpy()
    feature_rt = (
        filtered[rt_col].to_numpy() if rt_col in filtered.columns else [float("nan")] * len(filtered)
    )

    for i in range(len(filtered)):
        lo = mz_min[i]
        hi = mz_max[i]
        rt_f = feature_rt[i] if i < len(feature_rt) else float("nan")
        if pd.isna(lo) or pd.isna(hi):
            continue

        lo2 = lo - tolerance
        hi2 = hi + tolerance

        for substrate, mh, rt_p, pool_code in zip(names, mh_values, rt_values, pools_values):
            if not (lo2 <= mh <= hi2):
                continue

            if not pd.isna(rt_f) and not pd.isna(rt_p):
                delta = abs(float(rt_f) - float(rt_p))
            else:
                delta = float("nan")

            matches.append(
                Match(
                    feature_idx=i,
                    substrate=substrate,
                    pool_mh=float(mh),
                    pool_rt=float(rt_p) if not pd.isna(rt_p) else float("nan"),
                    feature_rt=float(rt_f) if not pd.isna(rt_f) else float("nan"),
                    delta_rt=float(delta) if not pd.isna(delta) else float("nan"),
                    pool=pool_code,
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
    sample_to_subpool: dict[str, str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    stop_idx = filtered.columns.get_loc(META_STOP_COL)
    if not isinstance(stop_idx, int):
        raise ValueError(
            f"Column '{META_STOP_COL}' appears multiple times in filtered table; cannot infer sample columns reliably."
        )
    sample_cols = list(filtered.columns)[stop_idx + 1 :]

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
        mz_feature = row.get("mz") if "mz" in filtered.columns else float("nan")

        for sample in sample_cols:
            sample_key = _normalize_sample_key(sample)
            sample_subpool = None
            if sample_to_subpool is not None:
                sample_subpool = sample_to_subpool.get(sample_key) or sample_to_subpool.get(Path(sample_key).stem)

            rt_hit = False
            if rt_tol and rt_tol > 0 and pd.notna(m.delta_rt):
                rt_hit = bool(float(m.delta_rt) <= float(rt_tol))

            pool_hit = False
            if sample_subpool is not None and m.pool is not None and str(m.pool).strip() != "":
                pool_hit = str(m.pool).strip() == str(sample_subpool).strip()

            if sample_subpool is None:
                tier = 1 if rt_hit else 2
            else:
                tier = 1 if (pool_hit and rt_hit) else 2 if pool_hit else 3 if rt_hit else 4

            if pool_hit and rt_hit:
                match_class = "mz+rt+pool"
            elif pool_hit:
                match_class = "mz+pool"
            elif rt_hit:
                match_class = "mz+rt"
            else:
                match_class = "mz"

            records.append(
                {
                    "sample": sample,
                    "substrate": m.substrate,
                    "area": row[sample],
                    "feature_id": feature_id,
                    "pool_mh": m.pool_mh,
                    "mz_range_min": lo,
                    "mz_range_max": hi,
                    "mz_feature": mz_feature,
                    "rt_feature": m.feature_rt,
                    "rt_pool": m.pool_rt,
                    "delta_rt": m.delta_rt,
                    "pool": m.pool,
                    "rt_hit": rt_hit,
                    "sample_subpool": sample_subpool,
                    "pool_hit": pool_hit if sample_subpool is not None else pd.NA,
                    "tier": tier,
                    "match_class": match_class,
                }
            )

    long_df = pd.DataFrame.from_records(records)

    # Aggregate within the *best* tier per (sample, substrate)
    best_tier = (
        long_df.groupby(["sample", "substrate"], as_index=False)
        .agg(best_tier=("tier", "min"))
    )
    tmp = long_df.merge(best_tier, on=["sample", "substrate"], how="left")
    tmp = tmp[tmp["tier"] == tmp["best_tier"]]

    grouped = (
        tmp.groupby(["sample", "substrate"], as_index=False)["area"]
        .apply(lambda s: _aggregate(s, aggregation))
    )

    matrix = grouped.pivot(index="sample", columns="substrate", values="area")
    # Preserve input order (samples as they appear in the feature table; substrates as configured)
    matrix = matrix.reindex(index=sample_cols)
    if include_all_substrates and expected_substrates:
        matrix = matrix.reindex(columns=expected_substrates)
    else:
        matrix = matrix.reindex(columns=sorted(matrix.columns.tolist()))

    return matrix, long_df


def build_info_sheets(
    long_df: pd.DataFrame,
    sample_order: list[str],
    substrate_order: list[str],
) -> dict[str, pd.DataFrame]:
    """Create per-field matrices (same shape as the area matrix).

    For metadata (rt/mz/etc), we select the feature with the highest area for each
    (sample, substrate) and report its metadata.
    """
    sheets: dict[str, pd.DataFrame] = {}
    if long_df.empty:
        return sheets

    # Pick a representative feature per (sample, substrate):
    #   - restrict to the best tier first (pool/rt hierarchy)
    #   - then pick the max-area feature
    tmp = long_df.dropna(subset=["area"]).copy()
    if tmp.empty:
        return sheets

    if "tier" in tmp.columns:
        min_tier = tmp.groupby(["sample", "substrate"])["tier"].transform("min")
        tmp = tmp[tmp["tier"] == min_tier]

    idx = tmp.groupby(["sample", "substrate"])["area"].idxmax()
    best = tmp.loc[idx]

    def pivot_field(field: str) -> pd.DataFrame:
        df = best.pivot(index="sample", columns="substrate", values=field)
        df = df.reindex(index=sample_order, columns=substrate_order)
        return df

    # One sheet per info field
    sheets["feature_id"] = pivot_field("feature_id")
    sheets["mz_feature"] = pivot_field("mz_feature")
    sheets["mz_range_min"] = pivot_field("mz_range_min")
    sheets["mz_range_max"] = pivot_field("mz_range_max")
    sheets["rt_feature"] = pivot_field("rt_feature")
    sheets["rt_pool"] = pivot_field("rt_pool")
    sheets["delta_rt"] = pivot_field("delta_rt")
    sheets["pool_mh"] = pivot_field("pool_mh")
    sheets["pool"] = pivot_field("pool")

    # Matching hierarchy/debugging
    if "tier" in best.columns:
        sheets["tier"] = pivot_field("tier")
    if "match_class" in best.columns:
        sheets["match_class"] = pivot_field("match_class")
    if "rt_hit" in best.columns:
        sheets["rt_hit"] = pivot_field("rt_hit")
    if "pool_hit" in best.columns:
        sheets["pool_hit"] = pivot_field("pool_hit")
    if "sample_subpool" in best.columns:
        sheets["sample_subpool"] = pivot_field("sample_subpool")

    return sheets


def write_excel_multi_sheet(
    out_xlsx: Path,
    area_matrix: pd.DataFrame,
    info_sheets: dict[str, pd.DataFrame],
    long_df: pd.DataFrame | None = None,
) -> None:
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        # First sheet: the current "short" area matrix
        area_matrix.to_excel(writer, sheet_name="area")

        # Optional: full long table for traceability in Excel
        if long_df is not None and not long_df.empty:
            long_df.to_excel(writer, sheet_name="matches_long", index=False)

        for name, df in info_sheets.items():
            df.to_excel(writer, sheet_name=_safe_sheet_name(name))


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Match filtered mzMine feature table peaks to pool substrates by M+H. "
            "A pool entry matches a feature if pool M+H is within [mz_range.min, mz_range.max] ± tolerance. "
            "RT and plate subpool are used to rank/prioritize matches (not to drop m/z hits). "
            "Outputs a sample×substrate area matrix plus traceability sheets."
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
        default=Path("mock_data/pools.csv"),
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
        default=0.05,
        help="Allowed absolute deviation in retention time (same units as filtered 'rt', most likely minutes). Default: 0.05",
    )
    parser.add_argument(
        "--plate",
        type=Path,
        default=Path("mock_data/plate_rep1.csv"),
        help="Optional plate metadata CSV (e.g. mock_data/plate5.csv) to prioritize matches by sample subpool.",
    )
    parser.add_argument(
        "--plate-filename-col",
        type=str,
        default="filename",
        help="Column in --plate that contains the filename (must match sample column names). Default: filename",
    )
    parser.add_argument(
        "--plate-subpool-col",
        type=str,
        default="Subpool",
        help="Column in --plate that contains the subpool code. Default: Subpool",
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
        "--out-xlsx",
        type=Path,
        default=Path(DEFAULT_XLSX_OUT),
        help="Output XLSX path containing multiple sheets (area first, then RT/mz/etc).",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="(Optional) Output CSV path for the sample×substrate area matrix.",
    )
    parser.add_argument(
        "--out-long",
        type=Path,
        default=None,
        help="(Optional) Output CSV path for the long match table (traceability).",
    )

    args = parser.parse_args()

    filtered = _read_filtered(args.filtered)
    pools = _read_pools(args.pools)

    sample_to_subpool = None
    if args.plate is not None:
        sample_to_subpool = _read_plate_subpools(
            args.plate,
            filename_col=args.plate_filename_col,
            subpool_col=args.plate_subpool_col,
        )

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
        sample_to_subpool=sample_to_subpool,
    )

    # Preserve explicit sample/substrate orders in all sheets
    stop_idx = filtered.columns.get_loc(META_STOP_COL)
    if not isinstance(stop_idx, int):
        raise ValueError(
            f"Column '{META_STOP_COL}' appears multiple times in filtered table; cannot infer sample columns reliably."
        )
    sample_cols = list(filtered.columns)[stop_idx + 1 :]
    if args.substrate_order == "custom" and args.substrate_order_file is not None:
        substrate_cols = [
            line.strip()
            for line in args.substrate_order_file.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
    elif args.substrate_order == "alphabetical":
        substrate_cols = sorted(_unique_preserve_order(pools["Name"].astype(str).tolist()))
    else:
        substrate_cols = _unique_preserve_order(pools["Name"].astype(str).tolist())

    info_sheets = build_info_sheets(
        long_df=long_df,
        sample_order=sample_cols,
        substrate_order=substrate_cols,
    )
    write_excel_multi_sheet(args.out_xlsx, matrix, info_sheets, long_df=long_df)
    print(f"Wrote Excel workbook: {args.out_xlsx} (area shape={matrix.shape})")

    if args.out is not None:
        matrix.to_csv(args.out, index=True)
        print(f"Wrote area CSV: {args.out}")
    if args.out_long is not None and not long_df.empty:
        long_df.to_csv(args.out_long, index=False)
        print(f"Wrote long CSV: {args.out_long} (rows={len(long_df)})")


if __name__ == "__main__":
    main()
