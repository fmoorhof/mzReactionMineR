## Plan: RT-Aware Pool Matching + Standardized Matrix

Update the pool-matching workflow so the output matrix is standardized (same substrate columns/order every run, even when missing) and so matches require both **m/z** (±0.0125) and **retention time** (RT within a tolerance around pools.csv “Retention time”). This prevents assigning one peak to multiple substrates when they share the same $m/z$ but differ in RT, and it stabilizes downstream plate-to-plate comparisons.

**Steps**
1. Extend pools parsing in `python/match_pools.py` to keep and parse additional columns:
   - Parse `Pool` and `Retention time` (decimal commas → floats); keep `Name`, `M+H`, `Retention time` in memory.
   - Decide behavior for missing pool RT: with your chosen “hard filter”, candidates with missing RT should not match (otherwise they bypass RT and can reintroduce duplicates).
2. Add RT constraints to matching:
   - Read feature RT from `filtered_for_mzmine.csv` column `rt` (already present).
   - A candidate (substrate) matches a feature only if:
     - `M+H` is within `[mz_range.min, mz_range.max] ± 0.0125`, AND
     - `abs(rt_feature - rt_pool) <= rt_tol` (new CLI flag `--rt-tol`).
3. Resolve remaining ambiguous matches deterministically:
   - If multiple substrates still pass both filters for the same feature, select the substrate with the smallest `abs(rt_feature - rt_pool)` (tie-breaker), and drop the others for that feature.
   - Keep traceability fields in the long output (e.g., `rt_feature`, `rt_pool`, `delta_rt`) so you can audit assignments.
4. Enforce a stable substrate column order and presence:
   - Add `--include-all-substrates` (default on for your use case) to ensure every substrate from pools.csv appears as a column even if nothing matched (filled as NaN/blank).
   - Add `--substrate-order pools-file|alphabetical|custom`:
     - Default for your standardization: `pools-file` order (as listed in pools.csv, optionally grouped by `Pool` if desired).
     - `custom` optionally takes `--substrate-order-file` (one substrate name per line) if you need a very specific plate standard.
5. Keep outputs compatible with your current workflow:
   - Wide matrix CSV: sample rows × substrate columns (NaN where missing), same as today but standardized.
   - Long CSV: one row per accepted match including `sample`, `substrate`, `feature_id`, `area`, `pool_mh`, `mz_range_min/max`, plus RT fields.
6. Update CLI help text and error messages:
   - If RT filtering is enabled but `--rt-tol` is not provided (or is non-positive), exit with a clear message.
   - If `Retention time` cannot be parsed for many pools, print a short warning summary.

**Verification**
- Run (example): `python python/match_pools.py --filtered filtered_for_mzmine.csv --pools pools.csv --tolerance 0.0125 --rt-tol 0.1 --include-all-substrates --substrate-order pools-file`
- Check:
  - The wide output has the same columns every run (even if no peaks match some substrates).
  - The long output shows no feature duplicated across multiple substrates (unless RT truly overlaps and tie-breaker selects one).
  - A known ambiguous $m/z$ (e.g., shared “365.1166”) stops producing duplicate substrate assignments unless RT supports it.

**Decisions**
- Uses a single global matrix with stable ordering.
- Missing entries remain NaN/blank.
- RT is a hard filter using pools.csv “Retention time” ± `--rt-tol`, plus “closest RT wins” if multiple candidates still pass.

If you tell me what RT tolerance you typically consider acceptable (e.g., `0.05`, `0.1`, `0.2` minutes), I’ll bake that in as the default for `--rt-tol` in the implementation.
