"""Materialize the cleaned layer: data/processed/*.parquet + issue log.

Run:  .venv/bin/python scripts/build.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ssc_coh.clean import build_processed          # noqa: E402
from ssc_coh.config import PROC_DIR                # noqa: E402
from ssc_coh.features import build_features        # noqa: E402


def main() -> None:
    processed_frames = build_processed()
    processed_frames["features"] = build_features(
        processed_frames["subjects"], processed_frames["ssc_subtype"], processed_frames["vitals"],
        processed_frames["labs"], processed_frames["pft"], processed_frames["mrss"],
        processed_frames["medications"], processed_frames["antibodies"],
        processed_frames["bal"], processed_frames["biopsies"],
        processed_frames["libraries"],
    )

    PROC_DIR.mkdir(parents=True, exist_ok=True)
    for name, frame in processed_frames.items():
        frame.to_parquet(PROC_DIR / f"{name}.parquet", index=False)
    processed_frames["issues"].to_csv(PROC_DIR / "issues.csv", index=False)

    print(f"wrote {len(processed_frames)} frames to {PROC_DIR}\n")
    for name, frame in processed_frames.items():
        print(f"  {name:24s} {frame.shape}")
    print("\n--- issue log ---")
    print(processed_frames["issues"].to_string(index=False, max_colwidth=58))


if __name__ == "__main__":
    main()
