#!/usr/bin/env python3
"""Install the Helical SDK on Windows by removing the problematic louvain pin.

Helical v1.4.5 depends on ``louvain==0.8.2`` which only ships Linux wheels and
fails to build from source on Windows.  This helper downloads the official
wheel, strips the louvain requirement from its metadata, repacks it, and then
lets ``pip`` install the patched wheel (which still installs every other
dependency declared by Helical).

Usage:
    python scripts/install_helical.py [--version 1.4.5]
"""

from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path


DEFAULT_VERSION = "1.4.5"


def run(cmd: list[str]) -> None:
    """Run a subprocess, streaming stdout/stderr, and fail loudly on error."""
    completed = subprocess.run(cmd, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"Command failed ({completed.returncode}): {' '.join(cmd)}")


def patch_wheel_metadata(extracted_dir: Path, version: str) -> None:
    """Remove the louvain requirement and update the RECORD hash."""
    dist_info = extracted_dir / f"helical-{version}.dist-info"
    metadata_path = dist_info / "METADATA"
    record_path = dist_info / "RECORD"

    original_lines = metadata_path.read_text().splitlines()
    filtered_lines = [
        line for line in original_lines if not line.startswith("Requires-Dist: louvain")
    ]
    metadata_path.write_text("\n".join(filtered_lines) + "\n")

    metadata_bytes = metadata_path.read_bytes()
    metadata_hash = base64.urlsafe_b64encode(
        hashlib.sha256(metadata_bytes).digest()
    ).rstrip(b"=").decode()

    # Update the RECORD entry so pip's integrity check passes.
    updated_rows: list[list[str]] = []
    with record_path.open(newline="") as handle:
        reader = csv.reader(handle)
        for row in reader:
            if row and row[0] == f"helical-{version}.dist-info/METADATA":
                updated_rows.append(
                    [row[0], f"sha256={metadata_hash}", str(len(metadata_bytes))]
                )
            else:
                updated_rows.append(row)

    with record_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerows(updated_rows)


def extract_wheel(wheel_path: Path, destination: Path) -> None:
    """Unpack the downloaded wheel into the temporary directory."""
    with zipfile.ZipFile(wheel_path) as zf:
        zf.extractall(destination)


def build_patched_wheel(src_dir: Path, dest_wheel: Path) -> None:
    """Zip the extracted tree back into a wheel."""
    with zipfile.ZipFile(dest_wheel, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(src_dir.rglob("*")):
            if path.is_dir():
                continue
            arcname = path.relative_to(src_dir).as_posix()
            zf.write(path, arcname)


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", default=DEFAULT_VERSION, help="Helical version")
    parser.add_argument(
        "--pip",
        default=sys.executable,
        help="Python executable whose pip should be used (defaults to current interpreter).",
    )
    args = parser.parse_args(argv)

    python_exe = Path(args.pip)
    pip_cmd = [str(python_exe), "-m", "pip"]

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        download_dir = tmp_path / "download"
        download_dir.mkdir(parents=True, exist_ok=True)

        print(f"Downloading Helical {args.version} wheel...")
        run(
            pip_cmd
            + [
                "download",
                f"helical=={args.version}",
                "--no-deps",
                "-d",
                str(download_dir),
            ]
        )

        wheels = list(download_dir.glob("helical-*.whl"))
        if not wheels:
            raise FileNotFoundError("Unable to find downloaded helical wheel")
        wheel_path = wheels[0]

        extracted_dir = tmp_path / "extracted"
        extract_wheel(wheel_path, extracted_dir)

        patch_wheel_metadata(extracted_dir, args.version)

        patched_wheel = tmp_path / wheel_path.name
        build_patched_wheel(extracted_dir, patched_wheel)

        print("Installing patched Helical wheel (louvain dependency removed)...")
        run(pip_cmd + ["install", str(patched_wheel)])

    print("Helical installation complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
