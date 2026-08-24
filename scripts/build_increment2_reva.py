#!/usr/bin/env python3
"""Regenerate Increment 2 Revision A checksums and rebuild the Helena archive."""

import hashlib
import os
import tarfile

WORKSPACE_ROOT = "/home/himanshu/bdb-qchipms"
PKG_DIR = os.path.join(WORKSPACE_ROOT, "DOCUMENTS", "Helena", "tcf7l2_regulatory_increment2")
TAR_OUT = os.path.join(WORKSPACE_ROOT, "DOCUMENTS", "Helena", "tcf7l2_regulatory_increment2_revA.tar.gz")
CHECKSUM_NAME = "checksums.sha256"


def compute_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(8192):
            h.update(chunk)
    return h.hexdigest()


def main():
    lines = []
    for root, dirs, files in os.walk(PKG_DIR):
        dirs.sort()
        for name in sorted(files):
            if name == CHECKSUM_NAME:
                continue
            full = os.path.join(root, name)
            rel = os.path.relpath(full, PKG_DIR)
            lines.append(f"{compute_sha256(full)}  {rel}")

    with open(os.path.join(PKG_DIR, CHECKSUM_NAME), "w") as f:
        f.write("\n".join(lines) + "\n")

    for line in lines:
        print(line)

    if os.path.exists(TAR_OUT):
        os.remove(TAR_OUT)
    with tarfile.open(TAR_OUT, "w:gz") as tar:
        tar.add(PKG_DIR, arcname="tcf7l2_regulatory_increment2")

    sha = compute_sha256(TAR_OUT)
    size_kb = os.path.getsize(TAR_OUT) / 1024.0
    print()
    print(f"archive: {TAR_OUT}")
    print(f"size_kb: {size_kb:.2f}")
    print(f"sha256:  {sha}")


if __name__ == "__main__":
    main()
