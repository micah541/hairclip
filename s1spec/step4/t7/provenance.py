"""Provenance helpers for the args blocks of T7 outputs: sha256 of input files, of directories of jets, and of code."""
import os, hashlib

STEP4 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()


def sha256_files(paths):
    """one fingerprint for a set of files: sha256 of the concatenated (basename + sha256) over the sorted paths"""
    h = hashlib.sha256()
    for p in sorted(paths): h.update((os.path.basename(p) + sha256(p)).encode())
    return h.hexdigest()


def rel(path):
    return os.path.relpath(os.path.abspath(path), STEP4)


def entry(path):
    """{path (relative to step4), sha256} of a file"""
    return dict(path=rel(path), sha256=sha256(path))


def jets_entry(d, slabs):
    fs = [os.path.join(d, f'jets_{s:03d}.npz') for s in slabs]
    return dict(path=rel(d), slabs=[min(slabs), max(slabs)], sha256_of_files=sha256_files(fs))


def code(paths):
    return {rel(p): sha256(p) for p in paths}
