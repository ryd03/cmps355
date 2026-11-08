"""Atomic artifacts, verified stage markers, and an OS-released writer lock."""
import contextlib
import hashlib
import json
import os
from pathlib import Path
import time
import uuid


def digest(data):
    return hashlib.sha256(data).hexdigest()


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def sha(path):
    return digest(Path(path).read_bytes())


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def atomic(path, content):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(uuid.uuid4().hex[:16] + ".partial")
    data = content.encode() if isinstance(content, str) else content
    with tmp.open("wb") as out:
        out.write(data)
        out.flush()
        os.fsync(out.fileno())
    os.replace(tmp, path)


def json_write(path, obj):
    atomic(path, json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + "\n")


def marker(path, files, dependencies, details=None):
    path = Path(path)
    json_write(path, {"files": {str(Path(f).relative_to(path.parent)): sha(f) for f in files},
                      "dependencies": dependencies, "details": details or {}, "utc_epoch": time.time()})


def complete(path, dependencies, validator=None):
    path = Path(path)
    try:
        obj = read_json(path)
        if obj["dependencies"] != dependencies or not obj["files"]:
            return False
        if any(sha(path.parent / f) != h for f, h in obj["files"].items()):
            return False
        if validator:
            validator()
        return True
    except (OSError, ValueError, KeyError, TypeError, AssertionError):
        return False


def quarantine(directory, paths):
    """Preserve invalid published files and abandoned partials before replacement."""
    directory = Path(directory)
    candidates = set(Path(p) for p in paths if Path(p).exists())
    candidates.update(directory.glob("*.partial"))
    if candidates:
        dest = directory / "quarantine" / uuid.uuid4().hex[:12]
        dest.mkdir(parents=True)
        for path in candidates:
            os.replace(path, dest / path.name)


@contextlib.contextmanager
def writer_lock(root):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    with (root / ".writer.lock").open("a+b") as handle:
        handle.seek(0, os.SEEK_END)
        if handle.tell() == 0:
            handle.write(b"0")
            handle.flush()
        handle.seek(0)
        try:
            if os.name == "nt":
                import msvcrt
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            raise RuntimeError("Another process holds the exclusive experiment writer lock") from exc
        try:
            yield
        finally:
            handle.seek(0)
            if os.name == "nt":
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(handle, fcntl.LOCK_UN)
