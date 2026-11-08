"""Subprocess MAFFT adapter; accepts a JSON argv prefix for portable wrappers."""
import json
import os
import re
from pathlib import Path
import shutil
import subprocess

RUNNER_VERSION = "mafft-fixed-flags-v2.1"


def executable_prefix():
    if os.environ.get("PILOT_MAFFT_COMMAND"):
        return json.loads(os.environ["PILOT_MAFFT_COMMAND"])
    installed = shutil.which("mafft")
    if installed:
        return [installed]
    bundled = Path(__file__).resolve().parents[1] / "tools" / "mafft-win" / "mafft.bat"
    if bundled.exists():
        return [str(bundled)]
    raise RuntimeError("MAFFT not found. Install MAFFT or set PILOT_MAFFT_COMMAND to a JSON argv prefix.")


def command(prefix, source, cap, config):
    m = config["mafft"]
    return [*prefix, "--amino", "--thread", str(m["thread"]), "--retree", str(m["retree"]),
            "--maxiterate", str(cap), str(Path(source).resolve())]


def version(prefix):
    result = subprocess.run([*prefix, "--version"], capture_output=True, timeout=30)
    if result.returncode:
        raise RuntimeError("MAFFT version check failed: " + result.stderr.decode(errors="replace"))
    output = (result.stdout + result.stderr).decode(errors="replace").strip()
    match = re.search(r"v\d+\.\d+[^\r\n]*", output)
    return match.group(0) if match else output
