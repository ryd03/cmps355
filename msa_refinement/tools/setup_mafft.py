"""Install official portable Windows MAFFT locally, preserving ZIP provenance."""
import hashlib
import json
from pathlib import Path
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parent
URL = "https://mafft.cbrc.jp/alignment/software/mafft-7.526-win64-signed.zip"


def main():
    archive = ROOT / "mafft-7.526-win64-signed.zip"
    if not archive.exists():
        urllib.request.urlretrieve(URL, archive)
    sha = hashlib.sha256(archive.read_bytes()).hexdigest()
    provenance = ROOT / "installation.json"
    if provenance.exists():
        expected = json.loads(provenance.read_text())["download_sha256"]
        if expected != sha:
            raise ValueError("Downloaded package differs from recorded SHA256")
    with zipfile.ZipFile(archive) as package:
        assert package.testzip() is None
        package.extractall(ROOT)
    provenance.write_text(json.dumps({"download_url":URL,"download_sha256":sha,
        "official_instructions":"https://mafft.cbrc.jp/alignment/software/windows_without_cygwin.html",
        "version":"v7.526 (2024/Apr/26)","platform":"official Windows all-in-one (MSYS)",
        "install_command":"python tools/setup_mafft.py"},indent=2)+"\n")
    print("Installed",sha)


if __name__ == "__main__":
    main()
