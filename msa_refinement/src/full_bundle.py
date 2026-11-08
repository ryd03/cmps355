"""Package complete Phase 1 evidence without touching the Phase 0 archive."""
import json
import os
from pathlib import Path
import uuid
import zipfile
from .full_pipeline import PROJECT
from .storage import digest, json_write, read_json, sha


def bundle(root):
    root=Path(root).resolve()
    for name in ["full_artifact_audit.json","full_resume_verification.json","phase0_preservation.json"]:
        if not read_json(root/name)["passed"]:
            raise ValueError("Required evidence failed: "+name)
    target=PROJECT/"full_review_bundle.zip"
    temp=target.with_name(uuid.uuid4().hex[:12]+".partial")
    files={}
    for folder in ["src","tests","config","data"]:
        for path in (PROJECT/folder).rglob("*"):
            if path.is_file() and "__pycache__" not in path.parts:
                files[path.relative_to(PROJECT).as_posix()]=path
    for path in PROJECT.glob("*.md"):
        files[("protocol/README_PHASE0.md" if path.name=="README.md" else path.name)]=path
    for path in PROJECT.glob("*.txt"):
        files[path.name]=path
    for name in ["setup_mafft.py","installation.json"]:
        files["tools/"+name]=PROJECT/"tools"/name
    for path in root.rglob("*"):
        if path.is_file() and path.name not in [".writer.lock","full_bundle_verification.json"]:
            files[path.relative_to(PROJECT).as_posix()]=path
    for path in (PROJECT/"runs/phase1_validation").glob("*"):
        if path.is_file():
            files[path.relative_to(PROJECT).as_posix()]=path
    provenance=PROJECT/"runs/phase1_provenance/phase0_before.json"
    files[provenance.relative_to(PROJECT).as_posix()]=provenance
    files["README.md"]=PROJECT/"PHASE1_REPRODUCE.md"
    for name in ["full_report.md","full_summary.csv","paired_comparisons.csv","confidence_intervals.csv","cell_summary.csv","full_job_manifest.csv"]:
        files[name]=root/name
    hashes={}
    with zipfile.ZipFile(temp,"w",zipfile.ZIP_DEFLATED) as archive:
        for name,path in sorted(files.items()):
            archive.write(path,name)
            hashes[name]=sha(path)
        archive.writestr("bundle_sha256_manifest.json",json.dumps(hashes,sort_keys=True,indent=2))
    with zipfile.ZipFile(temp) as archive:
        assert archive.testzip() is None
        assert all(digest(archive.read(name))==expected for name,expected in hashes.items())
        assert len([n for n in hashes if n.endswith('/aligned.fasta') and '/quarantine/' not in n])==900
    os.replace(temp,target)
    json_write(root/"full_bundle_verification.json",{"passed":True,"members_verified":len(hashes),"zip_sha256":sha(target),"zip_crc_and_member_sha256_verified":True})
    return target
