"""Separate Phase 1 CLI. The original pilot CLI and artifacts remain unchanged."""
import argparse
from pathlib import Path
import os
import subprocess
import sys
import time
from .storage import atomic,json_write,read_json,sha,writer_lock

PROJECT=Path(__file__).resolve().parents[1]
if (PROJECT/".deps").exists():
    sys.path.insert(0,str(PROJECT/".deps"))


def validate():
    from .full_pipeline import source_hashes
    directory=PROJECT/"runs/phase1_validation"
    directory.mkdir(parents=True,exist_ok=True)
    stamp=str(time.time_ns())
    (PROJECT/"vt").mkdir(exist_ok=True)
    argv=[sys.executable,"-m","pytest","-p","no:cacheprovider","-v","tests","--basetemp",str(PROJECT/"vt"/stamp[-6:]),
          "--junitxml",str(directory/(stamp+".junit.xml"))]
    env=dict(os.environ); env["PYTHONPATH"]=str(PROJECT/".deps")+os.pathsep+str(PROJECT)
    before=source_hashes()
    result=subprocess.run(argv,cwd=PROJECT,env=env,capture_output=True,text=True)
    output=result.stdout+result.stderr
    atomic(directory/(stamp+".log"),output); atomic(directory/"latest.log",output)
    summary=next((s for s in reversed(output.splitlines()) if " passed" in s or " failed" in s),"See log")
    code=result.returncode if before==source_hashes() else 1
    evidence={"exit_code":code,"summary":summary,"command":argv,"utc_epoch":time.time(),"python":sys.version,"source_sha256":before,
              "original_validation_files_modified":False}
    json_write(directory/"latest.json",evidence)
    json_write(directory/(stamp+".json"),evidence)
    print(output)
    return code


def require_validation():
    from .full_pipeline import source_hashes
    path=PROJECT/"runs/phase1_validation/latest.json"
    if not path.exists() or read_json(path)["exit_code"]!=0 or read_json(path)["source_sha256"]!=source_hashes():
        raise ValueError("Run python -m src.full_cli validate successfully on current source before Phase 1 execution")


def roots(workdir):
    for path in sorted(Path(workdir).glob("full_deletion_overlap_v2_*/experiment_manifest.json")):
        if read_json(path).get("phase")=="phase1":
            yield path.parent


def resume_check(root,config_path,workdir):
    from .full_pipeline import execute
    files=[p for p in root.rglob("*") if p.is_file() and (p.name in ["input.fasta","truth_alignment.fasta","ancestry.json","aligned.fasta","stderr.log","run_done.json","simulation_done.json","block_done.json","realization.json"] or
           any(part in p.parts for part in ["scores","comparisons"]) and p.suffix==".json")]
    def snapshot():
        return {p.relative_to(root).as_posix():{"sha256":sha(p),"mtime_ns":p.stat().st_mtime_ns} for p in files}
    before=snapshot(); summary_sha=sha(root/"full_summary.csv")
    resumed=execute(config_path,workdir)
    assert resumed.resolve()==root.resolve()
    assert snapshot()==before,"Completed artifacts changed on resume"
    assert sha(root/"full_summary.csv")==summary_sha
    invocation=read_json(root/"last_invocation.json")
    assert sum(b["alignments"].get("completed",0) for b in invocation["blocks"])==0
    skipped=sum(b["alignments"].get("skipped",0) for b in invocation["blocks"])
    assert skipped==900
    proof={"passed":True,"alignment_calls_made":0,"alignments_skipped":skipped,"artifacts_checked":len(files),
           "all_bytes_and_mtimes_unchanged":True,"summary_sha256_unchanged":summary_sha,"command":"python -m src.full_cli run --config config/full.json --workdir ./runs"}
    json_write(root/"full_resume_verification.json",proof)
    return proof


def main():
    parser=argparse.ArgumentParser(description="CMPS355 Phase 1 full experiment, separate from the pilot")
    sub=parser.add_subparsers(dest="action",required=True)
    sub.add_parser("validate")
    for name in ["freeze","run","status","report","audit","resume-check","bundle"]:
        cmd=sub.add_parser(name)
        cmd.add_argument("--workdir",default="./runs")
        if name in ["freeze","run","resume-check"]:
            cmd.add_argument("--config",default="config/full.json")
        if name=="run":
            cmd.add_argument("--retry-failed",action="store_true")
    args=parser.parse_args()
    if args.action=="validate":
        return validate()
    if args.action=="freeze":
        from .full_pipeline import prepare
        print(prepare(args.config,args.workdir)[0]); return 0
    if args.action=="run":
        require_validation()
        from .full_pipeline import execute
        from .full_reporting import inventory
        root=execute(args.config,args.workdir,args.retry_failed)
        jobs,rows=inventory(root)
        print(root)
        return 0 if len(rows)==900 and all(j["score_complete"] and j["comparison_complete"] for j in jobs) else 1
    selected=list(roots(args.workdir))
    if not selected:
        parser.error("No Phase 1 experiment found; freeze or run first")
    for root in selected:
        if args.action=="status":
            # Read-only status must remain available while the experiment owns the writer lock.
            from .full_reporting import inventory
            jobs,rows=inventory(root)
            print(f"{root.name}: {len(jobs)} planned; {sum(j['alignment_state']=='completed' for j in jobs)} completed; "
                  f"{sum(j['alignment_state']=='failed' for j in jobs)} failed; {len(rows)} validated score rows")
            continue
        if args.action=="resume-check":
            require_validation()
            print(resume_check(root,args.config,args.workdir)); continue
        with writer_lock(root):
            if args.action=="audit":
                from .full_audit import audit
                print(audit(root))
            elif args.action in ["report","bundle"]:
                from .full_reporting import build
                build(root)
                if args.action=="bundle":
                    from .full_audit import preserve_phase0
                    preserve_phase0(root)
                    from .full_bundle import bundle
                    print(bundle(root))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
