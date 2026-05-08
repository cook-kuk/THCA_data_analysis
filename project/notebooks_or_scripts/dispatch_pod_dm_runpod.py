#!/usr/bin/env python3
"""
RunPod Pod DM auto-dispatcher.

Steps:
  1. Create L40S 48GB pod with PyTorch + CUDA image
  2. Wait until SSH-ready
  3. Upload manifest + dispatch script via runpod cp / scp
  4. Execute bash dispatch in tmux on pod
  5. Print pod ID + SSH command for user monitoring

Usage:
  RUNPOD_API_KEY=rpa_... python3 dispatch_pod_dm_runpod.py
"""
from __future__ import annotations
import os, sys, time, subprocess
from pathlib import Path

import runpod

ROOT = Path(__file__).resolve().parent.parent.parent

# Read API key
KEY = os.environ.get("RUNPOD_API_KEY")
if not KEY:
    cfg = Path("/home/seungho/.runpod/config")
    if cfg.exists():
        for line in cfg.read_text().splitlines():
            if line.startswith("RUNPOD_API_KEY="):
                KEY = line.split("=", 1)[1].strip()
                break
if not KEY:
    print("RUNPOD_API_KEY missing"); sys.exit(1)
runpod.api_key = KEY

POD_NAME = "thca-spark-dm-l40s"
GPU_TYPE = "NVIDIA L40S"
DISK_GB = 100
VOLUME_GB = 50
IMAGE = "runpod/pytorch:2.4.0-py3.11-cuda12.4.1-devel-ubuntu22.04"
PORTS = "22/tcp,8888/http"

# Bash that runs on pod after SSH connect
DISPATCH_SCRIPT = ROOT / "project/notebooks_or_scripts/runpod_pod_DM_WSI_pathology.sh"
MANIFEST_TSV = ROOT / "project/data/manifests/tcga_thca_wsi_dm_balanced.tsv"
MANIFEST_JSON = ROOT / "project/data/manifests/tcga_thca_wsi_gdc_filter.json"


def main():
    print(f"=== creating pod {POD_NAME} on {GPU_TYPE} ===")
    pod = runpod.create_pod(
        name=POD_NAME,
        image_name=IMAGE,
        gpu_type_id=GPU_TYPE,
        cloud_type="ALL",
        gpu_count=1,
        volume_in_gb=VOLUME_GB,
        container_disk_in_gb=DISK_GB,
        ports=PORTS,
        env={
            "JUPYTER_PASSWORD": "thca-spark-dm",
        },
        support_public_ip=True,
        start_ssh=True,
    )
    pid = pod.get("id")
    print(f"  pod id: {pid}")
    print(f"  status: {pod.get('desiredStatus')}")
    Path("/tmp/pod_id.txt").write_text(pid)

    # wait for running
    print("\n=== waiting for pod RUNNING ===")
    for _ in range(60):
        time.sleep(10)
        info = runpod.get_pod(pid)
        st = info.get("desiredStatus")
        running = info.get("runtime", {}).get("ports", [])
        print(f"  status={st}, runtime ports={[p.get('publicPort') for p in running]}")
        if st == "RUNNING" and running:
            break
    else:
        print("timeout waiting for RUNNING")
        sys.exit(1)

    # find SSH host
    ssh_port = None
    ssh_ip = None
    for p in info.get("runtime", {}).get("ports", []):
        if p.get("privatePort") == 22:
            ssh_port = p.get("publicPort")
            ssh_ip = p.get("ip")
    print(f"\n=== SSH endpoint: {ssh_ip}:{ssh_port} ===")
    print(f"manual SSH: ssh root@{ssh_ip} -p {ssh_port}")
    print(f"\nSCP upload examples (using your local SSH key):")
    print(f"  scp -P {ssh_port} {DISPATCH_SCRIPT} root@{ssh_ip}:/workspace/")
    print(f"  scp -P {ssh_port} {MANIFEST_TSV} root@{ssh_ip}:/workspace/wsi_pathology_dm/manifests/")
    print(f"  scp -P {ssh_port} {MANIFEST_JSON} root@{ssh_ip}:/workspace/wsi_pathology_dm/manifests/")
    print(f"\nRun on pod (tmux for resilience):")
    print(f"  ssh root@{ssh_ip} -p {ssh_port} 'bash -lc \"")
    print(f"    apt-get install -y -qq tmux ; ")
    print(f"    mkdir -p /workspace/wsi_pathology_dm/manifests ; ")
    print(f"    tmux new -d -s podDM 'bash /workspace/runpod_pod_DM_WSI_pathology.sh > /workspace/podDM.log 2>&1' ; ")
    print(f"    sleep 2; tmux ls\"'")
    print(f"\nMonitor: ssh root@{ssh_ip} -p {ssh_port} 'tail -f /workspace/podDM.log'")
    print(f"\nStop pod: python3 -c 'import runpod; runpod.api_key=\"...\"; runpod.stop_pod(\"{pid}\")'")
    print(f"\nIMPORTANT: pod is RUNNING and accruing cost (~$0.50/h). Stop when done.")


if __name__ == "__main__":
    main()
