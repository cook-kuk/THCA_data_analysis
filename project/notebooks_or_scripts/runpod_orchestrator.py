#!/usr/bin/env python3
"""RunPod multi-pod orchestrator — D + B + C parallel.

Marathon mode: D = paper-blocking (D7-P3 closure); B + C = Phase 0 stretching.

Auto-shutdown: 24hr forced.
Spend limit: $5/hr advised (currently $80/hr).
Cost ceiling: $40 (D ~$5-7 + C ~$12-15 + B PoC ~$15-20).
"""
import os, sys, json, time
import requests

API_KEY = open(os.path.expanduser("~/.runpod/config")).read().split("=")[1].strip()
GQL = "https://api.runpod.io/graphql"
REST = "https://rest.runpod.io/v1"

HEADERS = {"Authorization": f"Bearer {API_KEY}", "Content-Type": "application/json"}

def gql(query, variables=None):
    r = requests.post(GQL, json={"query": query, "variables": variables or {}}, headers=HEADERS, timeout=30)
    r.raise_for_status()
    return r.json()

def rest(method, path, **kw):
    r = requests.request(method, f"{REST}{path}",
                          headers={"Authorization": f"Bearer {API_KEY}"},
                          timeout=60, **kw)
    return r

def balance():
    q = "query { myself { id email spendLimit clientBalance } }"
    return gql(q)["data"]["myself"]

def list_pods():
    return rest("GET", "/pods").json()

def pod_create(name, gpu_type_id, gpu_count, container_image, disk_gb, ports="22/tcp",
                env=None, volume_mount_path=None, network_volume_id=None, cloud_type="COMMUNITY"):
    """Create pod via REST."""
    body = {
        "name": name,
        "imageName": container_image,
        "gpuTypeIds": [gpu_type_id] if gpu_type_id else [],
        "gpuCount": gpu_count,
        "containerDiskInGb": disk_gb,
        "volumeInGb": 0 if not volume_mount_path else 50,
        "ports": ports,
        "env": env or {},
        "cloudType": cloud_type,
        "supportPublicIp": True,
        "interruptible": True,  # spot pricing
    }
    if volume_mount_path:
        body["volumeMountPath"] = volume_mount_path
    if network_volume_id:
        body["networkVolumeId"] = network_volume_id

    r = rest("POST", "/pods", json=body)
    print(f"  pod create response: {r.status_code} {r.text[:300]}")
    return r.json() if r.status_code == 200 else None

def pod_destroy(pod_id):
    r = rest("DELETE", f"/pods/{pod_id}")
    print(f"  pod {pod_id} destroyed: {r.status_code}")
    return r.status_code == 200

# ============================================================
# Pod specs for D + B + C
# ============================================================
POD_D = {
    "name": "thca-d-k2-star",
    "gpu_type_id": "NVIDIA RTX A4000",  # cheap GPU pod (~$0.20/hr), use CPU only
    "gpu_count": 1,
    "container_image": "runpod/pytorch:2.1.0-py3.10-cuda12.1.1-devel-ubuntu22.04",
    "disk_gb": 200,  # FASTQ + STAR index + outputs
    "task": "K2 STAR re-quantification (CPU-bound, GPU unused)",
    "expected_hr": 30,
    "expected_cost_usd": "5-10",
}

POD_C = {
    "name": "thca-c-alphafold",
    "gpu_type_id": "NVIDIA A100 80GB PCIe",  # $1.19/hr
    "gpu_count": 1,
    "container_image": "runpod/pytorch:2.1.0-py3.10-cuda12.1.1-devel-ubuntu22.04",
    "disk_gb": 100,
    "task": "AlphaFold HLA-DPB1*05:01 + TSHR/Tg peptide structure",
    "expected_hr": 10,
    "expected_cost_usd": "12-15",
}

POD_B = {
    "name": "thca-b-wsi-pathology",
    "gpu_type_id": "NVIDIA A40",  # $0.35/hr
    "gpu_count": 1,
    "container_image": "runpod/pytorch:2.1.0-py3.10-cuda12.1.1-devel-ubuntu22.04",
    "disk_gb": 200,  # WSI tiles
    "task": "TCGA WSI Hashimoto multimodal (PoC, 50 samples)",
    "expected_hr": 25,
    "expected_cost_usd": "10-15",
}

# ============================================================
# Master orchestration
# ============================================================
def main(mode="status"):
    print("=== RunPod Multi-Pod Orchestrator ===")
    print(f"  Spend limit: ${balance()['spendLimit']}/hr")
    print(f"  Balance: ${balance()['clientBalance']}")
    print(f"  Existing pods: {len(list_pods())}")

    if mode == "status":
        return

    if balance()["clientBalance"] < 20:
        print("\n❌ Balance < $20. User must add credit before proceeding.")
        print("   Recommended: $40 (D + C + B PoC, total ceiling)")
        print("   Min: $20 (D only)")
        return

    if mode == "spin_d":
        spec = POD_D
    elif mode == "spin_c":
        spec = POD_C
    elif mode == "spin_b":
        spec = POD_B
    elif mode == "spin_all":
        for s in [POD_D, POD_C, POD_B]:
            print(f"\n--- Spinning {s['name']} ---")
            pod = pod_create(
                name=s["name"],
                gpu_type_id=s["gpu_type_id"],
                gpu_count=s["gpu_count"],
                container_image=s["container_image"],
                disk_gb=s["disk_gb"],
                ports="22/tcp,8888/http",
                cloud_type="COMMUNITY",
            )
            if pod:
                print(f"  ✓ {s['name']}: {pod}")
        return
    else:
        print(f"unknown mode: {mode}")
        return

    print(f"\n--- Spinning {spec['name']} ---")
    print(f"  Expected: {spec['expected_hr']}hr, ~${spec['expected_cost_usd']}")
    pod = pod_create(
        name=spec["name"],
        gpu_type_id=spec["gpu_type_id"],
        gpu_count=spec["gpu_count"],
        container_image=spec["container_image"],
        disk_gb=spec["disk_gb"],
        ports="22/tcp,8888/http",
        cloud_type="COMMUNITY",
    )
    if pod:
        print(f"  ✓ {spec['name']} pod created: {json.dumps(pod, indent=2)}")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "status"
    main(mode)
