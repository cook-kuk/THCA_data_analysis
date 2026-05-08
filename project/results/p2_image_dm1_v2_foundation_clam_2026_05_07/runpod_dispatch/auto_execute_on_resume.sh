#!/bin/bash
# Pod resume 즉시 모든 단계 자동 실행:
#   1. runpodctl receive 0337-model-middle-finish-17  (3.4 GB tar)
#   2. tar xzf + cd workspace
#   3. env_setup.sh (CUDA + UNI 의존성)
#   4. HF_TOKEN 검증 (없으면 멈춤)
#   5. run_phase1.sh (UNI embedding + correlation)
#   6. PHASE1_REPORT.md 결과 받아오기
#   7. PASS 면 run_phase2.sh
#   8. PHASE2_REPORT.md 받아오기
#
# 모든 출력은 /tmp/p2_auto_execute.log

set +e
exec >> /tmp/p2_auto_execute.log 2>&1

POD_ID="50uem6t82i0t9g"
SEND_CODE="0337-model-middle-finish-17"
SLEEP_SEC=30

export RUNPOD_API_KEY=rpa_H807GFVSON867YLORS81R2GG1PVBOAW799FUSHQM15yg3o
export PATH="/tmp:$PATH"

echo ""
echo "==============================================="
echo "[$(date)] AUTO-EXECUTE WATCHER STARTED"
echo "  Pod: $POD_ID"
echo "  Code: $SEND_CODE"
echo "==============================================="

# ============================================================
# Step 1 — wait until pod is RUNNING + SSH info available
# ============================================================
echo "[wait] polling pod status..."
while true; do
    INFO=$(/tmp/runpodctl ssh info "$POD_ID" 2>&1)
    if echo "$INFO" | grep -q '"status".*RUNNING' && echo "$INFO" | grep -q '"port"'; then
        echo "[ok] pod is RUNNING with SSH info"
        echo "$INFO" | head -30
        break
    fi
    STATUS=$(echo "$INFO" | grep -oE '"status"[^,]*' | head -1)
    echo "  [$(date +%H:%M:%S)] not ready: $STATUS"
    sleep $SLEEP_SEC
done

# extract SSH host + port
SSH_HOST=$(echo "$INFO" | python3 -c "import json,sys;d=json.load(sys.stdin);
ports=d.get('ports',d.get('publicPorts',[]));
for p in ports:
    if p.get('privatePort')==22 or p.get('proxyType')=='ssh':
        print(p.get('ip','')); break" 2>/dev/null)
SSH_PORT=$(echo "$INFO" | python3 -c "import json,sys;d=json.load(sys.stdin);
ports=d.get('ports',d.get('publicPorts',[]));
for p in ports:
    if p.get('privatePort')==22 or p.get('proxyType')=='ssh':
        print(p.get('publicPort','')); break" 2>/dev/null)
echo "[ssh] host=$SSH_HOST port=$SSH_PORT"

# fallback to runpod proxy if direct fails
SSH_PROXY="${POD_ID}@ssh.runpod.io"

ssh_cmd() {
    local cmd="$1"
    if [[ -n "$SSH_HOST" && -n "$SSH_PORT" ]]; then
        ssh -tt -i ~/.ssh/id_ed25519 \
            -o StrictHostKeyChecking=accept-new \
            -o ConnectTimeout=15 \
            -p "$SSH_PORT" "root@$SSH_HOST" "$cmd"
    else
        ssh -tt -i ~/.ssh/id_ed25519 \
            -o StrictHostKeyChecking=accept-new \
            -o ConnectTimeout=15 \
            "$SSH_PROXY" "$cmd"
    fi
}

# wait additional 60s for sshd to be fully up
echo "[wait] +60s for sshd..."
sleep 60

# verify SSH alive
echo "[verify] ssh alive?"
SSH_TEST=$(ssh_cmd "echo ALIVE; nvidia-smi --query-gpu=name --format=csv,noheader" 2>&1)
echo "$SSH_TEST" | head -5
if ! echo "$SSH_TEST" | grep -q ALIVE; then
    echo "[ERROR] SSH not working — abort"
    exit 1
fi

# ============================================================
# Step 2 — install runpodctl on pod + receive
# ============================================================
echo ""
echo "==============================================="
echo "[$(date)] STEP 2: install runpodctl + receive"
echo "==============================================="
ssh_cmd "
cd /workspace || cd ~
which runpodctl || (
    wget -q -O /usr/local/bin/runpodctl https://github.com/runpod/runpodctl/releases/latest/download/runpodctl-linux-amd64
    chmod +x /usr/local/bin/runpodctl
)
runpodctl --version
echo 'starting receive...'
runpodctl receive $SEND_CODE
ls -lh p2_image_dm1_v2.tar.gz
" 2>&1 | tee -a /tmp/p2_step2.log

# ============================================================
# Step 3 — extract + env setup
# ============================================================
echo ""
echo "==============================================="
echo "[$(date)] STEP 3: tar extract + env_setup"
echo "==============================================="
ssh_cmd "
cd /workspace
tar xzf p2_image_dm1_v2.tar.gz
cd project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/
bash runpod_dispatch/env_setup.sh
" 2>&1 | tee -a /tmp/p2_step3.log

# ============================================================
# Step 4 — Phase 1
# ============================================================
echo ""
echo "==============================================="
echo "[$(date)] STEP 4: Phase 1 (UNI embedding + correlation)"
echo "==============================================="
if [[ -z "${HF_TOKEN:-}" ]]; then
    echo "[WARN] HF_TOKEN not set in environment — UNI gated; pod side may fall back to ImageNet ViT-L"
    echo "[WARN] To fix: export HF_TOKEN=hf_xxx before re-running this script"
fi

ssh_cmd "
cd /workspace/project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/
export HF_TOKEN='${HF_TOKEN:-}'
bash runpod_dispatch/run_phase1.sh
" 2>&1 | tee -a /tmp/p2_step4.log

# pull back Phase 1 results
echo "[pull] Phase 1 results back..."
ssh_cmd "
cd /workspace/project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/
tar czf /tmp/phase1_results.tar.gz phase1_gse250521/
" 2>&1 | head -5

# scp back (use tar pipe to avoid intermediate file)
if [[ -n "$SSH_HOST" && -n "$SSH_PORT" ]]; then
    scp -i ~/.ssh/id_ed25519 -P "$SSH_PORT" \
        "root@$SSH_HOST:/tmp/phase1_results.tar.gz" \
        /data/phase1_results.tar.gz
fi

if [[ -f /data/phase1_results.tar.gz ]]; then
    cd /home/seungho/personal/THCA_data_analysis
    tar xzf /data/phase1_results.tar.gz -C project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/
    echo "[ok] phase1 results extracted locally"
    cat project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase1_gse250521/PHASE1_REPORT.md 2>/dev/null | head -30
fi

echo ""
echo "==============================================="
echo "[$(date)] AUTO-EXECUTE WATCHER COMPLETE"
echo "==============================================="
echo "Phase 1 done. Phase 2 launch decision: review PHASE1_REPORT.md."
