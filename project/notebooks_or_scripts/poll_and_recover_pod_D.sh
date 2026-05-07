#!/bin/bash
# poll_and_recover_pod_D.sh — Pod D (K2 STAR) monitor + recovery
# Polls every 30 min via RunPod API + SSH probe. When run.log shows ALL DONE
# (or expected output files exist), SCPs results back and runs post_pod_D.
# Status file: /tmp/pod_D_status.txt  Log: /tmp/pod_D_recovery.log

set -u
PROJ="/home/seungho/personal/THCA_data_analysis"
LOG="/tmp/pod_D_recovery.log"
STATUS="/tmp/pod_D_status.txt"
RESULTS_DIR="$PROJ/project/results/d7p3_k2_starred"
POD_ID="xgk5y4oqwvo7te"
POLL_INTERVAL=1800   # 30 min
SSH_CONN_TIMEOUT=60  # banner exchange may be slow under STAR load
MAX_HOURS=30         # auto-shutdown 24hr, give 30hr cap
START_TIME=$(date +%s)

mkdir -p "$RESULTS_DIR"
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }
write_status() { echo "$1" > "$STATUS"; }

# Get fresh SSH endpoint via RunPod API (port mapping may change)
get_ssh_endpoint() {
    "$PROJ/.venv/bin/python3" - <<PY
import os, requests, json, sys
KEY = open(os.path.expanduser('~/.runpod/config')).read().split('=',1)[1].strip()
q = """query { pod(input: {podId: "$POD_ID"}) { desiredStatus runtime { ports { ip publicPort privatePort isIpPublic type } } }}"""
try:
    r = requests.post('https://api.runpod.io/graphql',
        headers={'Authorization': f'Bearer {KEY}', 'Content-Type': 'application/json'},
        json={'query': q}, timeout=20)
    pod = r.json().get('data', {}).get('pod')
    if not pod:
        print("DEAD"); sys.exit(0)
    if pod.get('desiredStatus') != 'RUNNING':
        print(f"NOT_RUNNING:{pod.get('desiredStatus')}"); sys.exit(0)
    rt = pod.get('runtime') or {}
    for p in (rt.get('ports') or []):
        if p.get('privatePort') == 22 and p.get('isIpPublic'):
            print(f"{p['ip']} {p['publicPort']}"); sys.exit(0)
    print("NO_SSH_PORT"); sys.exit(0)
except Exception as e:
    print(f"API_ERR:{e}"); sys.exit(0)
PY
}

ssh_pod() {
    local ip=$1 port=$2 cmd=$3
    timeout $((SSH_CONN_TIMEOUT + 30)) ssh \
        -o BatchMode=yes \
        -o StrictHostKeyChecking=accept-new \
        -o UserKnownHostsFile=$HOME/.ssh/known_hosts \
        -o ConnectTimeout=$SSH_CONN_TIMEOUT \
        -o ServerAliveInterval=15 \
        -o ServerAliveCountMax=4 \
        -p "$port" "root@$ip" "$cmd" 2>&1
}

scp_pod() {
    local ip=$1 port=$2 src=$3 dst=$4
    timeout 600 scp \
        -o BatchMode=yes \
        -o StrictHostKeyChecking=accept-new \
        -o UserKnownHostsFile=$HOME/.ssh/known_hosts \
        -o ConnectTimeout=$SSH_CONN_TIMEOUT \
        -o ServerAliveInterval=15 \
        -P "$port" "root@$ip:$src" "$dst" 2>&1
}

write_status "polling"
log "starting Pod D monitor (pod=$POD_ID, interval=${POLL_INTERVAL}s, max ${MAX_HOURS}h)"

ITER=0
while true; do
    ITER=$((ITER + 1))
    NOW=$(date +%s)
    ELAPSED_HR=$(( (NOW - START_TIME) / 3600 ))

    if [ "$ELAPSED_HR" -ge "$MAX_HOURS" ]; then
        log "monitor cap reached (${MAX_HOURS}h) — exiting without completion"
        write_status "timeout"
        exit 2
    fi

    EP=$(get_ssh_endpoint)
    log "iter=$ITER endpoint=$EP elapsed=${ELAPSED_HR}h"

    case "$EP" in
        DEAD|NOT_RUNNING:*)
            log "pod no longer running: $EP"
            write_status "pod_gone:$EP"
            exit 3
            ;;
        API_ERR:*|NO_SSH_PORT)
            log "transient API issue: $EP — sleep $POLL_INTERVAL and retry"
            sleep $POLL_INTERVAL
            continue
            ;;
    esac

    IP=$(echo "$EP" | awk '{print $1}')
    PORT=$(echo "$EP" | awk '{print $2}')

    PROBE=$(ssh_pod "$IP" "$PORT" "
        tail -40 /workspace/run.log 2>/dev/null
        echo '===STAGE_FILES==='
        ls -la /workspace/k2_star/results/k2_counts.tsv 2>/dev/null
        ls -la /workspace/k2_star/aligned/*.bam 2>/dev/null | head -5
        ls -d /workspace/k2_star/star_index/SAindex 2>/dev/null
        echo '===END==='
        echo === ALL DONE marker check ===
        grep -c 'ALL DONE' /workspace/run.log 2>/dev/null || echo 0
    ")
    PROBE_RC=$?

    if [ $PROBE_RC -ne 0 ]; then
        log "SSH probe failed (rc=$PROBE_RC) — pod may be busy with STAR I/O"
        log "  excerpt: $(echo "$PROBE" | tail -3 | tr '\n' ' ')"
        sleep $POLL_INTERVAL
        continue
    fi

    log "SSH probe OK"
    echo "$PROBE" | sed 's/^/  /' | tee -a "$LOG" > /dev/null

    if echo "$PROBE" | grep -q "ALL DONE"; then
        log "ALL DONE detected — proceeding with recovery"
        break
    fi

    if echo "$PROBE" | grep -qE "(Traceback|FATAL|Killed|Cannot allocate|No space left)"; then
        log "POD-SIDE ERROR detected in run.log:"
        echo "$PROBE" | grep -E "(Traceback|FATAL|Killed|Cannot allocate|No space left)" | tee -a "$LOG"
        write_status "pod_error_continuing"
        # Continue polling — user may decide to abort manually
    fi

    log "still running — sleep $POLL_INTERVAL"
    sleep $POLL_INTERVAL
done

# === Recovery phase ===
write_status "recovering"
log "recovery: SCP results from $IP:$PORT"

scp_pod "$IP" "$PORT" "/workspace/k2_star/results/k2_counts.tsv" "$RESULTS_DIR/k2_counts.tsv" \
    | tee -a "$LOG" >/dev/null && log "  ✓ k2_counts.tsv" || log "  ✗ k2_counts.tsv FAIL"

scp_pod "$IP" "$PORT" "/workspace/run.log" "$RESULTS_DIR/run.log" \
    | tee -a "$LOG" >/dev/null && log "  ✓ run.log" || log "  ✗ run.log FAIL"

# Try to grab any per-sample STAR Log.final.out for QC
ssh_pod "$IP" "$PORT" "tar czf /workspace/star_logs.tar.gz /workspace/k2_star/aligned/*Log.final.out 2>/dev/null && echo TARRED" >> "$LOG" 2>&1
scp_pod "$IP" "$PORT" "/workspace/star_logs.tar.gz" "$RESULTS_DIR/star_logs.tar.gz" \
    | tee -a "$LOG" >/dev/null && log "  ✓ star_logs.tar.gz" || log "  ✗ star_logs.tar.gz (optional, ok if missing)"

# Try featureCounts summary (column counts QC)
scp_pod "$IP" "$PORT" "/workspace/k2_star/results/k2_counts.tsv.summary" "$RESULTS_DIR/k2_counts.tsv.summary" \
    | tee -a "$LOG" >/dev/null && log "  ✓ counts.summary" || log "  (counts.summary not present)"

# Run post-processing
log "running post_pod_D_meta3cohort.py"
cd "$PROJ"
"$PROJ/.venv/bin/python3" "$PROJ/project/notebooks_or_scripts/post_pod_D_meta3cohort.py" 2>&1 \
    | tee -a "$LOG" \
    | tee "$RESULTS_DIR/post_pod_D_stdout.log" >/dev/null

POST_RC=${PIPESTATUS[0]}
if [ $POST_RC -eq 0 ]; then
    log "post_pod_D_meta3cohort.py succeeded"
    write_status "completed"
else
    log "post_pod_D_meta3cohort.py FAILED (rc=$POST_RC) — raw counts still SCP'd"
    write_status "completed_post_failed:rc=$POST_RC"
fi

log "Pod D recovery done"
exit 0
