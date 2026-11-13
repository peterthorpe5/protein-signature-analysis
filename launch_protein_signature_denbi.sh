#!/usr/bin/env bash
# Launch the signature GUI and independent worker on Sabrina or Stefan through SSH.
# Bash 3.2+; Python belongs in the tested package, never inside this script.
set -Eeuo pipefail

PS_ACTION="start"
PS_SERVER="sabrina"
PS_HOST=""
PS_LOCAL_PORT=""
PS_REMOTE_PORT="8505"
PS_OPEN="yes"
PS_RESTART="no"
PS_ENV=""
PS_CONFIG=""

while (($#)); do
    case "$1" in
        --help|-h)
            cat <<'USAGE'
Run on your Mac after installing the update on the chosen server:
  bash launch_protein_signature_denbi.sh
  bash launch_protein_signature_denbi.sh --no-browser
  bash launch_protein_signature_denbi.sh --status
  bash launch_protein_signature_denbi.sh --stop-tunnel
  bash launch_protein_signature_denbi.sh --restart-app
  bash launch_protein_signature_denbi.sh --server stefan

Optional named settings:
  --server NAME        sabrina (default) or stefan.
  --ssh-host HOST      Existing SSH alias or user@host (denbi / stefan).
  --local-port PORT    Mac port (Sabrina 8510, Stefan 8511).
  --remote-port PORT   Server loopback port (default 8505).
  --environment PATH   Absolute signature Conda environment path on that server.
  --config PATH        Absolute workbench TOML path on that server.

The separate protein-signatures tmux session keeps GUI and worker alive after
disconnecting. Stopping the tunnel leaves jobs running. --restart-app affects
only this launcher's GUI; the worker continues. See docs/DENBI_WORKBENCH.md.
USAGE
            exit 0 ;;
        --version) printf 'Protein signature server launcher 1.1.0\n'; exit 0 ;;
        --no-browser) PS_OPEN="no" ;;
        --restart-app) PS_RESTART="yes" ;;
        --status|--stop-tunnel)
            if [[ "$PS_ACTION" != "start" ]]; then
                printf 'Choose one launcher action.\n' >&2; exit 2
            fi
            PS_ACTION="${1#--}" ;;
        --server|--ssh-host|--local-port|--remote-port|--environment|--config)
            if (($# < 2)) || [[ -z "$2" ]]; then
                printf '%s requires a value.\n' "$1" >&2; exit 2
            fi
            case "$1" in
                --server)
                    if [[ "$2" != "sabrina" && "$2" != "stefan" ]]; then
                        printf '%s\n' '--server must be sabrina or stefan.' >&2; exit 2
                    fi
                    PS_SERVER="$2" ;;
                --ssh-host)
                    if [[ ! "$2" =~ ^[A-Za-z0-9][A-Za-z0-9._@-]*$ ]]; then
                        printf 'Invalid SSH host/alias.\n' >&2; exit 2
                    fi
                    PS_HOST="$2" ;;
                --local-port|--remote-port)
                    if [[ ! "$2" =~ ^[1-9][0-9]{0,4}$ ]] || (( $2 < 1024 || $2 > 65535 )); then
                        printf '%s must be between 1024 and 65535.\n' "$1" >&2; exit 2
                    fi
                    if [[ "$1" == "--local-port" ]]; then PS_LOCAL_PORT="$2"; else PS_REMOTE_PORT="$2"; fi ;;
                --environment|--config)
                    if [[ "$2" != /* || "$2" == *$'\n'* || "$2" == *$'\r'* ]]; then
                        printf '%s requires an absolute remote path.\n' "$1" >&2; exit 2
                    fi
                    if [[ "$1" == "--environment" ]]; then PS_ENV="$2"; else PS_CONFIG="$2"; fi ;;
            esac
            shift ;;
        *) printf 'Unknown option: %s. Use --help.\n' "$1" >&2; exit 2 ;;
    esac
    shift
done

if [[ "$PS_SERVER" == "stefan" ]]; then
    PS_HOST="${PS_HOST:-stefan}"
    PS_LOCAL_PORT="${PS_LOCAL_PORT:-8511}"
    PS_ENV="${PS_ENV:-/data/pete/software/miniforge3/envs/protein_signature_analysis}"
    PS_CONFIG="${PS_CONFIG:-/data/pete/protein_signature_research/config/workbench.toml}"
else
    PS_HOST="${PS_HOST:-denbi}"
    PS_LOCAL_PORT="${PS_LOCAL_PORT:-8510}"
    PS_ENV="${PS_ENV:-/home/ubuntu/data/pete/software/miniforge3/envs/protein_signature_analysis}"
    PS_CONFIG="${PS_CONFIG:-/home/ubuntu/data/protein_signature_research/config/workbench.toml}"
fi

if [[ "$PS_ACTION" != "start" && "$PS_RESTART" == "yes" ]]; then
    printf '%s\n' '--restart-app applies only to a normal launch.' >&2; exit 2
fi
for PS_BINARY in ssh curl; do
    if ! command -v "$PS_BINARY" >/dev/null 2>&1; then
        printf 'Required command not found: %s\n' "$PS_BINARY" >&2; exit 1
    fi
done
# Hash the host/path configuration to avoid cross-connecting other deployments.
if command -v shasum >/dev/null 2>&1; then
    PS_ID="$(printf '%s\n' "$PS_HOST" "$PS_ENV" "$PS_CONFIG" | shasum -a 256)"
else
    PS_ID="$(printf '%s\n' "$PS_HOST" "$PS_ENV" "$PS_CONFIG" | sha256sum)"
fi
PS_SOCKET="$HOME/.ssh/ps_${PS_ID:0:12}_${PS_LOCAL_PORT}_${PS_REMOTE_PORT}.sock"
PS_URL="http://127.0.0.1:$PS_LOCAL_PORT"

if [[ "$PS_ACTION" == "stop-tunnel" ]]; then
    if ssh -F /dev/null -S "$PS_SOCKET" -O check "$PS_HOST" >/dev/null 2>&1; then
        ssh -F /dev/null -S "$PS_SOCKET" -O exit "$PS_HOST"
        printf 'Tunnel disconnected. Server jobs remain running.\n'
    else
        printf 'No active tunnel owned by this launcher.\n'
    fi
    exit 0
fi

printf -v PS_REMOTE_COMMAND 'bash -s -- %q %q %q %q %q' \
    "$PS_ACTION" "$PS_REMOTE_PORT" "$PS_RESTART" "$PS_ENV" "$PS_CONFIG"
ssh -T -o ClearAllForwardings=yes -o ConnectTimeout=15 \
    -o ServerAliveInterval=15 -o ServerAliveCountMax=3 "$PS_HOST" "$PS_REMOTE_COMMAND" <<'PS_REMOTE'
set -Eeuo pipefail
PS_MODE="$1"
PS_PORT="$2"
PS_RESTART_GUI="$3"
PS_PREFIX="$4"
PS_PROFILE="$5"
PS_SESSION="protein-signatures"
PS_TMUX="$PS_PREFIX/bin/tmux"
PS_APP="$PS_PREFIX/bin/protein-signature-app"
PS_WORKER="$PS_PREFIX/bin/protein-signature-workbench"
PS_HEALTH="http://127.0.0.1:$PS_PORT/_stcore/health"
PS_APP_TARGET="$PS_SESSION:=app"
PS_WORKER_TARGET="$PS_SESSION:=worker"

for PS_BINARY in "$PS_TMUX" "$PS_APP" "$PS_WORKER"; do
    if [[ ! -x "$PS_BINARY" ]]; then
        printf 'Missing executable: %s. Complete the server installation first.\n' "$PS_BINARY" >&2
        exit 1
    fi
done
if [[ ! -r "$PS_PROFILE" ]]; then
    printf 'Missing profile: %s. Configure the installed release first.\n' "$PS_PROFILE" >&2; exit 1
fi
if ! command -v curl >/dev/null 2>&1; then
    printf 'curl is required on the server.\n' >&2; exit 1
fi
PS_WORKSPACE="$("$PS_WORKER" --config "$PS_PROFILE" --action workspace)"
if [[ "$PS_WORKSPACE" != /* || "$PS_WORKSPACE" == *$'\n'* ]]; then
    printf 'The profile must report one absolute workspace path.\n' >&2; exit 1
fi
PS_LOG_ROOT="$PS_WORKSPACE/logs"
PS_TEMP="$PS_WORKSPACE/temporary"
PS_CACHE="$PS_WORKSPACE/cache"
mkdir -p "$PS_LOG_ROOT" "$PS_TEMP" "$PS_CACHE"
export TMPDIR="$PS_TEMP" XDG_CACHE_HOME="$PS_CACHE" PYTHONUNBUFFERED=1
# tmux may have been started before a Conda environment was activated.
export PATH="$PS_PREFIX/bin:$PATH"

ps_pane_state() {
    # Return the state of exactly one owned service pane.
    local ps_value
    if ! ps_value="$("$PS_TMUX" list-panes -t "$1" -F '#{pane_dead}' 2>/dev/null)"; then
        printf 'absent\n'; return 0
    fi
    case "$ps_value" in
        0|1) printf '%s\n' "$ps_value" ;;
        *) printf 'Expected one pane in %s. Inspect the tmux session.\n' "$1" >&2; return 1 ;;
    esac
}

ps_service() {
    # Create/revive one owned service; preserve any live process.
    local ps_target="$1" ps_name="$2" ps_command="$3" ps_state
    ps_state="$(ps_pane_state "$ps_target")"
    if [[ "$ps_state" == "1" ]]; then
        "$PS_TMUX" respawn-pane -t "$ps_target" "$ps_command"
    elif [[ "$ps_state" == "absent" ]]; then
        if "$PS_TMUX" has-session -t "$PS_SESSION" 2>/dev/null; then
            "$PS_TMUX" new-window -d -t "$PS_SESSION:" -n "$ps_name" "$ps_command"
        else
            "$PS_TMUX" new-session -d -s "$PS_SESSION" -n "$ps_name" "$ps_command"
        fi
    fi
}

ps_verify_owner() {
    # Refuse to reuse a live service belonging to a different executable/profile.
    local ps_started ps_binary_quoted ps_profile_quoted
    ps_started="$("$PS_TMUX" list-panes -t "$1" -F '#{pane_start_command}')"
    printf -v ps_binary_quoted '%q' "$2"
    printf -v ps_profile_quoted '%q' "$PS_PROFILE"
    if [[ "$ps_started" != *"$2"* && "$ps_started" != *"$ps_binary_quoted"* ]] \
        || [[ "$ps_started" != *"$PS_PROFILE"* && "$ps_started" != *"$ps_profile_quoted"* ]]; then
        printf 'Live service in %s uses another executable/profile; inspect it before proceeding.\n' "$1" >&2
        return 1
    fi
}

if [[ "$(ps_pane_state "$PS_APP_TARGET")" == "0" ]]; then
    ps_verify_owner "$PS_APP_TARGET" "$PS_APP"
fi
if [[ "$(ps_pane_state "$PS_WORKER_TARGET")" == "0" ]]; then
    ps_verify_owner "$PS_WORKER_TARGET" "$PS_WORKER"
fi

if [[ "$PS_MODE" == "status" ]]; then
    "$PS_TMUX" list-windows -t "$PS_SESSION"
    "$PS_WORKER" --config "$PS_PROFILE" --action status
    [[ "$(ps_pane_state "$PS_APP_TARGET")" == "0" ]]
    [[ "$(ps_pane_state "$PS_WORKER_TARGET")" == "0" ]]
    curl --fail --silent --show-error --max-time 5 "$PS_HEALTH"
    printf '\n'; exit 0
fi

"$PS_WORKER" --config "$PS_PROFILE" --action doctor
PS_APP_STATE="$(ps_pane_state "$PS_APP_TARGET")"
if [[ "$PS_APP_STATE" == "0" && "$PS_RESTART_GUI" != "yes" ]]; then
    PS_STARTED="$("$PS_TMUX" list-panes -t "$PS_APP_TARGET" -F '#{pane_start_command}')"
    if [[ "$PS_STARTED" != *"--port $PS_PORT "* ]]; then
        printf 'The owned GUI uses another port. Use --restart-app with the requested port.\n' >&2
        exit 1
    fi
fi
if [[ "$PS_APP_STATE" == "0" && "$PS_RESTART_GUI" == "yes" ]]; then
    ps_verify_owner "$PS_APP_TARGET" "$PS_APP"
    "$PS_TMUX" kill-pane -t "$PS_APP_TARGET"
fi

printf -v PS_GUI_COMMAND '%q ' env "PATH=$PATH" "TMPDIR=$TMPDIR" \
    "XDG_CACHE_HOME=$XDG_CACHE_HOME" PYTHONUNBUFFERED=1 \
    "$PS_APP" --workbench-config "$PS_PROFILE" --port "$PS_PORT" --address 127.0.0.1
printf -v PS_GUI_LOG '%q' "$PS_LOG_ROOT/app.log"
PS_GUI_COMMAND="$PS_GUI_COMMAND >>$PS_GUI_LOG 2>&1"
printf -v PS_WORKER_COMMAND '%q ' env "PATH=$PATH" "TMPDIR=$TMPDIR" \
    "XDG_CACHE_HOME=$XDG_CACHE_HOME" PYTHONUNBUFFERED=1 \
    "$PS_WORKER" --config "$PS_PROFILE" --action worker
printf -v PS_WORKER_LOG '%q' "$PS_LOG_ROOT/worker.log"
PS_WORKER_COMMAND="$PS_WORKER_COMMAND >>$PS_WORKER_LOG 2>&1"
ps_service "$PS_WORKER_TARGET" worker "$PS_WORKER_COMMAND"
ps_service "$PS_APP_TARGET" app "$PS_GUI_COMMAND"

PS_READY="no"
for ((PS_ATTEMPT=0; PS_ATTEMPT<45; PS_ATTEMPT++)); do
    if [[ "$(ps_pane_state "$PS_APP_TARGET")" != "0" ]]; then break; fi
    if curl --fail --silent --max-time 1 "$PS_HEALTH" >/dev/null 2>&1; then
        PS_READY="yes"; break
    fi
    sleep 1
done
if [[ "$PS_READY" != "yes" || "$(ps_pane_state "$PS_WORKER_TARGET")" != "0" ]]; then
    tail -n 60 "$PS_LOG_ROOT/app.log" "$PS_LOG_ROOT/worker.log" >&2 || true
    printf 'Server services did not become ready. Check logs above.\n' >&2; exit 1
fi
printf 'Protein signature app and worker are running on %s.\n' "$(hostname)"
PS_REMOTE

if [[ "$PS_ACTION" == "status" ]]; then
    ssh -F /dev/null -S "$PS_SOCKET" -O check "$PS_HOST"
else
    if ! ssh -F /dev/null -S "$PS_SOCKET" -O check "$PS_HOST" >/dev/null 2>&1; then
        if [[ -S "$PS_SOCKET" ]]; then rm -- "$PS_SOCKET"; fi
        ssh -fNT -M -S "$PS_SOCKET" -o ControlPersist=yes -o ClearAllForwardings=yes \
            -o ConnectTimeout=15 -o ServerAliveInterval=15 -o ServerAliveCountMax=3 "$PS_HOST"
    fi
    ssh -F /dev/null -S "$PS_SOCKET" -O forward -o ExitOnForwardFailure=yes \
        -L "127.0.0.1:$PS_LOCAL_PORT:127.0.0.1:$PS_REMOTE_PORT" "$PS_HOST"
fi
curl --fail --silent --show-error --max-time 3 --retry 5 --retry-connrefused \
    --retry-delay 1 --retry-max-time 15 "$PS_URL/_stcore/health" >/dev/null
printf '\nProtein signature workbench: %s\n' "$PS_URL"
if [[ "$PS_ACTION" == "start" && "$PS_OPEN" == "yes" && "$(uname -s)" == "Darwin" ]]; then
    open "$PS_URL" || printf 'Open the URL above in your browser.\n' >&2
fi
