#!/usr/bin/env sh
# ---------------------------------------------------------------------------
# Portable Python launcher for Claude Code hooks.
#
# Why this exists: hooks that hardcode `conda run -n <env>` break on any other
# machine/fork that does not have an env by that exact name (e.g. the original
# template shipped `-n [YOUR_CONDA_ENV]`, which fails with
# EnvironmentLocationNotFound). This launcher resolves an interpreter WITHOUT
# committing a machine-specific env name, so the shared config is universal.
#
# Usage:  run-python.sh <script.py> [args...]
#
# Resolution order (first that works wins):
#   1. $CLAUDE_HOOK_PYTHON  -- explicit override; may be a full command, e.g.
#                             "conda run -n myenv python" or "/usr/bin/python3"
#   2. conda env named in $CONDA_DEFAULT_ENV  (only if conda + that env exist)
#   3. python3 / python on PATH
#   4. nothing found -> exit 0 (a hook must never block on a missing interpreter)
#
# To pin a specific interpreter on one machine without touching committed
# config, set CLAUDE_HOOK_PYTHON or CONDA_DEFAULT_ENV in the gitignored
# .claude/settings.local.json "env" block.
# ---------------------------------------------------------------------------
SCRIPT="$1"
[ -n "$SCRIPT" ] || exit 0
shift

# 1. Explicit override (intentionally unquoted so it can be a multi-word command)
if [ -n "$CLAUDE_HOOK_PYTHON" ]; then
  exec $CLAUDE_HOOK_PYTHON "$SCRIPT" "$@"
fi

# 2. Configured conda env, only if it actually exists
if [ -n "$CONDA_DEFAULT_ENV" ] && command -v conda >/dev/null 2>&1 \
   && conda run -n "$CONDA_DEFAULT_ENV" true >/dev/null 2>&1; then
  exec conda run --no-capture-output -n "$CONDA_DEFAULT_ENV" python "$SCRIPT" "$@"
fi

# 3. Plain Python on PATH (universal default)
if command -v python3 >/dev/null 2>&1; then exec python3 "$SCRIPT" "$@"; fi
if command -v python  >/dev/null 2>&1; then exec python  "$SCRIPT" "$@"; fi

# 4. No interpreter available -- do not block the hook
exit 0
