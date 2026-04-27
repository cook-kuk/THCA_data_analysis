"""ClaudeCodeExecutor -- subprocess wrapper for shelling out to Claude Code.

Used when one agent dispatches subgoals to a separate Claude Code session
(the multi-terminal pattern from v3 -> v17p35).
"""
from __future__ import annotations

import asyncio
import json
import shutil
from dataclasses import dataclass


@dataclass
class ClaudeCodeRun:
    cmd: list[str]
    stdout: str
    stderr: str
    return_code: int


class ClaudeCodeExecutor:
    def __init__(self, claude_binary: str = "claude"):
        self.binary = shutil.which(claude_binary) or claude_binary

    async def dispatch(self, prompt: str, cwd: str | None = None, extra_args: list[str] | None = None) -> ClaudeCodeRun:
        cmd = [self.binary, "-p", prompt]
        if extra_args:
            cmd.extend(extra_args)
        proc = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            cwd=cwd,
        )
        stdout, stderr = await proc.communicate()
        return ClaudeCodeRun(
            cmd=cmd,
            stdout=stdout.decode(errors="ignore"),
            stderr=stderr.decode(errors="ignore"),
            return_code=proc.returncode or 0,
        )
