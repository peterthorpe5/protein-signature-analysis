"""Shell syntax and isolated SSH argument tests; no real server connection is made."""

import os
import subprocess
from pathlib import Path

import pytest

LAUNCHER = Path(__file__).parents[2] / "launch_protein_signature_denbi.sh"


@pytest.mark.parametrize(
    "arguments",
    [
        ["--help"],
        ["--version"],
        ["--local-port", "850"],
        ["--remote-port", "65536"],
        ["--config", "relative/path"],
        ["--ssh-host", "-bad"],
        ["--unknown"],
        ["--status", "--stop-tunnel"],
        ["--status", "--restart-app"],
        ["--environment"],
    ],
)
def test_launcher_named_option_validation(arguments):
    """Invalid options are rejected before any SSH connection is attempted."""
    result = subprocess.run(
        args=["bash", str(LAUNCHER), *arguments], capture_output=True, text=True, check=False
    )
    assert result.returncode == (0 if arguments[0] in {"--help", "--version"} else 2)


def test_launcher_preserves_remote_paths_and_dedicated_forward(tmp_path):
    """SSH setup uses quoted paths and the requested loopback ports, with no code upload."""
    binary = tmp_path / "bin"
    binary.mkdir()
    ssh = binary / "ssh"
    ssh.write_text(
        data=(
            "#!/usr/bin/env bash\n"
            'printf "%s\\n" "$@" >> "$PS_FAKE_LOG"\n'
            'if [[ " $* " == *" -O check "* ]]; then exit 1; fi\n'
            'if [[ "$1" == "-T" ]]; then cat > "$PS_FAKE_REMOTE"; fi\n'
            "exit 0\n"
        ),
        encoding="utf-8",
    )
    ssh.chmod(mode=0o755)
    curl = binary / "curl"
    curl.write_text(data="#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
    curl.chmod(mode=0o755)
    log, remote = tmp_path / "ssh_arguments.txt", tmp_path / "remote.sh"
    environment = {
        **os.environ,
        "PATH": f"{binary}:{os.environ['PATH']}",
        "PS_FAKE_LOG": str(log),
        "PS_FAKE_REMOTE": str(remote),
    }
    result = subprocess.run(
        args=[
            "bash",
            str(LAUNCHER),
            "--no-browser",
            "--local-port",
            "8521",
            "--remote-port",
            "8522",
            "--config",
            "/data/profile with spaces.toml",
        ],
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    text = log.read_text()
    assert "127.0.0.1:8521:127.0.0.1:8522" in text
    assert "profile\\ with\\ spaces.toml" in text
    assert "http://127.0.0.1:8521" in result.stdout
    assert subprocess.run(args=["bash", "-n", str(remote)], check=False).returncode == 0
    assert "python -" not in remote.read_text()
    assert "protein-signatures" in remote.read_text()
