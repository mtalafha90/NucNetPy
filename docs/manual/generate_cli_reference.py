"""Regenerate reference/cli.md from the command-line parser.

The page is the parser's own help text, so it can only match the code if it
is rebuilt whenever an option changes:

    python docs/manual/generate_cli_reference.py

``make -C docs/manual pdf`` and ``html`` run this first.  Use Python 3.10,
3.11 or 3.12: argparse words its help differently on 3.9 and 3.13.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

INTRO = """# Command line

Installing the package puts a `nucnetpy` command on the path. This page is
generated from the argument parser by `docs/manual/generate_cli_reference.py`,
so it matches the code.
"""


def main() -> None:
    os.environ["COLUMNS"] = "80"          # argparse wraps to the terminal width
    from nucnetpy.cli import build_parser

    parser = build_parser()
    commands = next(a for a in parser._actions
                    if isinstance(a, argparse._SubParsersAction)).choices
    parts = [INTRO, "```\n" + parser.format_help().rstrip() + "\n```\n"]
    for name, sub in commands.items():
        parts.append(f"## `{name}`\n\n```\n{sub.format_help().rstrip()}\n```\n")
    target = Path(__file__).parent / "reference" / "cli.md"
    target.write_text("\n".join(parts))
    print(f"wrote {target} ({len(commands)} commands)")


if __name__ == "__main__":
    main()
