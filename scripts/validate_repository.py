"""Run dependency-free structural checks before publishing the repository."""

from __future__ import annotations

import ast
import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
ABSOLUTE_PATH = re.compile(
    r"(?:(?<![A-Za-z0-9_])[A-Za-z]:[/\\]|/(?:home|Users|Volumes|mnt)/)"
)


def add_error(errors: list[str], path: Path, message: str) -> None:
    errors.append(f"{path.relative_to(ROOT)}: {message}")


def validate_r_delimiters(source: str, path: Path, label: str, errors: list[str]) -> None:
    """Check balanced R delimiters while ignoring strings and comments."""
    opening = {"(": ")", "[": "]", "{": "}"}
    closing = {value: key for key, value in opening.items()}
    stack: list[tuple[str, int]] = []
    quote = None
    escaped = False
    in_comment = False
    line_number = 1

    for character in source:
        if character == "\n":
            line_number += 1
            in_comment = False
            continue
        if in_comment:
            continue
        if escaped:
            escaped = False
            continue
        if quote is not None:
            if character == "\\":
                escaped = True
            elif character == quote:
                quote = None
            continue
        if character in ("\"", "'", "`"):
            quote = character
        elif character == "#":
            in_comment = True
        elif character in opening:
            stack.append((character, line_number))
        elif character in closing:
            if not stack or stack[-1][0] != closing[character]:
                add_error(errors, path, f"{label} has an unmatched {character} near line {line_number}")
                return
            stack.pop()

    if quote is not None:
        add_error(errors, path, f"{label} has an unterminated quoted string")
    elif stack:
        character, opening_line = stack[-1]
        add_error(errors, path, f"{label} has an unclosed {character} from line {opening_line}")


def validate_notebook(path: Path, errors: list[str]) -> None:
    try:
        notebook = json.loads(path.read_text(encoding="utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        add_error(errors, path, f"invalid notebook JSON ({error})")
        return

    kernel = notebook.get("metadata", {}).get("kernelspec", {}).get("name")
    for cell_number, cell in enumerate(notebook.get("cells", []), start=1):
        source = "".join(cell.get("source", []))
        if ABSOLUTE_PATH.search(source):
            add_error(errors, path, f"cell {cell_number} contains an absolute path")

        if cell.get("cell_type") != "code":
            continue

        if cell.get("outputs"):
            add_error(errors, path, f"cell {cell_number} contains saved output")
        if cell.get("execution_count") is not None:
            add_error(errors, path, f"cell {cell_number} has an execution counter")
        if any(ord(character) > 127 for character in source):
            add_error(errors, path, f"cell {cell_number} contains a non-ASCII code/comment character")

        if kernel == "ir":
            validate_r_delimiters(source, path, f"cell {cell_number}", errors)

        if kernel == "python3" and source.strip():
            # IPython line magics are valid in a notebook but not in Python's
            # standard abstract-syntax-tree parser.
            python_source = "\n".join(
                line for line in source.splitlines()
                if not line.lstrip().startswith(("%", "!"))
            )
            try:
                ast.parse(python_source)
            except SyntaxError as error:
                add_error(errors, path, f"cell {cell_number} has invalid Python syntax ({error.msg})")


def validate_text_source(path: Path, errors: list[str]) -> None:
    try:
        source = path.read_text(encoding="utf-8")
    except UnicodeDecodeError as error:
        add_error(errors, path, f"is not valid UTF-8 ({error})")
        return

    if ABSOLUTE_PATH.search(source):
        add_error(errors, path, "contains an absolute local path")
    if any(ord(character) > 127 for character in source):
        add_error(errors, path, "contains a non-ASCII code/comment character")
    if path.is_relative_to(ANALYSIS) and re.search(r"\bsetwd\s*\(", source):
        add_error(errors, path, "changes the working directory")
    if path.is_relative_to(ANALYSIS) and re.search(
        r"\b(?:install\.packages|install_github|py_install)\s*\(", source
    ):
        add_error(errors, path, "installs software during an analysis")
    if path.suffix.lower() == ".r":
        validate_r_delimiters(source, path, "R source", errors)


def main() -> int:
    errors: list[str] = []
    required_paths = (
        ROOT / "README.md",
        ROOT / "config" / "sample_manifest.csv",
        ROOT / "data" / "README.md",
        ROOT / "docs" / "analysis_workflow.md",
        ROOT / "docs" / "reproducibility_checklist.md",
        ROOT / "environment" / "README.md",
        ROOT / "analysis" / "01_preprocessing" / "01_preprocess_scrna_samples.R",
    )
    for required_path in required_paths:
        if not required_path.exists():
            add_error(errors, required_path, "required repository file is missing")

    for path in sorted(ANALYSIS.rglob("*")):
        if not path.is_file():
            continue
        if path.suffix == ".ipynb":
            validate_notebook(path, errors)
        elif path.suffix.lower() in {".r", ".py", ".sh"}:
            validate_text_source(path, errors)

    if errors:
        print(f"Repository validation failed with {len(errors)} issue(s):")
        for error in errors:
            print(f"- {error}")
        return 1

    notebook_count = len(list(ANALYSIS.rglob("*.ipynb")))
    source_count = sum(
        1 for path in ANALYSIS.rglob("*")
        if path.is_file() and path.suffix.lower() in {".r", ".py", ".sh"}
    )
    print(
        f"Repository validation passed: {source_count} scripts and "
        f"{notebook_count} notebooks checked."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
