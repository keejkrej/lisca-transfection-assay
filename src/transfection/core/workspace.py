"""Assay-specific analysis/results helpers (`analysis/PosN/chC.csv`, sample packs).

Folder names come from lisca when exported; path helpers for `roi/PosN` live in
`transfection.core.roi` and call `lisca.core.bbox` / `lisca.core.workspace`.
"""

from __future__ import annotations

import re
from pathlib import Path

from transfection.core.constants import ANALYSIS_DIRNAME, RESULTS_DIRNAME
from transfection.core.sample import SampleMapping

_TRACE_ALPHA = 0.1
_WORKSPACE_METRICS_STEM = re.compile(r"^ch\d+$")
_POS_DIR = re.compile(r"^Pos(\d+)$")
_CH_STEM = re.compile(r"^ch(\d+)$")


def trace_color_alpha_from_fluor_name(name: str) -> tuple[str, float]:
    haystack = name.lower()
    if "egfp" in haystack:
        color = "green"
    elif "mcherry" in haystack:
        color = "red"
    elif "gfp" in haystack:
        color = "green"
    elif "yfp" in haystack:
        color = "yellow"
    elif "bfp" in haystack:
        color = "blue"
    else:
        color = "gray"
    return (color, _TRACE_ALPHA)


def is_trace_csv(path: Path) -> bool:
    return bool(_WORKSPACE_METRICS_STEM.fullmatch(path.stem))


def workspace_analysis_dir(workspace: Path) -> Path:
    return workspace.resolve() / ANALYSIS_DIRNAME


def workspace_results_dir(workspace: Path) -> Path:
    return workspace.resolve() / RESULTS_DIRNAME


def analysis_position_dir(workspace: Path, position: int) -> Path:
    return workspace_analysis_dir(workspace) / f"Pos{position}"


def default_position_trace_csv_path(
    workspace: Path,
    position: int,
    signal_channel: int,
) -> Path:
    return (analysis_position_dir(workspace, position) / f"ch{signal_channel}.csv").resolve()


def analysis_position_table_csv(workspace: Path, position: int, kind: str) -> Path:
    return (analysis_position_dir(workspace, position) / f"{kind}.csv").resolve()


def discover_trace_csvs(analysis_dir: Path) -> list[Path]:
    if not analysis_dir.is_dir():
        raise ValueError(
            f"Expected {ANALYSIS_DIRNAME}/ directory at {analysis_dir}. "
            "Run transfection traces first."
        )
    csvs = sorted(
        analysis_dir.glob("Pos*/ch*.csv"),
        key=lambda path: (path.parent.name, path.name),
    )
    if not csvs:
        raise ValueError(f"No CSV metrics files in {analysis_dir}")
    metrics = [path for path in csvs if is_trace_csv(path)]
    if not metrics:
        raise ValueError(
            f"No position metrics CSV files (expected Pos{{position}}/ch{{channel}}.csv) in {analysis_dir}"
        )
    return metrics


def discover_analysis_table_csvs(workspace: Path, kind: str) -> list[Path]:
    analysis_dir = workspace_analysis_dir(workspace)
    if not analysis_dir.is_dir():
        raise ValueError(
            f"Expected {ANALYSIS_DIRNAME}/ directory at {analysis_dir}. "
            f"Run transfection {kind} first."
        )
    csvs = sorted(
        analysis_dir.glob(f"Pos*/{kind}.csv"),
        key=lambda path: path.parent.name,
    )
    if not csvs:
        raise ValueError(
            f"No {kind}.csv files in {analysis_dir}/PosN/. Run transfection {kind} first."
        )
    return csvs


def parse_trace_path(csv_path: Path) -> tuple[int, int]:
    """Return ``(position, signal_channel)`` from ``analysis/Pos{n}/ch{n}.csv``."""
    parent_match = _POS_DIR.fullmatch(csv_path.parent.name)
    stem_match = _CH_STEM.fullmatch(csv_path.stem)
    if parent_match is None or stem_match is None:
        raise ValueError(
            f"Expected analysis path Pos{{position}}/ch{{channel}}.csv, got {csv_path}"
        )
    return int(parent_match.group(1)), int(stem_match.group(1))


def parse_analysis_position_dir(path: Path) -> int:
    match = _POS_DIR.fullmatch(path.name)
    if match is None:
        raise ValueError(f"Expected analysis/Pos{{n}} directory, got {path}")
    return int(match.group(1))


def build_position_signal_sample_lookup(mapping: SampleMapping) -> dict[tuple[int, int], str]:
    lookup: dict[tuple[int, int], str] = {}
    for sample, entry in mapping.items():
        for position in entry.positions:
            for signal_channel in entry.signal_channels:
                key = (position, signal_channel)
                if key in lookup and lookup[key] != sample:
                    raise ValueError(
                        f"Ambiguous sample for position {position} "
                        f"signal channel {signal_channel}: "
                        f"{lookup[key]!r} and {sample!r}"
                    )
                lookup[key] = sample
    return lookup


def resolve_sample(csv_path: Path, mapping: SampleMapping) -> str:
    """Sample name owning ``analysis/Pos{n}/ch{c}.csv``."""
    position, signal_channel = parse_trace_path(csv_path)
    lookup = build_position_signal_sample_lookup(mapping)
    key = (position, signal_channel)
    if key not in lookup:
        raise ValueError(
            f"No assay mapping entry for Pos{position} signal channel {signal_channel} ({csv_path})"
        )
    return lookup[key]


def infer_workspace_for_plot_csv(csv_file: Path) -> Path:
    parent = csv_file.parent.resolve()
    if parent.name == RESULTS_DIRNAME:
        return parent.parent
    if parent.parent.name == RESULTS_DIRNAME:
        return parent.parent.parent
    if parent.parent.name == ANALYSIS_DIRNAME:
        return parent.parent.parent
    if parent.name == ANALYSIS_DIRNAME:
        return parent.parent
    return parent


def infer_workspace_for_analysis_dir(analysis_dir: Path) -> Path:
    resolved = analysis_dir.resolve()
    if resolved.name == ANALYSIS_DIRNAME:
        return resolved.parent
    return resolved.parent.resolve()


def infer_workspace_root(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if resolved.is_file():
        return infer_workspace_for_plot_csv(resolved)
    if (resolved / "assay.json").is_file() or (resolved / ANALYSIS_DIRNAME).is_dir():
        return resolved
    if resolved.name == ANALYSIS_DIRNAME:
        return resolved.parent
    if resolved.parent.name == RESULTS_DIRNAME:
        return resolved.parent.parent
    return resolved


def boxplot_tick_labels(samples: list[str], trace_counts: list[int]) -> list[str]:
    # Single-line labels so tilted x-ticks stay readable.
    return [
        f"{sample} (n={n})"
        for sample, n in zip(samples, trace_counts, strict=True)
    ]
