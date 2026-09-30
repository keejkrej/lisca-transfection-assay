"""Workspace `assay.json` — single config for the transfection pipeline.

Schema matches LiSCA Studio / `@lisca/contracts` (`AssayJsonFile`). Agents write
this file directly; there is no generator CLI in this package. See `AGENTS.md`.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from lisca.migrations import migrate_workspace

from transfection.core.sample import (
    SampleMapping,
    build_sample_mapping,
    parse_default_channels,
    parse_sample_channel_overrides,
)

ASSAY_FILENAME = "assay.json"
# Defaults when assay.json omits fields (transfection assay only — this package).
DEFAULT_INTERVAL_MINUTES = 10.0
# Second-pass onset-time (t0) search cap (minutes). Explicit 0 still means onset fixed at 0.
# Basic translation–degradation model only (no protein maturation).
DEFAULT_MAX_ONSET_MINUTES = 120.0


MISSING_SAMPLES_FOR_PLOT = (
    "plot/results stages require assay.json samples[] to group analysis/ into "
    "results/<sample>/. traces, auc, and fit do not need samples."
)


@dataclass(frozen=True)
class AssayConfig:
    path: Path
    assay_type: str
    name: str
    data_path: str
    mapping: SampleMapping
    segmentation_channel: int
    signal_channels: tuple[int, ...]
    interval_minutes: float | None
    max_onset_minutes: float
    skip_segment: bool


def resolve_assay_path(workspace: Path, assay: Path | None = None) -> Path:
    if assay is None:
        return (workspace / ASSAY_FILENAME).resolve()
    return assay.expanduser().resolve()


def load_assay(path: Path | str) -> AssayConfig:
    path = Path(path).expanduser().resolve()
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ValueError(f"missing assay.json at {path}") from exc
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid assay.json {path}: {exc}") from exc
    if not isinstance(raw, dict):
        raise ValueError(f"assay.json must be a JSON object: {path}")
    return _parse_assay(raw, path=path)


def load_assay_for_workspace(workspace: Path, assay: Path | None = None) -> AssayConfig:
    """Load the workspace assay.json after migrating the workspace on disk.

    An explicit ``assay`` path outside the workspace is read as-is.
    """
    migrate_workspace(Path(workspace))
    return load_assay(resolve_assay_path(workspace, assay))


def parse_interval_minutes(amount: object, unit: object) -> float | None:
    if amount is None:
        return None
    try:
        value = float(amount)
    except (TypeError, ValueError):
        return None
    if value <= 0:
        return None
    unit_str = "minute" if unit is None else str(unit)
    factor = {
        "second": 1.0 / 60.0,
        "minute": 1.0,
        "hour": 60.0,
    }.get(unit_str)
    if factor is None:
        return None
    return value * factor


def require_interval_minutes(config: AssayConfig, *, override: float | None = None) -> float:
    if override is not None:
        if override <= 0:
            raise ValueError(f"--interval must be > 0, got {override}")
        return override
    if config.interval_minutes is None:
        raise ValueError(
            f"missing --interval and could not read a positive interval.value from {config.path}"
        )
    return config.interval_minutes


def resolve_interval_minutes(
    workspace: Path,
    *,
    assay: Path | None = None,
    override: float | None = None,
) -> float:
    """Interval from ``--interval`` or assay.json; errors if neither is available."""
    if override is not None:
        if override <= 0:
            raise ValueError(f"--interval must be > 0, got {override}")
        return override
    assay_path = resolve_assay_path(workspace, assay)
    if not assay_path.is_file():
        raise ValueError(
            f"missing --interval and no assay.json at {assay_path} "
            "(pass --interval when replotting a copied sample folder)"
        )
    return require_interval_minutes(load_assay(assay_path))


def require_samples(config: AssayConfig) -> SampleMapping:
    """Plot/results grouping. Fails when assay.json defines no ``samples[]``."""
    if not config.mapping:
        raise ValueError(f"{config.path}: {MISSING_SAMPLES_FOR_PLOT}")
    return config.mapping


def _parse_assay(raw: dict[str, Any], *, path: Path) -> AssayConfig:
    assay_type = str(raw.get("type") or "").strip() or "unknown"
    name = str(raw.get("name") or "").strip() or assay_type

    data = raw.get("data") if isinstance(raw.get("data"), dict) else {}
    data_path = str(data.get("path") or "").strip() if isinstance(data, dict) else ""

    samples = raw.get("samples")
    analysis = raw.get("analysis") if isinstance(raw.get("analysis"), dict) else {}
    if samples is None:
        mapping = {}
    elif not isinstance(samples, list):
        raise ValueError(f"{path}: samples must be an array")
    elif not samples:
        mapping = {}
    else:
        mapping = build_sample_mapping(
            samples, analysis if isinstance(analysis, dict) else None, source=path
        )

    default_segmentation, default_signal = parse_default_channels(
        analysis if isinstance(analysis, dict) else None, source=path
    )
    if default_segmentation is None or default_signal is None:
        raise ValueError(f"{path}: missing analysis.channels")
    signal_channels = tuple(default_signal)
    extra_signals = parse_sample_channel_overrides(
        analysis if isinstance(analysis, dict) else None, source=path
    )
    if extra_signals:
        merged = list(signal_channels)
        for _segmentation, signals in extra_signals.values():
            for channel in signals:
                if channel not in merged:
                    merged.append(channel)
        signal_channels = tuple(merged)

    interval_obj = raw.get("interval") if isinstance(raw.get("interval"), dict) else {}
    interval = parse_interval_minutes(
        interval_obj.get("value") if isinstance(interval_obj, dict) else None,
        interval_obj.get("unit") if isinstance(interval_obj, dict) else None,
    )
    if interval is None:
        interval = DEFAULT_INTERVAL_MINUTES

    max_onset = DEFAULT_MAX_ONSET_MINUTES
    if isinstance(analysis, dict) and analysis.get("maxOnsetMinutes") is not None:
        try:
            max_onset = float(analysis["maxOnsetMinutes"])
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{path}: analysis.maxOnsetMinutes must be a number") from exc
        if max_onset < 0:
            raise ValueError(f"{path}: analysis.maxOnsetMinutes must be >= 0")

    skip_segment = False
    if isinstance(analysis, dict) and analysis.get("skipSegment") is not None:
        skip_raw = analysis["skipSegment"]
        if not isinstance(skip_raw, bool):
            raise ValueError(f"{path}: analysis.skipSegment must be a boolean")
        skip_segment = skip_raw

    return AssayConfig(
        path=path,
        assay_type=assay_type,
        name=name,
        data_path=data_path,
        mapping=mapping,
        segmentation_channel=default_segmentation,
        signal_channels=signal_channels,
        interval_minutes=interval,
        max_onset_minutes=max_onset,
        skip_segment=skip_segment,
    )
