"""Workspace `assay.json` — single config for the transfection pipeline.

Schema matches LiSCA Studio / `@lisca/contracts` (`AssayJsonFile`). Agents write
this file directly; there is no generator CLI in this package. See `AGENTS.md`.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from transfection.core.slide import (
    SlideChannelMapping,
    SlideMapping,
    parse_position_spec,
    validate_slide_mapping,
)

ASSAY_FILENAME = "assay.json"
# Defaults when assay.json omits fields (transfection assay only — this package).
DEFAULT_INTERVAL_MINUTES = 10.0
# Second-pass onset-time (t0) search cap (minutes). Explicit 0 still means onset fixed at 0.
# Basic translation–degradation model only (no protein maturation).
DEFAULT_MAX_ONSET_MINUTES = 120.0


@dataclass(frozen=True)
class AssayConfig:
    path: Path
    assay_type: str
    name: str
    data_path: str
    mapping: SlideMapping
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
    return load_assay(resolve_assay_path(workspace, assay))


def build_slide_mapping_from_samples(samples: list[Any], *, source: Path | str) -> SlideMapping:
    if not isinstance(samples, list):
        raise ValueError(f"{source}: samples must be an array")

    mapping: SlideMapping = {}
    for index, row in enumerate(samples):
        if not isinstance(row, dict):
            raise ValueError(f"{source}: samples[{index}] must be an object")

        sample_name = str(row.get("name", "")).strip()
        if not sample_name:
            continue

        slide = _require_nonneg_int_field(row, "slide", source=source, index=index)
        signal_channel = _require_nonneg_int_field(row, "fluorescence", source=source, index=index)
        mask_channel = _require_nonneg_int_field(row, "brightfield", source=source, index=index)
        positions_raw = row.get("positions")
        if positions_raw is None or str(positions_raw).strip() == "":
            raise ValueError(f"{source}: samples[{index}] missing positions")
        try:
            positions = parse_position_spec(str(positions_raw))
        except ValueError as exc:
            raise ValueError(f"{source}: samples[{index}] positions: {exc}") from exc

        mapping[slide] = SlideChannelMapping(
            positions=positions,
            signal_channel=signal_channel,
            mask_channel=mask_channel,
            sample_name=sample_name,
        )

    return validate_slide_mapping(mapping)


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


def _parse_assay(raw: dict[str, Any], *, path: Path) -> AssayConfig:
    assay_type = str(raw.get("type") or "").strip() or "unknown"
    name = str(raw.get("name") or "").strip() or assay_type

    data = raw.get("data") if isinstance(raw.get("data"), dict) else {}
    data_path = str(data.get("path") or "").strip() if isinstance(data, dict) else ""

    samples = raw.get("samples")
    if samples is None:
        raise ValueError(f"{path}: missing samples array")
    mapping = build_slide_mapping_from_samples(samples, source=path)

    interval_obj = raw.get("interval") if isinstance(raw.get("interval"), dict) else {}
    interval = parse_interval_minutes(
        interval_obj.get("value") if isinstance(interval_obj, dict) else None,
        interval_obj.get("unit") if isinstance(interval_obj, dict) else None,
    )
    if interval is None:
        interval = DEFAULT_INTERVAL_MINUTES

    analysis = raw.get("analysis") if isinstance(raw.get("analysis"), dict) else {}
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
        interval_minutes=interval,
        max_onset_minutes=max_onset,
        skip_segment=skip_segment,
    )


def _require_nonneg_int_field(
    row: dict[str, Any],
    field: str,
    *,
    source: Path | str,
    index: int,
) -> int:
    if field not in row:
        raise ValueError(f"{source}: samples[{index}] missing {field}")
    raw = row[field]
    try:
        value = int(str(raw).strip())
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{source}: samples[{index}].{field} must be an integer, got {raw!r}") from exc
    if value < 0:
        raise ValueError(f"{source}: samples[{index}].{field} must be non-negative, got {value}")
    return value
