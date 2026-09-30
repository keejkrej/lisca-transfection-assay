"""In-memory Sample mapping used by analysis stages.

Built from Studio `assay.json` (`samples[]` + `analysis.channels` /
`analysis.sampleChannels`). A Sample is identified by its name, unique within
the assay. Position ranges use inclusive Studio semantics (`1:12` → 1…12),
matching `crates/lisca`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class SampleAnalysis:
    name: str
    positions: list[int]
    signal_channels: list[int]
    segmentation_channel: int


# Keyed by sample name, insertion-ordered in assay order (never sorted).
type SampleMapping = dict[str, SampleAnalysis]


def parse_position_token(token: str) -> list[int]:
    """Expand one position token. Ranges are inclusive on both ends (`1:3` → 1,2,3)."""
    raw = token.strip()
    if not raw:
        raise ValueError("Empty position token")

    if ":" not in raw:
        try:
            value = int(raw)
        except ValueError as exc:
            raise ValueError(f"Invalid position token: {raw!r}") from exc
        if value < 0:
            raise ValueError(f"Positions must be non-negative, got {value}")
        return [value]

    parts = [part.strip() for part in raw.split(":")]
    if len(parts) not in {2, 3}:
        raise ValueError(f"Invalid slice token: {raw!r}")
    if any(part == "" for part in parts[:2]):
        raise ValueError(f"Slices must include explicit start and stop: {raw!r}")

    try:
        start = int(parts[0])
        stop = int(parts[1])
        step = int(parts[2]) if len(parts) == 3 else 1
    except ValueError as exc:
        raise ValueError(f"Invalid slice token: {raw!r}") from exc

    if start < 0 or stop < 0:
        raise ValueError(f"Positions must be non-negative in slice {raw!r}")
    if step <= 0:
        raise ValueError(f"Slice step must be > 0 in {raw!r}")
    if stop < start:
        raise ValueError(f"Invalid empty position range: {raw!r}")

    values = list(range(start, stop + 1, step))
    if not values:
        raise ValueError(f"Slice produced no positions: {raw!r}")
    return values


def parse_position_spec(spec: str) -> list[int]:
    tokens = [token.strip() for token in spec.split(",")]
    if not any(tokens):
        raise ValueError("Position spec is empty")

    positions: list[int] = []
    for token in tokens:
        if not token:
            raise ValueError("Position spec contains an empty token")
        positions.extend(parse_position_token(token))

    return sorted(set(positions))


def validate_sample_mapping(mapping: SampleMapping) -> SampleMapping:
    """Check each Sample and normalize positions. Preserves assay order."""
    if not mapping:
        raise ValueError("sample mapping defines no samples")
    ordered: SampleMapping = {}
    for key, entry in mapping.items():
        name = entry.name.strip()
        if not name:
            raise ValueError("sample name must be non-empty")
        if name != key:
            raise ValueError(f"sample mapping key {key!r} does not match name {name!r}")
        if not entry.positions:
            raise ValueError(f"sample {name!r} defines no positions")
        if not entry.signal_channels:
            raise ValueError(f"sample {name!r}: signal channel list must be non-empty")
        for signal_channel in entry.signal_channels:
            if signal_channel < 0:
                raise ValueError(f"signal channel must be non-negative, got {signal_channel}")
        if entry.segmentation_channel < 0:
            raise ValueError(
                f"segmentation channel must be non-negative, got {entry.segmentation_channel}"
            )
        ordered[name] = SampleAnalysis(
            name=name,
            positions=sorted(set(entry.positions)),
            signal_channels=list(entry.signal_channels),
            segmentation_channel=entry.segmentation_channel,
        )
    return ordered


def build_sample_mapping(
    samples: list[Any],
    analysis: dict[str, Any] | None,
    *,
    source: Path | str,
) -> SampleMapping:
    """Build the mapping from ``samples[]`` and ``analysis`` in assay order.

    Each Sample takes its ``analysis.sampleChannels`` row (matched by name) or
    falls back to ``analysis.channels``.
    """
    if not isinstance(samples, list):
        raise ValueError(f"{source}: samples must be an array")

    default_segmentation, default_signal = parse_default_channels(analysis, source=source)
    overrides = parse_sample_channel_overrides(analysis, source=source)

    mapping: SampleMapping = {}
    for index, row in enumerate(samples):
        if not isinstance(row, dict):
            raise ValueError(f"{source}: samples[{index}] must be an object")

        raw_name = row.get("name")
        name = raw_name.strip() if isinstance(raw_name, str) else ""
        if not name:
            raise ValueError(f"{source}: samples[{index}]: sample name must be non-empty")
        if name in mapping:
            raise ValueError(f"{source}: duplicate sample name {name!r} in samples[]")

        if name in overrides:
            segmentation_channel, signal_channels = overrides[name]
        elif default_segmentation is not None and default_signal is not None:
            segmentation_channel, signal_channels = default_segmentation, default_signal
        else:
            raise ValueError(
                f"{source}: missing analysis.channels (and no sampleChannels override) "
                f"for sample {name!r}"
            )

        positions_raw = row.get("positions")
        if positions_raw is None or str(positions_raw).strip() == "":
            raise ValueError(f"{source}: samples[{index}] missing positions")
        try:
            positions = parse_position_spec(str(positions_raw))
        except ValueError as exc:
            raise ValueError(f"{source}: samples[{index}] positions: {exc}") from exc

        mapping[name] = SampleAnalysis(
            name=name,
            positions=positions,
            signal_channels=list(signal_channels),
            segmentation_channel=segmentation_channel,
        )

    unknown = [name for name in overrides if name not in mapping]
    if unknown:
        raise ValueError(f"{source}: analysis.sampleChannels: unknown sample {unknown[0]!r}")

    return validate_sample_mapping(mapping)


def parse_default_channels(
    analysis: dict[str, Any] | None,
    *,
    source: Path | str,
) -> tuple[int | None, list[int] | None]:
    """``analysis.channels`` as ``(segmentation, signal)``; ``(None, None)`` if absent."""
    if not isinstance(analysis, dict):
        return None, None
    channels = analysis.get("channels")
    if channels is None:
        return None, None
    if not isinstance(channels, dict):
        raise ValueError(f"{source}: analysis.channels must be an object")
    segmentation = _require_nonneg_int_value(
        channels.get("segmentation"), field="analysis.channels.segmentation", source=source
    )
    signal = _require_signal_list(
        channels.get("signal"), field="analysis.channels.signal", source=source
    )
    return segmentation, signal


def parse_sample_channel_overrides(
    analysis: dict[str, Any] | None,
    *,
    source: Path | str,
) -> dict[str, tuple[int, list[int]]]:
    """``analysis.sampleChannels`` keyed by sample name, in file order."""
    if not isinstance(analysis, dict):
        return {}
    rows = analysis.get("sampleChannels")
    if rows is None:
        return {}
    if not isinstance(rows, list):
        raise ValueError(f"{source}: analysis.sampleChannels must be an array")

    overrides: dict[str, tuple[int, list[int]]] = {}
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            raise ValueError(f"{source}: analysis.sampleChannels[{index}] must be an object")
        raw_sample = row.get("sample")
        if not isinstance(raw_sample, str):
            raise ValueError(f"{source}: analysis.sampleChannels[{index}] missing sample")
        sample = raw_sample.strip()
        if sample in overrides:
            raise ValueError(f"{source}: analysis.sampleChannels: duplicate sample {sample!r}")
        segmentation = _require_nonneg_int_field(
            row, "segmentation", source=source, index=index, where="analysis.sampleChannels"
        )
        signal = _require_signal_list(
            row.get("signal"),
            field=f"analysis.sampleChannels[{index}].signal",
            source=source,
        )
        overrides[sample] = (segmentation, signal)
    return overrides


def _require_signal_list(raw: object, *, field: str, source: Path | str) -> list[int]:
    if not isinstance(raw, list) or len(raw) == 0:
        raise ValueError(f"{source}: {field} must be a non-empty array of integers")
    values: list[int] = []
    for item in raw:
        try:
            value = int(item)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{source}: {field} must contain integers, got {item!r}") from exc
        if value < 0:
            raise ValueError(f"{source}: {field} values must be non-negative, got {value}")
        values.append(value)
    return values


def _require_nonneg_int_value(raw: object, *, field: str, source: Path | str) -> int:
    if raw is None:
        raise ValueError(f"{source}: missing {field}")
    try:
        value = int(raw) if not isinstance(raw, str) else int(raw.strip())
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{source}: {field} must be an integer, got {raw!r}") from exc
    if value < 0:
        raise ValueError(f"{source}: {field} must be non-negative, got {value}")
    return value


def _require_nonneg_int_field(
    row: dict[str, Any],
    field: str,
    *,
    source: Path | str,
    index: int,
    where: str,
) -> int:
    if field not in row:
        raise ValueError(f"{source}: {where}[{index}] missing {field}")
    return _require_nonneg_int_value(
        row[field],
        field=f"{where}[{index}].{field}",
        source=source,
    )
