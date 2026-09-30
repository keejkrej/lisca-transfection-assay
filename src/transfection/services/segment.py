from __future__ import annotations

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from transfection.core import (
    SampleAnalysis,
    SampleMapping,
    compute_roi_mask_stack,
    default_mask_path,
    discover_roi_positions,
    load_assay_for_workspace,
    position_dir,
    read_position_index,
    validate_channel_index,
    write_mask_tif,
)
from transfection.core.parallel import worker_count



MaskWrittenCallback = Callable[[str, Path, int], None]


@dataclass(frozen=True)
class SegmentationRunResult:
    # (sample name, mask directory, mask count) in assay order.
    written_outputs: list[tuple[str, Path, int]]
    skipped_positions: dict[str, list[int]]


def _run_position_segmentation(
    workspace: Path,
    *,
    sample: str,
    segmentation_channel: int,
    resolved_pos: int,
    variation_radius: int,
    gaussian_sigma: float,
    force: bool,
) -> tuple[str, int, int, int, Path | None]:
    try:
        pos_dir = position_dir(workspace, resolved_pos)
    except ValueError:
        return (sample, segmentation_channel, resolved_pos, 0, None)

    index = read_position_index(pos_dir)
    validate_channel_index(index, segmentation_channel)
    mask_count = 0
    first_output: Path | None = None

    for roi in index.rois:
        output_path = default_mask_path(
            workspace,
            position=index.position,
            roi_file_name=roi.file_name,
        )
        if output_path.exists() and not force:
            mask_count += 1
            if first_output is None:
                first_output = output_path
            continue

        mask_stack = compute_roi_mask_stack(
            pos_dir,
            index,
            roi,
            channel=segmentation_channel,
            variation_radius=variation_radius,
            gaussian_sigma=gaussian_sigma,
        )
        write_mask_tif(mask_stack, output_path)
        mask_count += 1
        if first_output is None:
            first_output = output_path

    return (sample, segmentation_channel, resolved_pos, mask_count, first_output)


def _position_segmentation_task(
    payload: tuple[str, str, int, int, int, float, bool],
) -> tuple[str, int, int, int, Path | None]:
    workspace_str, sample, segmentation_channel, resolved_pos, variation_radius, gaussian_sigma, force = payload
    return _run_position_segmentation(
        Path(workspace_str),
        sample=sample,
        segmentation_channel=segmentation_channel,
        resolved_pos=resolved_pos,
        variation_radius=variation_radius,
        gaussian_sigma=gaussian_sigma,
        force=force,
    )


def _position_tasks(
    workspace: Path,
    samples: SampleMapping,
    *,
    variation_radius: int,
    gaussian_sigma: float,
    force: bool,
) -> list[tuple[str, str, int, int, int, float, bool]]:
    return [
        (
            str(workspace),
            sample,
            entry.segmentation_channel,
            resolved_pos,
            variation_radius,
            gaussian_sigma,
            force,
        )
        for sample, entry in samples.items()
        for resolved_pos in entry.positions
    ]


def skipped_positions_summary(skipped_positions: dict[str, list[int]]) -> str:
    return "; ".join(
        f"sample {sample!r} -> {', '.join(str(pos) for pos in positions)}"
        for sample, positions in skipped_positions.items()
    )


def run_segmentation_for_mapping(
    workspace: Path,
    *,
    mapping: SampleMapping | None,
    segmentation_channel: int | None = None,
    variation_radius: int = 2,
    gaussian_sigma: float = 1.0,
    force: bool = False,
    on_mask_written: MaskWrittenCallback | None = None,
) -> SegmentationRunResult:
    """Segment every Position of each Sample. Without samples[], every roi/PosN."""
    if variation_radius < 0:
        raise ValueError(f"--variation-radius must be >= 0, got {variation_radius}")
    if gaussian_sigma < 0:
        raise ValueError(f"--gaussian-sigma must be >= 0, got {gaussian_sigma}")

    workspace = workspace.resolve()
    if mapping:
        samples = mapping
    else:
        if segmentation_channel is None:
            segmentation_channel = load_assay_for_workspace(workspace).segmentation_channel
        positions = discover_roi_positions(workspace)
        # Unnamed placeholder: no samples[] means one group of every roi/PosN.
        samples = {
            "": SampleAnalysis(
                name="",
                positions=positions,
                signal_channels=[0],
                segmentation_channel=segmentation_channel,
            )
        }
    tasks = _position_tasks(
        workspace,
        samples,
        variation_radius=variation_radius,
        gaussian_sigma=gaussian_sigma,
        force=force,
    )
    if not tasks:
        raise ValueError("assay mapping defines no valid positions")

    skipped_positions: dict[str, list[int]] = defaultdict(list)
    written_by_sample: dict[str, tuple[Path, int]] = {}

    def consume(row: tuple[str, int, int, int, Path | None]) -> None:
        sample, _segmentation_channel, resolved_pos, written_count, first_output = row
        if first_output is None:
            skipped_positions[sample].append(resolved_pos)
            return
        current_output, current_count = written_by_sample.get(sample, (first_output, 0))
        written_by_sample[sample] = (current_output, current_count + written_count)

    max_workers = worker_count(len(tasks))
    if max_workers == 1:
        for task in tasks:
            consume(_position_segmentation_task(task))
    else:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(_position_segmentation_task, task) for task in tasks]
            for future in as_completed(futures):
                consume(future.result())

    written_outputs: list[tuple[str, Path, int]] = []
    for sample in samples:
        if sample not in written_by_sample:
            continue
        first_output, mask_count = written_by_sample[sample]
        written_outputs.append((sample, first_output.parent, mask_count))
        if on_mask_written is not None:
            on_mask_written(sample, first_output.parent, mask_count)

    ordered_skipped = {
        sample: sorted(set(skipped_positions[sample]))
        for sample in samples
        if skipped_positions.get(sample)
    }
    if not written_outputs:
        if ordered_skipped:
            raise ValueError(
                f"No ROI directories found for positions in assay mapping. "
                f"Skipped positions: {skipped_positions_summary(ordered_skipped)}"
            )
        raise ValueError("assay mapping defines no valid positions")

    return SegmentationRunResult(
        written_outputs=written_outputs,
        skipped_positions=ordered_skipped,
    )


def format_written_masks_message(sample: str, output_dir: Path, mask_count: int) -> str:
    noun = "mask" if mask_count == 1 else "masks"
    target = f"sample {sample!r}" if sample else "roi/"
    return f"Prepared {mask_count} {noun} for {target} under: {output_dir}"


def format_skipped_positions_message(skipped_positions: dict[str, list[int]]) -> str:
    total_skipped_positions = sum(len(positions) for positions in skipped_positions.values())
    return (
        f"Skipped {total_skipped_positions} missing positions from sample mapping: "
        f"{skipped_positions_summary(skipped_positions)}"
    )


def run_segment(
    *,
    workspace: Path,
    assay: Path | None = None,
    mapping: SampleMapping | None = None,
    variation_radius: int = 2,
    gaussian_sigma: float = 1.0,
    force: bool = False,
    on_mask_written: MaskWrittenCallback | None = None,
) -> SegmentationRunResult:
    config = load_assay_for_workspace(workspace, assay)
    if mapping is None:
        mapping = config.mapping or None
    return run_segmentation_for_mapping(
        workspace,
        mapping=mapping,
        segmentation_channel=config.segmentation_channel,
        variation_radius=variation_radius,
        gaussian_sigma=gaussian_sigma,
        force=force,
        on_mask_written=on_mask_written,
    )
