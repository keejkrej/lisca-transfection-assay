"""Acquisition frame labels become minutes by multiplying the interval."""

from __future__ import annotations

import numpy as np

from transfection.core.frame_time import frame_minutes


def test_dense_labels_match_label_times_interval() -> None:
    labels = np.array([0.0, 1.0, 2.0, 3.0])
    np.testing.assert_allclose(frame_minutes(labels, 10.0), labels * 10.0)


def test_stepped_acquisition_labels_scale_with_interval() -> None:
    labels = np.array([0.0, 6.0, 12.0, 18.0])
    np.testing.assert_allclose(frame_minutes(labels, 10.0), [0.0, 60.0, 120.0, 180.0])


def test_unsorted_labels_keep_their_minutes() -> None:
    labels = np.array([12.0, 0.0, 6.0])
    np.testing.assert_allclose(frame_minutes(labels, 10.0), [120.0, 0.0, 60.0])
