from transfection.core import (
    figure_size_for_grid,
    figure_size_for_panels,
    resolve_subplot_grid,
    subplot_grid_shape,
)


def test_subplot_grid_shape_one_to_twelve() -> None:
    expected = {
        1: (1, 1),
        2: (1, 2),
        3: (2, 2),
        4: (2, 2),
        5: (2, 3),
        6: (2, 3),
        7: (3, 3),
        8: (3, 3),
        9: (3, 3),
        10: (3, 4),
        11: (3, 4),
        12: (3, 4),
    }
    for count, shape in expected.items():
        assert subplot_grid_shape(count) == shape
        assert resolve_subplot_grid(count) == shape


def test_subplot_grid_shape_above_twelve_is_near_square() -> None:
    assert subplot_grid_shape(13) == (4, 4)
    assert subplot_grid_shape(16) == (4, 4)


def test_resolve_subplot_grid_respects_explicit_columns() -> None:
    assert resolve_subplot_grid(4, columns=3) == (2, 3)
    assert resolve_subplot_grid(12, columns=6) == (2, 6)


def test_figure_size_scales_with_grid() -> None:
    assert figure_size_for_panels(1) == figure_size_for_grid(1, 1)
    size_4 = figure_size_for_panels(4)
    size_12 = figure_size_for_panels(12)
    assert size_4 == (2 * 5.5, 2 * 4.0)
    # 3×4 uses denser scale so dual-slide figures stay manageable.
    assert size_12[0] > size_4[0]
    assert size_12[1] > size_4[1]
