use std::path::Path;

use mplot::prelude::{Figure, FigureBuilder, LineStyle, SaveOptions, Size};
use mplot::Color;

/// Single-panel charts (inches).
/// Height is a bit taller than the transfection 4.5 in default so tilted
/// category x-labels keep a usable axes box after the bottom gutter.
pub const FIGURE_SINGLE_WIDTH_IN: f64 = 6.5;
pub const FIGURE_SINGLE_HEIGHT_IN: f64 = 5.5;

/// Matches transfection `FIGURE_DPI`.
pub const FIGURE_DPI: u32 = 100;
/// Matches transfection `font.size` / `axes.titlesize` / `axes.labelsize`.
pub const LABEL_FONTSIZE: f64 = 18.0;
pub const TITLE_FONTSIZE: f64 = 18.0;
/// Matches transfection `xtick.labelsize` / `ytick.labelsize`.
pub const TICK_FONTSIZE: f64 = 17.0;

/// One axes (traces, summary, area, scatter, AUC / fit parameter boxplots).
pub fn figure_builder_single() -> FigureBuilder {
    Figure::builder()
        .size(Size::inches(
            FIGURE_SINGLE_WIDTH_IN,
            FIGURE_SINGLE_HEIGHT_IN,
        ))
        .label_fontsize(LABEL_FONTSIZE)
        .tick_fontsize(TICK_FONTSIZE)
        .title_fontsize(TITLE_FONTSIZE)
        .gaps(0.2, 0.2)
}

/// Square-ish joint plot (scatter + attached x/y histograms). Small gaps so the
/// three axes read as one framed panel.
pub const FIGURE_JOINT_WIDTH_IN: f64 = 6.5;
pub const FIGURE_JOINT_HEIGHT_IN: f64 = 6.5;

pub fn figure_builder_joint() -> FigureBuilder {
    Figure::builder()
        .size(Size::inches(FIGURE_JOINT_WIDTH_IN, FIGURE_JOINT_HEIGHT_IN))
        .label_fontsize(LABEL_FONTSIZE)
        .tick_fontsize(TICK_FONTSIZE)
        .title_fontsize(TITLE_FONTSIZE)
        .gaps(0.04, 0.04)
        .constrained_layout(true)
}

pub fn trace_line_style(color_name: &str, alpha: f64) -> LineStyle {
    LineStyle::new()
        .color(Color::hex(color_name))
        .alpha(alpha)
        .width(1.5) // matplotlib lines.linewidth default / transfection traces
}

/// Outer figure pad for single-panel charts (`SaveOptions::pad_inches`).
/// Larger than matplotlib's 0.1 default so tilted category labels and spines
/// aren't flush with the PNG edge.
pub const SAVE_PAD_SINGLE_INCHES: f64 = 0.25;

/// Save a figure, passing `pad_inches` through to mplot (`SaveOptions::pad_inches`).
pub fn save_figure(figure: &Figure, path: &Path, pad_inches: f64) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|error| error.to_string())?;
    }
    figure
        .save(
            path,
            SaveOptions::new().dpi(FIGURE_DPI).pad_inches(pad_inches),
        )
        .map_err(|error| error.to_string())
}
