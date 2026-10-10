//! NumPy-style array helpers for ROI frames, masks, and shared numeric kernels.
//!
//! Morphology filters live in `image_ops` (`ndarray-ndimage`, reflect borders,
//! the same kernel the Python segmenter uses). Quantile interpolation is
//! `marray::np::percentile` and trapezoidal integrals are
//! `marray::sp::integrate::trapz`. Non-finite samples are dropped here first.
//! Do not switch the quantile path to `ndarray-stats` quickselect: a tied plate
//! of traces overflows the worker stack.

use ndarray::{Array1, ArrayView2};

#[derive(Debug, Clone)]
pub struct Frame2D {
    pub width: usize,
    pub height: usize,
    data: Vec<f64>,
}

impl Frame2D {
    pub fn from_vec(data: Vec<f64>, width: usize, height: usize) -> Result<Self, String> {
        if width == 0 || height == 0 {
            return Err("frame dimensions must be positive".to_string());
        }
        if data.len() != width * height {
            return Err(format!(
                "frame length {} does not match {width}x{height}",
                data.len()
            ));
        }
        Ok(Self {
            width,
            height,
            data,
        })
    }

    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }

    pub fn as_view(&self) -> ArrayView2<'_, f64> {
        ArrayView2::from_shape((self.height, self.width), &self.data)
            .expect("frame shape matches backing storage")
    }

    pub fn into_vec(self) -> Vec<f64> {
        self.data
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MaskedRoiStats {
    pub area: u32,
    pub intensity: f64,
    pub background: f64,
    pub corrected: f64,
}

/// Full-frame ROI reduction: whole crop as foreground; bg = 10th percentile.
pub const FULL_FRAME_BACKGROUND_QUANTILE: f64 = 0.1;

/// Masked ROI reduction matching transfection `compute_masked_roi_metrics`
/// (background = median of non-mask pixels).
pub fn masked_roi_stats(frame: &[f64], mask: &[bool]) -> Result<MaskedRoiStats, String> {
    if frame.len() != mask.len() {
        return Err(format!(
            "frame/mask length mismatch: {} vs {}",
            frame.len(),
            mask.len()
        ));
    }
    if frame.is_empty() {
        return Ok(MaskedRoiStats {
            area: 0,
            intensity: 0.0,
            background: 0.0,
            corrected: 0.0,
        });
    }

    let mut background_pixels = Vec::new();
    let mut intensity = 0.0;
    let mut area = 0u32;
    for (value, &is_fg) in frame.iter().zip(mask.iter()) {
        if is_fg {
            intensity += *value;
            area += 1;
        } else {
            background_pixels.push(*value);
        }
    }
    let background = if background_pixels.is_empty() {
        0.0
    } else {
        quantile(&background_pixels, 0.5)
    };
    let corrected = intensity - f64::from(area) * background;

    Ok(MaskedRoiStats {
        area,
        intensity,
        background,
        corrected,
    })
}

/// Full-frame ROI reduction matching transfection `analysis.skipSegment`.
pub fn full_frame_roi_stats(frame: &[f64]) -> MaskedRoiStats {
    if frame.is_empty() {
        return MaskedRoiStats {
            area: 0,
            intensity: 0.0,
            background: 0.0,
            corrected: 0.0,
        };
    }
    let area = frame.len() as u32;
    let intensity: f64 = frame.iter().sum();
    let background = quantile(frame, FULL_FRAME_BACKGROUND_QUANTILE);
    MaskedRoiStats {
        area,
        intensity,
        background,
        corrected: intensity - f64::from(area) * background,
    }
}

/// Linear interpolation quantile on unsorted `f64` values (`numpy.quantile` default).
pub fn quantile(values: &[f64], q: f64) -> f64 {
    quantile_linear(values, q)
}

/// Linear interpolation quantile on unsorted `f64` values (`numpy.quantile` linear).
///
/// Non-finite values are dropped. A single input is returned as-is, including
/// non-finite values. `q` outside `[0, 1]` is clamped; a non-finite `q` yields `0.0`.
pub fn quantile_linear(values: &[f64], q: f64) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    if values.len() == 1 {
        return values[0];
    }
    let Some(q) = finite_quantile_fraction(q) else {
        return 0.0;
    };
    let mut finite: Vec<f64> = values
        .iter()
        .copied()
        .filter(|value| value.is_finite())
        .collect();
    if finite.is_empty() {
        return 0.0;
    }
    if finite.len() == 1 {
        return finite[0];
    }
    finite.sort_unstable_by(f64::total_cmp);
    quantile_linear_sorted(&finite, q)
}

/// Linear interpolation quantile on a pre-sorted finite slice (`q` in `[0, 1]`).
///
/// `marray::np::percentile` uses index `q * (n - 1)`, matching
/// `numpy.quantile(..., method="linear")`.
pub fn quantile_linear_sorted(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    if sorted.len() == 1 {
        return sorted[0];
    }
    let Some(q) = finite_quantile_fraction(q) else {
        return 0.0;
    };
    marray::np::percentile(&marray::np::array(sorted.to_vec()), q * 100.0)
}

fn finite_quantile_fraction(q: f64) -> Option<f64> {
    q.is_finite().then_some(q.clamp(0.0, 1.0))
}

/// Percentile on unsorted `f64` values (`pct` in `[0, 100]`, linear interpolation).
pub fn percentile(values: &[f64], pct: f64) -> f64 {
    quantile(values, pct / 100.0)
}

/// Floor-index quantile on pre-sorted data (viewer contrast semantics).
pub fn quantile_floor_sorted(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let clamped = q.clamp(0.0, 1.0);
    let index = (clamped * (sorted.len().saturating_sub(1)) as f64).floor() as usize;
    sorted[index.min(sorted.len() - 1)]
}

/// Trapezoidal integration of `(times, values)` pairs (`numpy.trapz` parity).
pub fn trapezoidal_integral(times: &[f64], values: &[f64]) -> f64 {
    if times.len() < 2 || times.len() != values.len() {
        return 0.0;
    }
    let times = marray::np::array(times.to_vec());
    let values = marray::np::array(values.to_vec());
    marray::sp::integrate::trapz(&values, Some(&times), None)
}

/// Coefficients for the basic translation–degradation model (Müller et al. 2024
/// Eq. 3; **no** protein maturation). Used by the optimizer and by `traces_fit`
/// reconstruction. Written CSV/XLSX store paper observables only
/// (`onset_time`, `expression_rate`, lifetimes, `baseline_intensity`);
/// recover `protein_degradation_rate` β = ln(2)/protein_lifetime,
/// `mrna_degradation_rate` δ = ln(2)/mrna_lifetime, and
/// `expression_amplitude` = expression_rate / (δ − β) at plot time.
///
/// `I(t) = baseline_intensity + expression_amplitude * (e^{-β Δt} − e^{-δ Δt})`
/// for `t ≥ onset_time` (`t0`), else `baseline_intensity`.
///
/// Paper terms: `onset_time` = onset time \(t_0\);
/// `expression_rate = expression_amplitude * (δ − β)` = \(m_0 k_{TL}\);
/// `ln(2)/δ` = mRNA lifetime (half-life); `ln(2)/β` = protein lifetime (half-life).
/// `baseline_intensity` is a baseline nuisance, not a kinetic rate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KineticFitCoeffs {
    pub baseline_intensity: f64,
    pub protein_degradation_rate: f64,
    pub mrna_degradation_rate: f64,
    pub onset_time: f64,
    pub expression_amplitude: f64,
}

/// Paper lifetime: τ = ln(2)/rate (half-life). Stored in minutes when `rate` is per minute.
pub fn half_life_minutes(degradation_rate_per_minute: f64) -> f64 {
    std::f64::consts::LN_2 / degradation_rate_per_minute
}

/// Recover β or δ (per minute) from a stored half-life in minutes.
pub fn degradation_rate_per_minute(half_life_minutes_value: f64) -> f64 {
    std::f64::consts::LN_2 / half_life_minutes_value
}

/// Recover m0·kTL/(δ−β) from stored paper observables.
pub fn expression_amplitude_from_observables(
    expression_rate: f64,
    mrna_lifetime: f64,
    protein_lifetime: f64,
) -> f64 {
    let mrna_degradation_rate = degradation_rate_per_minute(mrna_lifetime);
    let protein_degradation_rate = degradation_rate_per_minute(protein_lifetime);
    expression_rate / (mrna_degradation_rate - protein_degradation_rate)
}

/// Paper Eq. (4): AUC = (ln 2)^2 · m0 k_TL · τ_mRNA · τ_EGFP.
/// Lifetimes must be half-lives. This is not the trace-integrated fluorescence AUC.
pub fn auc_from_fit_half_lives(
    expression_rate: f64,
    mrna_lifetime: f64,
    protein_lifetime: f64,
) -> f64 {
    std::f64::consts::LN_2.powi(2) * expression_rate * mrna_lifetime * protein_lifetime
}

pub fn kinetic_basis_value(
    time: f64,
    protein_degradation_rate: f64,
    mrna_degradation_rate: f64,
    onset_time: f64,
) -> f64 {
    if time < onset_time {
        return 0.0;
    }
    let dt = time - onset_time;
    (-protein_degradation_rate * dt).exp() - (-mrna_degradation_rate * dt).exp()
}

pub fn fitted_trace_value(time: f64, coeffs: &KineticFitCoeffs) -> f64 {
    coeffs.baseline_intensity
        + coeffs.expression_amplitude
            * kinetic_basis_value(
                time,
                coeffs.protein_degradation_rate,
                coeffs.mrna_degradation_rate,
                coeffs.onset_time,
            )
}

pub fn evaluate_kinetic_candidate(
    times: &[f64],
    values: &[f64],
    protein_degradation_rate: f64,
    mrna_degradation_rate: f64,
    onset_time: f64,
) -> Option<(f64, KineticFitCoeffs)> {
    if times.len() != values.len() || times.is_empty() {
        return None;
    }
    let times = Array1::from_iter(times.iter().copied());
    let values = Array1::from_iter(values.iter().copied());
    let basis = times.mapv(|time| {
        kinetic_basis_value(
            time,
            protein_degradation_rate,
            mrna_degradation_rate,
            onset_time,
        )
    });
    if !basis.iter().all(|value| value.is_finite()) {
        return None;
    }
    let (baseline_intensity, expression_amplitude) = lstsq_affine(
        basis.as_slice().unwrap_or(&[]),
        values.as_slice().unwrap_or(&[]),
    )?;
    if !baseline_intensity.is_finite()
        || !expression_amplitude.is_finite()
        || expression_amplitude <= 0.0
    {
        return None;
    }
    let predicted = &basis * expression_amplitude + baseline_intensity;
    let residuals = &predicted - &values;
    let sse = residuals.mapv(|delta| delta * delta).sum();
    if !sse.is_finite() {
        return None;
    }
    Some((
        sse,
        KineticFitCoeffs {
            baseline_intensity,
            protein_degradation_rate,
            mrna_degradation_rate,
            onset_time,
            expression_amplitude,
        },
    ))
}

pub fn lstsq_affine(basis: &[f64], values: &[f64]) -> Option<(f64, f64)> {
    if basis.len() != values.len() || basis.is_empty() {
        return None;
    }
    let basis = Array1::from_iter(basis.iter().copied());
    let values = Array1::from_iter(values.iter().copied());
    let n = basis.len() as f64;
    let sum_1 = n;
    let sum_x = basis.sum();
    let sum_xx = basis.mapv(|value| value * value).sum();
    let sum_y = values.sum();
    let sum_xy = (&basis * &values).sum();

    let det = sum_1 * sum_xx - sum_x * sum_x;
    if det.abs() <= f64::EPSILON {
        return None;
    }
    let offset = (sum_y * sum_xx - sum_x * sum_xy) / det;
    let amplitude = (sum_1 * sum_xy - sum_x * sum_y) / det;
    Some((offset, amplitude))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn masked_roi_stats_uses_background_median() {
        // bg pixels [20, 100, 40] → median 40 (mean would be ~53.33)
        let frame = [10.0, 20.0, 100.0, 40.0];
        let mask = [true, false, false, false];
        let stats = masked_roi_stats(&frame, &mask).unwrap();
        assert_eq!(stats.area, 1);
        assert!((stats.intensity - 10.0).abs() < 1e-9);
        assert!((stats.background - 40.0).abs() < 1e-9);
        assert!((stats.corrected - (10.0 - 40.0)).abs() < 1e-9);
    }

    #[test]
    fn full_frame_roi_stats_uses_tenth_percentile_background() {
        let frame = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0];
        let stats = full_frame_roi_stats(&frame);
        assert_eq!(stats.area, 10);
        assert!((stats.intensity - 450.0).abs() < 1e-9);
        assert!((stats.background - quantile(&frame, 0.1)).abs() < 1e-9);
    }

    #[test]
    fn percentile_matches_quantile_scale() {
        let values = [1.0, 2.0, 3.0, 4.0];
        assert!((percentile(&values, 50.0) - quantile(&values, 0.5)).abs() < 1e-9);
    }

    #[test]
    fn quantile_linear_matches_numpy_linear() {
        // np.quantile([1, 4, 2, 3], q, method="linear")
        let values = [1.0, 4.0, 2.0, 3.0];
        assert!((quantile_linear(&values, 0.1) - 1.3).abs() < 1e-12);
        assert!((quantile_linear(&values, 0.5) - 2.5).abs() < 1e-12);
        assert!((quantile_linear(&values, 0.0) - 1.0).abs() < 1e-12);
        assert!((quantile_linear(&values, 1.0) - 4.0).abs() < 1e-12);
    }

    #[test]
    fn quantile_linear_skips_non_finite_values() {
        let values = [1.0, f64::NAN, 3.0, f64::INFINITY];
        assert!((quantile_linear(&values, 0.5) - 2.0).abs() < 1e-12);
        assert!(quantile_linear(&[f64::NAN], 0.5).is_nan());
        assert_eq!(quantile_linear(&[], 0.5), 0.0);
        assert_eq!(quantile_linear(&[1.0, 2.0], f64::NAN), 0.0);
    }

    #[test]
    fn quantile_linear_tied_plate_does_not_recurse() {
        // ndarray-stats quickselect removed one tied pivot per frame and
        // overflowed the analyze worker at ~10k identical trace points.
        let values = vec![3.5_f64; 50_000];
        assert!((quantile_linear(&values, 0.01) - 3.5).abs() < 1e-12);
        assert!((quantile_linear(&values, 0.99) - 3.5).abs() < 1e-12);
    }

    #[test]
    fn quantile_floor_uses_lower_index() {
        let sorted = [1.0, 2.0, 3.0, 4.0];
        assert_eq!(quantile_floor_sorted(&sorted, 0.75), 3.0);
    }

    #[test]
    fn kinetic_candidate_matches_scalar_formula() {
        let times = [0.0, 1.0, 2.0, 3.0];
        let values = [1.0, 2.0, 2.5, 2.0];
        let result = evaluate_kinetic_candidate(&times, &values, 0.5, 1.0, 0.0).unwrap();
        let predicted = times
            .iter()
            .map(|time| fitted_trace_value(*time, &result.1))
            .collect::<Vec<_>>();
        for (actual, expected) in values.iter().zip(predicted.iter()) {
            assert!((actual - expected).abs() < 0.5);
        }
    }

    #[test]
    fn trapezoidal_integral_matches_reference() {
        let times = [0.0, 1.0, 2.0];
        let values = [0.0, 2.0, 4.0];
        assert!((trapezoidal_integral(&times, &values) - 4.0).abs() < 1e-9);
    }

    #[test]
    fn lstsq_recovers_affine_coefficients() {
        let basis = vec![0.0, 1.0, 2.0, 3.0];
        let values = vec![1.0, 3.0, 5.0, 7.0];
        let (offset, amplitude) = lstsq_affine(&basis, &values).unwrap();
        assert!((offset - 1.0).abs() < 1e-9);
        assert!((amplitude - 2.0).abs() < 1e-9);
    }

    #[test]
    fn half_life_is_ln2_over_rate_not_reciprocal() {
        let rate = 0.1;
        let tau = half_life_minutes(rate);
        assert!((tau - std::f64::consts::LN_2 / rate).abs() < 1e-12);
        assert!((tau - 1.0 / rate).abs() > 1.0);
    }

    #[test]
    fn auc_from_fit_uses_half_life_eq4() {
        let expression_rate = 2.0;
        let mrna_lifetime = 30.0;
        let protein_lifetime = 120.0;
        let auc = auc_from_fit_half_lives(expression_rate, mrna_lifetime, protein_lifetime);
        let expected =
            std::f64::consts::LN_2.powi(2) * expression_rate * mrna_lifetime * protein_lifetime;
        assert!((auc - expected).abs() < 1e-12);
    }

    #[test]
    fn rates_and_amplitude_reconstruct_from_observables() {
        let protein_rate = 0.1;
        let mrna_rate = 0.5;
        let amplitude = 100.0;
        let protein_lifetime = half_life_minutes(protein_rate);
        let mrna_lifetime = half_life_minutes(mrna_rate);
        let expression_rate = amplitude * (mrna_rate - protein_rate);
        assert!((degradation_rate_per_minute(protein_lifetime) - protein_rate).abs() < 1e-12);
        assert!((degradation_rate_per_minute(mrna_lifetime) - mrna_rate).abs() < 1e-12);
        assert!(
            (expression_amplitude_from_observables(
                expression_rate,
                mrna_lifetime,
                protein_lifetime
            ) - amplitude)
                .abs()
                < 1e-12
        );
    }
}
