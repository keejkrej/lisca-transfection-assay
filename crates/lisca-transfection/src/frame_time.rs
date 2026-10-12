//! Minutes from trace column `t`.
//!
//! `t` stores `index.json` `timeIndices`: the source acquisition frame.
//! Interval is the time between consecutive acquisition frames, so
//! minutes = frame × interval. Labels `0, 6, 12` with a 10 minute interval
//! are `0, 60, 120`.

/// Map acquisition frame labels to minutes.
pub fn frame_minutes(labels: &[f64], interval: f64) -> Vec<f64> {
    labels.iter().map(|label| label * interval).collect()
}

#[cfg(test)]
mod tests {
    use super::frame_minutes;

    #[test]
    fn dense_labels_match_label_times_interval() {
        let labels = [0.0, 1.0, 2.0, 3.0];
        assert_eq!(frame_minutes(&labels, 10.0), vec![0.0, 10.0, 20.0, 30.0]);
    }

    #[test]
    fn stepped_acquisition_labels_scale_with_interval() {
        // A saved stack can keep source frames 0, 6, 12. Each step of 6 is
        // six intervals, not one.
        let labels = [0.0, 6.0, 12.0, 18.0];
        assert_eq!(frame_minutes(&labels, 10.0), vec![0.0, 60.0, 120.0, 180.0]);
    }

    #[test]
    fn unsorted_labels_keep_their_minutes() {
        let labels = [12.0, 0.0, 6.0];
        assert_eq!(frame_minutes(&labels, 10.0), vec![120.0, 0.0, 60.0]);
    }
}
