//! On-disk folder names for a LiSCA workspace.
//!
//! Canonical schema is lisca (`docs/analysis/schema.md`). This crate must not
//! depend on the `lisca` crate (that crate depends on this one → cycle). Python
//! `transfection` imports `lisca` for the same paths.

use std::path::{Path, PathBuf};

pub const ROI_DIRNAME: &str = "roi";
pub const MASK_DIRNAME: &str = "mask";
pub const BBOX_DIRNAME: &str = "bbox";
pub const ANALYSIS_DIRNAME: &str = "analysis";
pub const RESULTS_DIRNAME: &str = "results";

pub fn roi_dir(workspace: &Path) -> PathBuf {
    workspace.join(ROI_DIRNAME)
}

pub fn roi_pos_dir(workspace: &Path, pos: u32) -> PathBuf {
    roi_dir(workspace).join(format!("Pos{pos}"))
}

pub fn mask_pos_dir(workspace: &Path, pos: u32) -> PathBuf {
    workspace.join(MASK_DIRNAME).join(format!("Pos{pos}"))
}

pub fn analysis_dir(workspace: &Path) -> PathBuf {
    workspace.join(ANALYSIS_DIRNAME)
}

pub fn analysis_pos_dir(workspace: &Path, pos: u32) -> PathBuf {
    analysis_dir(workspace).join(format!("Pos{pos}"))
}

pub fn results_dir(workspace: &Path) -> PathBuf {
    workspace.join(RESULTS_DIRNAME)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn folder_names_match_lisca_tree() {
        let root = Path::new("/ws");
        assert_eq!(roi_pos_dir(root, 7), Path::new("/ws/roi/Pos7"));
        assert_eq!(analysis_pos_dir(root, 7), Path::new("/ws/analysis/Pos7"));
        assert_eq!(results_dir(root), Path::new("/ws/results"));
        assert_eq!(mask_pos_dir(root, 7), Path::new("/ws/mask/Pos7"));
        let _ = BBOX_DIRNAME;
    }
}
