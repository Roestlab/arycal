use std::{env, path::{Path, PathBuf}};


/// Extract a basename with **single extension removal**.
/// Rules:
/// - Remove only the last `.<ext>` if there is a dot **not at position 0**.
/// - If the file name consists solely of dots (".", "..", "...") -> return unchanged.
/// - Hidden files with no additional dot (".bashrc" has one dot at pos 0 only) are returned unchanged.
/// - If the resulting stem would be empty (e.g. input was ".ext" or ".a") we keep the original.
/// - Returns empty String for paths that have no file name component (e.g. a directory path).
pub fn extract_basename<P: AsRef<Path>>(path: P) -> String {
    let p = path.as_ref();
    let file = match p.file_name().and_then(|s| s.to_str()) {
        Some(f) => f,
        None => return String::new(), // directory or invalid UTF-8 name
    };

    // All dots (".", "..", "...") -> leave unchanged
    if !file.is_empty() && file.chars().all(|c| c == '.') {
        return file.to_string();
    }

    // Find last dot
    if let Some(idx) = file.rfind('.') {
        if idx > 0 {
            // Candidate stem (part before last dot)
            let stem = &file[..idx];
            // If stem became empty (rare given idx>0) fall back to original
            if stem.is_empty() {
                return file.to_string();
            }
            return stem.to_string();
        }
        // idx == 0 => leading dot only (hidden file with no extension) -> keep whole
    }

    file.to_string()
}

/// Lexically extract the directory component (no filesystem access).
/// Rules:
/// - "/path/to/file.txt"  -> Some("/path/to")
/// - "/path/to/"          -> Some("/path/to")
/// - "file.txt"           -> Some(".")
/// - "."                  -> Some(".")
/// - ".."                 -> Some("..")
/// - "/"                  -> Some("/")
/// - ""                   -> None
pub fn extract_directory<P: AsRef<Path>>(path: P) -> Option<PathBuf> {
    let original = path.as_ref();

    if original.as_os_str().is_empty() {
        return None;
    }
    if original == Path::new("/") {
        return Some(PathBuf::from("/"));
    }
    if original == Path::new(".") {
        return Some(PathBuf::from("."));
    }
    if original == Path::new("..") {
        return Some(PathBuf::from(".."));
    }

    // Detect *syntactic* directory: trailing separator(s) (excluding root alone)
    let mut s = original.to_string_lossy().to_string();
    let had_trailing_sep = s.ends_with(std::path::MAIN_SEPARATOR);

    if had_trailing_sep {
        // Trim all trailing separators except leave a single "/" alone
        while s.len() > 1 && s.ends_with(std::path::MAIN_SEPARATOR) {
            s.pop();
        }
        // After trimming, s now names the directory itself; return it directly.
        if s.is_empty() {
            return Some(PathBuf::from("/")); // unlikely, but safe
        }
        // Re-handle special cases (like "." or ".." with trailing slash)
        if s == "." || s == ".." {
            return Some(PathBuf::from(s));
        }
        return Some(PathBuf::from(s));
    }

    // No trailing separator: treat last component as file name (or single component).
    let p = Path::new(&s);

    if let Some(parent) = p.parent() {
        if parent.as_os_str().is_empty() {
            return Some(PathBuf::from(".")); // single component (e.g. "file.txt")
        }
        return Some(parent.to_path_buf());
    }

    // Single component (not ".", ".." because handled above) -> "."
    Some(PathBuf::from("."))
}

/// Construct a filename from a directory, a "stem", an extra suffix, and an extension,
/// returning the full `PathBuf`.
///
/// # Behavior
///
/// - `stem` is the base name (without extension) you want to start from.
/// - `suffix` is appended **directly** after the stem (e.g. `"_targets"`, `"_decoys"`).
/// - `ext` is the *final* extension (without leading dot). If `ext` is empty, no dot or
///   extension is added. If `ext` **does** contain a leading dot (e.g. `".pqp"`), it is
///   normalized by stripping that dot.
/// - If `stem` is empty **or** consists only of dots (e.g. `""`, `"."`, `".."`), a
///   fallback stem `"output"` is used to avoid generating ambiguous or unsafe file names.
/// - The components are combined using `PathBuf::join`, so path separators are handled
///   correctly for the current platform.
///
/// # Examples
///
/// ```no_run
/// use std::path::PathBuf;
///
/// let p = add_suffix_and_ext("/data/out", "library", "_targets", "pqp");
/// assert_eq!(p, PathBuf::from("/data/out/library_targets.pqp"));
///
/// // Leading dot in extension is normalized
/// let p = add_suffix_and_ext("/data/out", "library", "_decoys", ".pqp");
/// assert_eq!(p, PathBuf::from("/data/out/library_decoys.pqp"));
///
/// // Empty extension -> no trailing dot
/// let p = add_suffix_and_ext("/data/out", "dataset", "_raw", "");
/// assert_eq!(p, PathBuf::from("/data/out/dataset_raw"));
///
/// // Empty / dot-only stem falls back to "output"
/// let p = add_suffix_and_ext("/data/out", "", "_targets", "pqp");
/// assert_eq!(p, PathBuf::from("/data/out/output_targets.pqp"));
/// let p = add_suffix_and_ext("/data/out", "..", "_targets", "pqp");
/// assert_eq!(p, PathBuf::from("/data/out/output_targets.pqp"));
/// ```
///
/// # Parameters
///
/// * `base_dir` - Directory into which the file will be placed. Accepts any type
///   implementing `AsRef<Path>` (e.g. `&str`, `PathBuf`, `&Path`).
/// * `stem` - The initial filename stem (without suffix / extension).
/// * `suffix` - A string appended immediately after `stem` (no separator added automatically).
/// * `ext` - Desired final extension; may optionally start with a dot; empty string means no extension.
///
/// # Returns
///
/// A `PathBuf` representing the combined path: `<base_dir>/<safe_stem><suffix>[.<ext>]`.
///
/// # Notes
///
/// This function performs only minimal sanitization (empty / dot-only stem). If you
/// need stricter validation (removing whitespace, path separators, etc.), sanitize
/// `stem` and/or `suffix` before calling.
pub fn add_suffix_and_ext<P: AsRef<Path>>(
    base_dir: P,
    stem: &str,
    suffix: &str,
    ext: &str
) -> PathBuf {
    let safe_stem = if stem.is_empty() || stem.chars().all(|c| c == '.') {
        "output"
    } else {
        stem
    };
    let filename = if ext.is_empty() {
        format!("{safe_stem}{suffix}")
    } else {
        format!("{safe_stem}{suffix}.{}", ext.trim_start_matches('.'))
    };
    base_dir.as_ref().join(filename)
}


/// Look for an executable named `binary_name` in:
/// 1. `extra_dir` (if provided),  
/// 2. every directory on the `PATH` environment variable.  
///
/// On Unix, only files with any of the owner/group/other execute bits set will match.
/// On Windows, we also try `.exe`, `.bat`, and `.cmd` suffixes automatically.
///
/// # Examples
///
/// ```rust
/// // First try `/usr/local/bin/mytool`, then fall back to $PATH
/// let exe = find_executable("mytool", Some(Path::new("/usr/local/bin")))
///     .expect("couldn't find `mytool`");
/// ```
pub fn find_executable(binary_name: &str, extra_dir: Option<&Path>) -> Option<PathBuf> {
    // 1) Try the extra directory first
    if let Some(dir) = extra_dir {
        let cand = dir.join(binary_name);
        if is_executable(&cand) {
            return Some(cand);
        }
        #[cfg(windows)]
        {
            for ext in &[".exe", ".bat", ".cmd"] {
                let cand = dir.join(format!("{}{}", binary_name, ext));
                if cand.is_file() {
                    return Some(cand);
                }
            }
        }
    }

    // 2) Then search each entry in $PATH
    if let Ok(paths) = env::var("PATH") {
        for dir in env::split_paths(&paths) {
            let cand = dir.join(binary_name);
            if is_executable(&cand) {
                return Some(cand);
            }
            #[cfg(windows)]
            for ext in &[".exe", ".bat", ".cmd"] {
                let cand = dir.join(format!("{}{}", binary_name, ext));
                if cand.is_file() {
                    return Some(cand);
                }
            }
        }
    }

    None
}

/// Returns `true` if `path` exists, is a file, and (on Unix) has at least one exec bit set.
fn is_executable(path: &Path) -> bool {
    if !path.is_file() {
        return false;
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        match path.metadata() {
            Ok(meta) => meta.permissions().mode() & 0o111 != 0,
            Err(_) => false,
        }
    }
    #[cfg(not(unix))]
    {
        // On Windows and others, existence is enough
        true
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_basename() {
        assert_eq!(extract_basename("file.txt"), "file");
        assert_eq!(extract_basename("archive.tar.gz"), "archive.tar");
        assert_eq!(extract_basename("no_extension"), "no_extension");
        assert_eq!(extract_basename(".hiddenfile"), ".hiddenfile");
    }

    #[test]
    fn test_extract_directory() {
        assert_eq!(extract_directory("/path/to/file.txt"), Some(PathBuf::from("/path/to")));
        assert_eq!(extract_directory("/path/to/"), Some(PathBuf::from("/path/to")));
        assert_eq!(extract_directory("file.txt"), Some(PathBuf::from(".")));
        assert_eq!(extract_directory("."), Some(PathBuf::from(".")));
        assert_eq!(extract_directory(".."), Some(PathBuf::from("..")));
        assert_eq!(extract_directory("/"), Some(PathBuf::from("/")));
    }
}