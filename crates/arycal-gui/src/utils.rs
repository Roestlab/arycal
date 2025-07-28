use std::{
    env,
    error::Error,
    fs,
    io::{self, Cursor, Read},
    path::{Path, PathBuf},
    process::Command,
};

use flate2::read::GzDecoder;
use reqwest::blocking::Client;
use serde::Deserialize;
use tar::Archive;
use zip::ZipArchive;

/// A subset of the GitHub release JSON we need.
#[derive(Deserialize)]
struct Release {
    assets: Vec<Asset>,
}

#[derive(Deserialize)]
struct Asset {
    name: String,
    browser_download_url: String,
}

/// Download the latest `arycal` release for the current platform, unpack it,
/// and install the `arycal` binary into `~/.local/bin` (or on Windows `%USERPROFILE%\.local\bin`).
pub fn install_latest_arycal() -> Result<(), Box<dyn Error>> {
    // 1) Determine our target triple substring
    let arch = env::consts::ARCH; // e.g. "x86_64"
    let os   = env::consts::OS;   // e.g. "linux", "macos", "windows"
    let target = format!("{}-{}", arch, os);

    let target = if cfg!(windows) {
        // On Windows, we use "x86_64-pc-windows-msvc" or similar
        "arycal-x86_64-pc-windows-msvc.zip"
    } else if cfg!(target_os = "macos") {
        // On macOS, we use "x86_64-apple-darwin" or similar
        "arycal-x86_64-apple-darwin.tar.gz"
    } else {
        // For Linux and other Unix-like systems, we use "x86_64-unknown-linux-gnu" or similar
        "arycal-x86_64-unknown-linux-musl.tar.gz"
    };

    // 2) Fetch the latest release JSON
    let client = Client::builder()
        .user_agent("arycal-installer")
        .build()?;
    let release: Release = client
        .get("https://api.github.com/repos/singjc/arycal/releases/latest")
        .send()?
        .error_for_status()?
        .json()?;

    // 3) Find the matching asset
    let asset = release
        .assets
        .into_iter()
        .find(|a| a.name.contains(&target))
        .ok_or_else(|| format!("No release asset found for target `{}`", target))?;

    println!("Found asset `{}` → {}", asset.name, asset.browser_download_url);

    // 4) Download into memory
    let resp = client
        .get(&asset.browser_download_url)
        .send()?
        .error_for_status()?;
    let buf = resp.bytes()?.to_vec();

    // 5) Prepare install directory
    let install_dir: PathBuf = if cfg!(windows) {
        dirs::home_dir()
            .ok_or("Could not locate home directory")?
            .join(".local")
            .join("bin")
    } else {
        dirs::home_dir()
            .ok_or("Could not locate home directory")?
            .join(".local")
            .join("bin")
    };
    fs::create_dir_all(&install_dir)?;

    // 6) Unpack
    if asset.name.ends_with(".tar.gz") {
        let tar = GzDecoder::new(Cursor::new(buf));
        let mut archive = Archive::new(tar);
        for entry in archive.entries()? {
            let mut entry = entry?;
            let path = entry.path()?;
            // Look for the `arycal` binary at the top level
            if let Some(file_name) = path.file_name() {
                if file_name == "arycal" || file_name == "arycal.exe" {
                    let out_path = install_dir.join(file_name);
                    println!("Installing `{}` to `{}`", file_name.to_string_lossy(), out_path.display());
                    entry.unpack(&out_path)?;
                    // make executable on Unix
                    #[cfg(unix)]
                    {
                        use std::os::unix::fs::PermissionsExt;
                        let mut perms = fs::metadata(&out_path)?.permissions();
                        perms.set_mode(0o755);
                        fs::set_permissions(&out_path, perms)?;
                    }
                    break;
                }
            }
        }
    } else if asset.name.ends_with(".zip") {
        let reader = Cursor::new(buf);
        let mut zip = ZipArchive::new(reader)?;
        for i in 0..zip.len() {
            let mut file = zip.by_index(i)?;
            let name = file.name().rsplit('/').next().unwrap_or("");
            if name == "arycal.exe" {
                let out_path = install_dir.join(name);
                println!("Installing `{}` to `{}`", name, out_path.display());
                let mut outfile = fs::File::create(&out_path)?;
                io::copy(&mut file, &mut outfile)?;
                break;
            }
        }
    } else {
        return Err("Unknown asset archive format".into());
    }

    println!("✔ Installed into `{}`", install_dir.display());
    println!("⚠️  Make sure `{}` is on your PATH", install_dir.display());
    Ok(())
}

/// Get system based install directory for binaries.
/// Returns `~/.local/bin` on Unix-like systems or `%USERPROFILE%\.local\bin` on Windows.
/// If the directory does not exist, it will be created.
/// Returns the path to the install directory.
/// # Errors
/// If the home directory cannot be determined or if the directory cannot be created.
pub fn get_install_dir() -> Result<PathBuf, Box<dyn Error>> {
    let home = dirs::home_dir().ok_or("Could not locate home directory")?;
    let install_dir = if cfg!(windows) {
        home.join(".local").join("bin")
    } else {
        home.join(".local").join("bin")
    };
    fs::create_dir_all(&install_dir)?;
    Ok(install_dir)
}

#[cfg(test)]
mod tests {
    use super::install_latest_arycal;

    /// This is a live‐network test.  It will actually hit github.com,
    /// download, unpack and install the latest arycal binary into ~/.local/bin.
    /// Marked `ignore` so it doesn’t run in your normal `cargo test`.
    #[test]
    #[ignore]
    fn test_install_latest_arycal() {
        // If this returns Err, something went wrong fetching or unpacking.
        install_latest_arycal().expect("failed to install latest arycal");
    }
}