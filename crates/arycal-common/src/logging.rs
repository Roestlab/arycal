use std::sync::{Arc, Mutex, mpsc};
use std::thread;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::ops::Range;
use tqdm::Tqdm;
use tqdm::tqdm;

/// A thread-safe progress bar implementation using `tqdm`.
///
/// This struct manages a progress bar that updates asynchronously in a separate thread.
/// It ensures safe concurrent updates by using an atomic counter and a message-passing
/// channel to send update signals.
///
/// # Fields
/// * `total` - The total number of steps in the progress bar.
/// * `progress` - A shared, mutex-protected `tqdm` progress bar.
/// * `count` - An atomic counter to track progress updates.
/// * `sender` - A channel sender to send progress update messages.
/// * `progress_thread` - An optional handle for the background progress update thread.
/// * `description` - A description displayed alongside the progress bar.
pub struct Progress {
    total: usize,
    progress: Arc<Mutex<Tqdm<Range<usize>>>>, 
    count: AtomicUsize,               // Atomic counter for tracking progress
    sender: mpsc::Sender<usize>,      // Channel to send updates
    progress_thread: Option<thread::JoinHandle<()>>, // Background thread to update tqdm
    description: String,              // Description for the progress bar
}

impl Progress {
    /// Creates a new progress bar with the given total steps and description.
    ///
    /// This function initializes a `tqdm` progress bar and spawns a background thread
    /// to update the bar asynchronously.
    ///
    /// # Arguments
    /// * `total` - The total number of steps for the progress bar.
    /// * `description` - A string describing the progress bar.
    ///
    /// # Returns
    /// * A new `Progress` instance.
    ///
    /// # Example
    /// ```
    /// let progress = Progress::new(100, "Processing data");
    /// ```
    pub fn new(total: usize, description: &str) -> Self {
        let progress = Arc::new(Mutex::new(tqdm(0..total).desc(Some(description)))); // Initialize Tqdm
        let count = AtomicUsize::new(0);

        let (tx, rx) = mpsc::channel();
        let progress_clone = Arc::clone(&progress);

        // Spawn a thread to handle progress updates
        let handle = thread::spawn(move || {
            for _ in rx {
                let _ = progress_clone.lock().unwrap().pbar.update(1); // Always update by 1
            }
        });

        Self {
            total,
            progress,
            count,
            sender: tx,
            progress_thread: Some(handle),
            description: description.to_string(),
        }
    }

    /// Increments the progress counter safely and updates the progress bar.
    ///
    /// This function uses an atomic counter to track updates and sends an update signal
    /// to the progress thread via a channel.
    ///
    /// If the total progress exceeds the expected `total`, a warning is printed to avoid overflow.
    ///
    /// # Example
    /// ```
    /// let progress = Progress::new(100, "Loading");
    /// progress.inc();
    /// ```
    pub fn inc(&self) {
        let new_count = self.count.fetch_add(1, Ordering::AcqRel) + 1;
    
        if new_count > self.total {
            println!("⚠️ WARNING: Extra update detected! Skipping...");
            return; // Prevent overflow
        }

        let _ = self.sender.send(1); // Always send 1 instead of new_count
    }

    /// Updates the progress bar's description.
    pub fn update_description(&self, new_desc: &str) {
        let mut progress = self.progress.lock().unwrap();
        progress.set_desc(Some(new_desc)); // Update the description dynamically
    }

    /// Finalizes the progress bar by ensuring all updates are completed.
    ///
    /// This function drops the sender to close the channel and waits for the background thread
    /// to finish execution, ensuring no updates are lost.
    ///
    /// # Example
    /// ```
    /// let progress = Progress::new(100, "Downloading files");
    /// for _ in 0..100 {
    ///     progress.inc();
    /// }
    /// progress.finish();
    /// ```
    pub fn finish(self) {
        drop(self.sender); // Close the channel to signal the thread to finish
        if let Some(handle) = self.progress_thread {
            let _ = handle.join(); // Wait for the progress thread to exit
        }
    }
}