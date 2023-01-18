use colored::Colorize;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum FileError {
    #[error("failed to open {}", .filename.bright_yellow().bold())]
    CannotOpen {
        source: std::io::Error,
        filename: String,
    },

    #[error("{} not found", .0.bright_yellow().bold())]
    NotFound(String),

    #[error("while creating {filename}")]
    WhileCreating {
        source: std::io::Error,
        filename: String,
    },

    #[error("invalid filename: {}", .0.yellow().bold())]
    InvalidFilename(String),
}

#[derive(Error, Debug)]
pub enum DataError {
    #[error("ID {} not found in the specified database", .0.yellow().bold())]
    UnknownId(String),

    #[error("failed to connect to database {}", .filename.yellow().bold())]
    FailedToConnect {
        source: rusqlite::Error,
        filename: String,
    },
}
