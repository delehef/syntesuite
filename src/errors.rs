use colored::Colorize;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum FileError {
    #[error("failed to open {}", .filename.bright_yellow().bold())]
    CannotOpen {
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

    #[error("inline gene books can not be accessed mutably")]
    ImmutableBook,
}

#[derive(Error, Debug)]
pub enum ParseError {
    #[error("wrongly formatted GFF file")]
    GffError(crate::gff::GffError),

    #[error("wrongly formatted BED file")]
    BedError(crate::bed::BedError),

    #[error("invalid phase value: {0}")]
    InvalidPhase(String),

    #[error("invalid strand value: {0}")]
    InvalidStrand(String),
}
