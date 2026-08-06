#![allow(dead_code)]
use std::io::{BufRead, BufReader, Lines, Read};
use thiserror::Error;

use crate::Strand;

#[derive(Debug, Error)]
pub enum BedError {
    #[error("BED entry with missing fields: {0}")]
    RecordTooShort(String),

    #[error("invalid integer field: {0}")]
    InvalidInteger(String),

    #[error("IO error: {0}")]
    IoError(#[source] std::io::Error),
}

#[derive(Debug)]
pub struct BedRecord {
    chr: String,
    start: usize,
    end: usize,
    id: Option<String>,
    score: Option<f32>,
    strand: Option<Strand>,
}
impl BedRecord {
    pub fn id(&self) -> Option<&str> {
        self.id.as_deref()
    }

    pub fn chr(&self) -> &str {
        &self.chr
    }

    pub fn start(&self) -> usize {
        self.start
    }

    pub fn end(&self) -> usize {
        self.end
    }

    pub fn strand(&self) -> Option<Strand> {
        self.strand
    }
}

pub struct BedReader<T> {
    buffer_lines: Lines<BufReader<T>>,
}
impl<T: Read> BedReader<T> {
    pub fn new(file: T) -> BedReader<T> {
        BedReader {
            buffer_lines: BufReader::new(file).lines(),
        }
    }
}
impl<T: Read> Iterator for BedReader<T> {
    type Item = Result<BedRecord, BedError>;

    fn next(&mut self) -> Option<Self::Item> {
        fn make_record(line: &str) -> Result<BedRecord, BedError> {
            let mut s = line.split_whitespace();

            Ok(BedRecord {
                chr: s
                    .next()
                    .map(|s| s.to_string())
                    .ok_or_else(|| BedError::RecordTooShort(line.to_owned()))?,
                start: s
                    .next()
                    .ok_or_else(|| BedError::RecordTooShort(line.to_owned()))?
                    .parse()
                    .map_err(|_| BedError::InvalidInteger(line.to_owned()))?,
                end: s
                    .next()
                    .ok_or_else(|| BedError::RecordTooShort(line.to_owned()))?
                    .parse()
                    .map_err(|_| BedError::InvalidInteger(line.to_owned()))?,
                id: s.next().map(|s| s.to_string()),
                score: s.next().map(|x| x.parse().unwrap_or_default()),
                strand: s.next().and_then(|x| x.try_into().ok()),
            })
        }

        loop {
            match self.buffer_lines.next()? {
                Err(e) => return Some(Err(BedError::IoError(e))),
                Ok(line) if line.starts_with('#') || line.is_empty() => continue,
                Ok(line) => return Some(make_record(&line)),
            }
        }
    }
}
