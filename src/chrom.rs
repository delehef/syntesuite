//! A parser for the ChromTable format:
//! ```
//! Chrom [TAB] Start [TAB] Stop [TAB] Strand [TAB] geneid
//! ```

use std::io::{BufRead, BufReader, Lines, Read};

use thiserror::Error;

use crate::Strand;

#[derive(Debug, Error)]
pub enum ChromError {
    #[error("ChromTable entry with missing fields: {0}")]
    RecordTooShort(String),
    #[error("Unrecognized strand format: {0}")]
    UnknownStrand(String),
}

#[derive(Debug)]
pub struct ChromRecord {
    chr: String,
    start: usize,
    end: usize,
    id: String,
    strand: Strand,
}

impl ChromRecord {
    pub fn id(&self) -> &str {
        &self.id
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

    pub fn strand(&self) -> Strand {
        self.strand
    }
}

pub struct ChromReader<T> {
    buffer_lines: Lines<BufReader<T>>,
}
impl<T: Read> ChromReader<T> {
    pub fn new(file: T) -> ChromReader<T> {
        ChromReader {
            buffer_lines: BufReader::new(file).lines(),
        }
    }
}
impl<T: Read> Iterator for ChromReader<T> {
    type Item = Result<ChromRecord, ChromError>;

    fn next(&mut self) -> Option<Self::Item> {
        fn make_record(line: &str) -> Result<ChromRecord, ChromError> {
            let mut s = line.split('\t');

            Ok(ChromRecord {
                chr: s
                    .next()
                    .map(|s| s.to_string())
                    .ok_or_else(|| ChromError::RecordTooShort(line.to_owned()))?,
                start: s
                    .next()
                    .ok_or_else(|| ChromError::RecordTooShort(line.to_owned()))?
                    .parse()
                    .unwrap(),
                end: s
                    .next()
                    .ok_or_else(|| ChromError::RecordTooShort(line.to_owned()))?
                    .parse()
                    .unwrap(),
                strand: s
                    .next()
                    .ok_or_else(|| ChromError::RecordTooShort(line.to_owned()))?
                    .try_into()
                    .map_err(|_| ChromError::UnknownStrand(line.to_owned()))?,
                id: s
                    .next()
                    .ok_or_else(|| ChromError::RecordTooShort(line.to_owned()))?
                    .to_owned(),
            })
        }

        self.buffer_lines
            .by_ref()
            .map(|l| l.unwrap())
            .find(|line| !line.starts_with('#') && !line.is_empty())
            .map(|l| make_record(&l))
    }
}
