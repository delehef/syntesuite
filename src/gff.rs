#![allow(dead_code)]
use std::collections::HashMap;
use std::io::prelude::*;
use std::io::{BufReader, Lines};
use thiserror::Error;

use crate::{Phase, Strand};

#[derive(Debug, Error)]
pub enum GffError {
    #[error("GFF entry with missing fields: {0}")]
    RecordTooShort(String),

    #[error("attribute entry contains more than one `=`: {0}")]
    IncorrectAttribute(String),
}

/// A key to a GFF3 record attribute, as defined in http://gmod.org/wiki/GFF3
#[derive(PartialEq, Eq, Hash, Clone, Debug)]
pub enum Key {
    ID,
    Name,
    Alias,
    Parent,
    Target,
    Gap,
    DerivesFrom,
    Note,
    Dbxref,
    OntologyTerm,
    K(String),
}
impl From<&str> for Key {
    fn from(s: &str) -> Self {
        match s.to_lowercase().as_ref() {
            "id" => Key::ID,
            "name" => Key::Name,
            "alias" => Key::Alias,
            "parent" => Key::Parent,
            "target" => Key::Target,
            "gap" => Key::Gap,
            "derives_from" => Key::DerivesFrom,
            "note" => Key::Note,
            "dbxref" => Key::Dbxref,
            "ontology_term" => Key::OntologyTerm,
            _ => Key::K(s.to_string()),
        }
    }
}

type Attributes = HashMap<Key, Vec<String>>;
#[derive(Debug)]
pub struct GffRecord {
    chr: String,
    source: Option<String>,
    class: Option<String>,
    start: usize,
    end: usize,
    score: Option<f32>,
    strand: Option<Strand>,
    phase: Option<Phase>,
    attributes: Attributes,
}
impl GffRecord {
    pub fn chr(&self) -> &str {
        &self.chr
    }
    pub fn id(&self) -> Option<&str> {
        self.attributes
            .get(&Key::ID)
            .and_then(|x| x.first())
            .map(|x| x.as_str())
    }
    pub fn source(&self) -> Option<&String> {
        self.source.as_ref()
    }
    pub fn class(&self) -> Option<&String> {
        self.class.as_ref()
    }
    pub fn start(&self) -> usize {
        self.start
    }
    pub fn end(&self) -> usize {
        self.end
    }
    pub fn score(&self) -> Option<f32> {
        self.score
    }
    pub fn strand(&self) -> Option<Strand> {
        self.strand
    }
    pub fn phase(&self) -> Option<Phase> {
        self.phase
    }
    pub fn attributes(&self) -> &Attributes {
        &self.attributes
    }
    /// If the record has a Parent attribute, return its first value
    pub fn parent(&self) -> Option<&String> {
        self.parents().and_then(|v| v.first())
    }
    /// If the record has a Parent attribute, return all its values
    pub fn parents(&self) -> Option<&Vec<String>> {
        self.attributes.get(&Key::Parent)
    }
    /// If the record has a Target attribute, return its first value
    pub fn target(&self) -> Option<&String> {
        self.targets().and_then(|v| v.first())
    }
    /// If the record has a Target attribute, return all its values
    pub fn targets(&self) -> Option<&Vec<String>> {
        self.attributes.get(&Key::Target)
    }
}

pub struct GffReader<T> {
    buffer_lines: Lines<BufReader<T>>,
}
impl<T: Read> GffReader<T> {
    pub fn new(file: T) -> GffReader<T> {
        GffReader {
            buffer_lines: BufReader::new(file).lines(),
        }
    }
}
impl<T: Read> Iterator for GffReader<T> {
    type Item = Result<GffRecord, GffError>;

    fn next(&mut self) -> Option<Self::Item> {
        fn make_record(line: &str) -> Result<GffRecord, GffError> {
            let mut s = line.split('\t');

            Ok(GffRecord {
                chr: s
                    .next()
                    .map(|s| s.to_string())
                    .ok_or_else(|| GffError::RecordTooShort(line.to_owned()))?,
                source: s
                    .next()
                    .map(|x| if x == "." { None } else { Some(x.to_string()) })
                    .ok_or_else(|| GffError::RecordTooShort(line.to_owned()))?,
                class: s
                    .next()
                    .ok_or_else(|| GffError::RecordTooShort(line.to_owned()))
                    .map(|x| if x == "." { None } else { Some(x.to_string()) })?,
                start: s
                    .next()
                    .ok_or_else(|| GffError::RecordTooShort(line.to_owned()))?
                    .parse()
                    .map_err(|_| GffError::RecordTooShort(line.to_owned()))?,
                end: s
                    .next()
                    .ok_or_else(|| GffError::RecordTooShort(line.to_owned()))?
                    .parse()
                    .map_err(|_| GffError::RecordTooShort(line.to_owned()))?,
                score: s
                    .next()
                    .map(|x| {
                        if x == "." {
                            None
                        } else {
                            Some(x.parse().unwrap())
                        }
                    })
                    .ok_or_else(|| GffError::RecordTooShort(line.to_owned()))?,
                strand: s
                    .next()
                    .map(|x| {
                        if x == "." {
                            None
                        } else {
                            Some(x.try_into().unwrap())
                        }
                    })
                    .ok_or_else(|| GffError::RecordTooShort(line.to_owned()))?,
                phase: s
                    .next()
                    .map(|x| {
                        if x == "." {
                            None
                        } else {
                            Some(x.try_into().unwrap())
                        }
                    }) // TODO remove the unwrap
                    .ok_or_else(|| GffError::RecordTooShort(line.to_owned()))?,
                attributes: s
                    .next()
                    .ok_or_else(|| GffError::RecordTooShort(line.to_owned()))?
                    .split(';')
                    .map(|pair| {
                        let s = pair.split('=').collect::<Vec<_>>();
                        if s.len() != 2 {
                            return Err(GffError::IncorrectAttribute(pair.to_string()));
                        }
                        Ok((
                            Key::from(s[0]),
                            s[1].split(',').map(|x| x.to_string()).collect(),
                        ))
                    })
                    .collect::<Result<Attributes, GffError>>()?,
            })
        }

        self.buffer_lines
            .by_ref()
            .map(|l| l.unwrap())
            .find(|line| !line.starts_with('#') && !line.is_empty())
            .map(|l| make_record(&l))
    }
}
