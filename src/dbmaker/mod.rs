use anyhow::*;
use colored::Colorize;
use flate2::bufread::GzDecoder;
use log::*;
use regex::Regex;
use rusqlite::Connection;
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufRead, BufReader, Seek},
};
use thiserror::*;

use crate::{
    bed,
    errors::{DataError, FileError, ParseError},
    gff, Record, Strand,
};

#[derive(Error, Debug)]
pub enum Error {
    #[error("{} is not a valid regex", .re.yellow().bold())]
    InvalidRegex { source: regex::Error, re: String },

    #[error("capture group {} missing in {}", .cap.yellow().bold(), .re.blue().bold())]
    MissingCaptureGroup { cap: String, re: String },

    #[error("the provided regex did not match any species name in {}", .0.bold().yellow())]
    SpeciesNotFound(String),

    #[error("the provided regex did not match any ID in {}", .0.bold().yellow())]
    IdNotFound(String),

    #[error("record {} has no ID", .0.bold().yellow())]
    RecordWithoutId(String),
}

struct Annotation {
    id: String,
    dir: Strand,
    start: usize,
    stop: usize,
    ancestral_id: usize,
}

fn parse_family(
    f: &str,
    current_ancestral_id: &mut usize,
    id2ancestral: &mut HashMap<String, usize>,
) -> Result<()> {
    trace!("Processing {}", f.bright_white().bold());
    for l in BufReader::new(File::open(f).map_err(|e| FileError::CannotOpen {
        source: e,
        filename: f.to_owned(),
    })?)
    .lines()
    {
        for id in l?.split_whitespace() {
            id2ancestral.insert(id.into(), *current_ancestral_id);
        }
    }
    *current_ancestral_id += 1;

    Ok(())
}

fn parse_genome_gff3(f: &str) -> Result<Box<dyn Iterator<Item = Result<Record, ParseError>>>> {
    let mut f = File::open(f).map_err(|e| FileError::CannotOpen {
        source: e,
        filename: f.to_owned(),
    })?;
    let gz = GzDecoder::new(BufReader::new(f.try_clone().unwrap()));

    Ok(match gz.header() {
        Some(_) => Box::new(
            gff::GffReader::new(gz)
                .map(|r| r.map(|r| r.into()).map_err(|e| ParseError::GffError(e))),
        ),
        None => {
            f.rewind()?;
            Box::new(
                gff::GffReader::new(BufReader::new(f))
                    .map(|r| r.map(|r| r.into()).map_err(|e| ParseError::GffError(e))),
            )
        }
    })
}

fn parse_genome_bed(f: &str) -> Result<Box<dyn Iterator<Item = Result<Record, ParseError>>>> {
    let mut f = File::open(f).map_err(|e| FileError::CannotOpen {
        source: e,
        filename: f.to_owned(),
    })?;
    let gz = GzDecoder::new(BufReader::new(f.try_clone().unwrap()));

    Ok(match gz.header() {
        Some(_) => Box::new(
            bed::BedReader::new(gz)
                .map(|r| r.map(|r| r.into()).map_err(|e| ParseError::BedError(e))),
        ),
        None => {
            f.rewind()?;
            Box::new(
                bed::BedReader::new(BufReader::new(f))
                    .map(|r| r.map(|r| r.into()).map_err(|e| ParseError::BedError(e))),
            )
        }
    })
}

fn parse_file(
    filename: &str,
    species_pattern: &str,
) -> Result<(
    String,
    impl Iterator<Item = Result<crate::Record, ParseError>>,
)> {
    info!("Processing {}", filename.bright_white().bold());
    let species_regex = Regex::new(species_pattern).map_err(|e| Error::InvalidRegex {
        source: e,
        re: species_pattern.to_string(),
    })?;
    if !species_regex
        .capture_names()
        .any(|n| n.map(|n| n == "species").unwrap_or(false))
    {
        return Err(Error::MissingCaptureGroup {
            cap: "species".into(),
            re: species_pattern.into(),
        }
        .into());
    }
    let species = species_regex
        .captures(
            std::path::Path::new(filename)
                .file_name()
                .ok_or_else(|| FileError::InvalidFilename(filename.to_string()))?
                .to_str()
                .ok_or_else(|| FileError::InvalidFilename(filename.to_string()))?,
        )
        .ok_or_else(|| Error::SpeciesNotFound(filename.to_string()))?["species"]
        .to_string();
    info!("Species: {}", species);
    let records = if filename.ends_with("gff")
        || filename.ends_with("gff3")
        || filename.ends_with("gff.gz")
        || filename.ends_with("gff3.gz")
    {
        parse_genome_gff3(filename)?
    } else if filename.ends_with("bed") || filename.ends_with("bed.gz") {
        parse_genome_bed(filename)?
    } else {
        bail!(
            "unable to process {}: unknown filetype",
            filename.yellow().bold()
        )
    };
    Ok((species, records))
}

fn parse_genome(
    f: &str,
    species_pattern: &str,
    id_type: &str,
    id_pattern: &str,
    genomes: &mut HashMap<String, HashMap<String, Vec<Annotation>>>,
    id2ancestral: &HashMap<String, usize>,
) -> Result<()> {
    let id_regex = Regex::new(id_pattern).map_err(|e| Error::InvalidRegex {
        source: e,
        re: id_pattern.to_string(),
    })?;
    if !id_regex
        .capture_names()
        .any(|n| n.map(|n| n == "id").unwrap_or(false))
    {
        return Err(Error::MissingCaptureGroup {
            cap: "id".into(),
            re: species_pattern.into(),
        }
        .into());
    }

    let mut seen = HashSet::new();
    let (species, records) = parse_file(f, species_pattern)?;
    for record in records {
        let record = record?;
        if record.is_class(id_type) {
            let id = record.id().ok_or_else(|| {
                Error::RecordWithoutId(format!(
                    "{}:{}-{}",
                    record.chr(),
                    record.start(),
                    record.end()
                ))
            })?;
            let id = id_regex
                .captures(id)
                .ok_or_else(|| Error::IdNotFound(id.into()))?["id"]
                .to_string();
            if let Some(ancestral_id) = id2ancestral.get(&id) {
                if seen.insert(id.clone()) {
                    genomes
                        .entry(species.clone())
                        .or_default()
                        .entry(record.chr().into())
                        .or_default()
                        .push(Annotation {
                            id: id.to_string(),
                            dir: record.strand(),
                            start: record.start(),
                            stop: record.end(),
                            ancestral_id: *ancestral_id,
                        });
                }
            } else {
                debug!("Skipping ID {} not found in families", id.bold().yellow());
            }
            trace!(
                "{}:{}/{} - {}",
                id,
                record.chr(),
                record.start(),
                record.end()
            );
        }
    }

    if let Some(genome) = genomes.get_mut(&species) {
        for (_, ids) in genome.iter_mut() {
            ids.sort_by_key(|a| a.start);
        }
    } else {
        warn!("{} appears to be empty", species.yellow().bold());
    }
    Ok(())
}

pub fn db_from_files(
    families: &[String],
    gffs: &[String],
    db_file: &str,
    species_pattern: &str,
    id_type: &str,
    id_pattern: &str,
    window: isize,
) -> Result<()> {
    let mut current_ancestral_id = 1;
    let mut id2ancestral = HashMap::new();
    info!("Parsing families...");
    for name in families.iter() {
        let path = std::path::Path::new(name);
        if path.is_dir() {
            for f in path
                .read_dir()
                .with_context(|| anyhow!("while reading {}", name))?
                .map(|e| {
                    e.map(|e| e.path().to_str().unwrap().to_owned())
                        .map_err(|_| todo!())
                })
            {
                parse_family(
                    f.unwrap().as_str(),
                    &mut current_ancestral_id,
                    &mut id2ancestral,
                )?;
            }
        } else {
            parse_family(name, &mut current_ancestral_id, &mut id2ancestral)?;
        }
    }

    info!("Parsing GFF3s...");
    let mut genomes = HashMap::new();
    for name in gffs.iter() {
        let path = std::path::Path::new(name);
        if path.is_dir() {
            for f in path
                .read_dir()
                .with_context(|| anyhow!("while reading {}", name))?
                .map(|e| {
                    e.map(|e| e.path().to_str().unwrap().to_owned())
                        .map_err(|_| todo!())
                })
            {
                parse_genome(
                    f.unwrap().as_str(),
                    species_pattern,
                    id_type,
                    id_pattern,
                    &mut genomes,
                    &id2ancestral,
                )?;
            }
        } else {
            parse_genome(
                name,
                species_pattern,
                id_type,
                id_pattern,
                &mut genomes,
                &id2ancestral,
            )?;
        }
    }

    info!("Creating database...");
    let mut conn = Connection::open(db_file).map_err(|e| DataError::FailedToConnect {
        source: e,
        filename: db_file.into(),
    })?;
    conn.execute("DROP TABLE IF EXISTS genomes;", [])
        .with_context(|| "while dropping table")?;
    conn.execute(
        "CREATE TABLE genomes (
            species text, chr text, ancestral_id integer, id text,
            start integer, stop integer, direction char,
            left_tail_ids text, right_tail_ids text
        )",
        [],
    )
    .with_context(|| "while creating database")?;
    info!("Filling database...");
    conn.execute("pragma temp_store = memory;", [])
        .with_context(|| "while setting temp_store")?;
    for (species, genome) in genomes.iter() {
        debug!("Inserting {}", species.bold());
        for (chr, ids) in genome.iter() {
            trace!("Inserting {}", chr.bold());
            let tx = conn.transaction()?;
            for (j, id) in ids.iter().enumerate() {
                let j = j as isize;
                let i = (0.max(j - window)) as usize;
                let k = ((ids.len() as isize - 1).min(j + window)) as usize;
                let left_landscape_ids = ids[i..j as usize]
                    .iter()
                    .map(|a| format!("{}{}", a.dir, a.ancestral_id))
                    .collect::<Vec<_>>();
                let right_landscape_ids = ids[j as usize + 1..=k]
                    .iter()
                    .map(|a| format!("{}{}", a.dir, a.ancestral_id))
                    .collect::<Vec<_>>();
                let insert = format!(
                    "INSERT INTO genomes (species, chr, ancestral_id, id, start, stop, direction, left_tail_ids, right_tail_ids) VALUES ('{}','{}','{}','{}','{}','{}','{}','{}','{}')",
                    species,
                    chr,
                    id.ancestral_id,
                    id.id,
                    id.start,
                    id.stop,
                    String::from(id.dir),
                    left_landscape_ids
                        .into_iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<_>>()
                        .join("."),
                    right_landscape_ids
                        .into_iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<_>>()
                        .join("."),
                );
                tx.execute(&insert, [])?;
            }
            tx.commit()?;
        }
    }

    info!("Creating DB indices...");
    conn.execute_batch(
        "CREATE INDEX genomes_species ON genomes(species);
         CREATE INDEX genomes_chr     ON genomes(chr);
         CREATE INDEX genomes_id      ON genomes(id);
         CREATE INDEX genomes_start   ON genomes(start);",
    )
    .with_context(|| "while creating indices")?;

    Ok(())
}
