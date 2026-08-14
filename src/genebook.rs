use anyhow::*;
use derive_more::{Deref, Display, From, Into};
use log::*;
use rusqlite::Connection;
use std::collections::HashMap;
use std::num::TryFromIntError;
use std::sync::Mutex;

use crate::{dbmaker::LANDSCAPE_DELIMITER, errors, Strand};

#[derive(Default, Clone, Copy, From, Into, Deref, Display, PartialEq, Eq, Hash)]
#[display("{}", self.0)]
pub struct FamilyId(usize);
impl TryFrom<i64> for FamilyId {
    type Error = TryFromIntError;

    fn try_from(value: i64) -> Result<Self, Self::Error> {
        usize::try_from(value).map(FamilyId)
    }
}

#[allow(dead_code)]
pub enum GeneBook {
    InMemory {
        genes: HashMap<String, Gene>,
        species: Vec<String>,
    },
    Inline {
        conn: Mutex<Connection>,
        window: usize,
        id_column: String,
    },
}

#[derive(Clone, Copy)]
pub struct TailGene {
    pub family: FamilyId,
    pub strand: Strand,
}
impl std::fmt::Debug for TailGene {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}/{}", self.family, self.strand)
    }
}
impl std::cmp::PartialEq for TailGene {
    // NOTE: strand is deliberately ignored
    fn eq(&self, other: &Self) -> bool {
        self.family == other.family
    }
}
impl std::cmp::Eq for TailGene {}

#[derive(Clone, Default)]
pub struct Gene {
    pub id: String,
    pub species: String,
    pub family: FamilyId,
    pub chr: String,
    pub pos: usize,
    pub strand: Strand,
    pub left_landscape: Vec<TailGene>,
    pub right_landscape: Vec<TailGene>,
}
impl Gene {
    pub fn landscape(&self) -> impl Iterator<Item = TailGene> + '_ {
        self.left_landscape
            .iter()
            .cloned()
            .chain(std::iter::once(TailGene {
                family: self.family,
                strand: self.strand,
            }))
            .chain(self.right_landscape.iter().cloned())
    }
}

impl GeneBook {
    fn parse_landscape(landscape: &str) -> Result<Vec<TailGene>> {
        fn parse_tailgene(g: &str) -> Result<TailGene> {
            let strand = g
                .chars()
                .next()
                .and_then(|c| c.try_into().ok())
                .unwrap_or_default();
            let family_id = g
                .strip_prefix(['+', '-', '.'])
                .unwrap_or(g)
                .parse::<usize>()
                .with_context(|| format!("invalid family ID in landscape entry: {g:?}"))?;
            Ok(TailGene {
                family: family_id.into(),
                strand,
            })
        }

        if landscape.is_empty() {
            Ok(Vec::new())
        } else {
            landscape
                .split(LANDSCAPE_DELIMITER)
                .map(parse_tailgene)
                .collect()
        }
    }

    fn get_rows<P: rusqlite::Params>(
        mut query: rusqlite::Statement,
        params: P,
        window: usize,
    ) -> Result<HashMap<String, Gene>> {
        let genes = query
            .query_map(params, |r| {
                std::result::Result::Ok((
                    r.get::<_, String>(0)?, // id
                    r.get::<_, String>(1)?, // left tail
                    r.get::<_, String>(2)?, // right tail
                    r.get::<_, i64>(3)?
                        .try_into()
                        .expect("SQLite integer should fit in usize"), // ancestral id
                    r.get::<_, String>(4)?, // species
                    r.get::<_, String>(5)?, // chr
                    r.get::<_, i64>(6)?
                        .try_into()
                        .expect("SQLite integer should fit in usize"), // position
                    r.get::<_, String>(7)?, // direction
                ))
            })?
            .collect::<Result<Vec<_>, _>>()?;

        genes
            .into_iter()
            .map(|g| {
                let id = g.0.to_string();
                let mut left_landscape = Self::parse_landscape(&g.1)?;
                left_landscape.reverse();
                left_landscape.truncate(window);
                left_landscape.reverse();

                let mut right_landscape = Self::parse_landscape(&g.2)?;
                right_landscape.truncate(window);

                let strand =
                    g.7.as_str()
                        .try_into()
                        .with_context(|| format!("invalid strand {:?} in database", g.7))?;

                Ok((
                    g.0.clone(),
                    Gene {
                        id,
                        species: g.4,
                        family: g.3,
                        chr: g.5,
                        pos: g.6,
                        strand,
                        left_landscape,
                        right_landscape,
                    },
                ))
            })
            .collect()
    }

    pub fn in_memory(filename: &str, window: usize, id_column: &str) -> Result<Self> {
        info!("Caching the database...");

        let conn = Connection::open(filename).map_err(|e| errors::DataError::FailedToConnect {
            source: e,
            filename: filename.into(),
        })?;
        let query = conn.prepare(&format!(
            "SELECT {id_column}, left_tail_ids, right_tail_ids, ancestral_id, species, chr, start, direction FROM genomes"
        ))?;
        let genes = Self::get_rows(query, [], window)?;
        let species = conn
            .prepare("SELECT DISTINCT species FROM genomes")?
            .query_map([], |row| row.get::<_, String>(0))?
            .collect::<Result<Vec<_>, _>>()?;

        info!("Done.");
        Ok(GeneBook::InMemory { genes, species })
    }

    pub fn cached<S: AsRef<str>>(
        filename: &str,
        window: usize,
        id_column: &str,
        ids: &[S],
    ) -> Result<Self> {
        info!("Caching the database...");

        let conn = Connection::open(filename).map_err(|e| errors::DataError::FailedToConnect {
            source: e,
            filename: filename.into(),
        })?;

        let query = conn.prepare(&format!(
            "SELECT {id_column}, left_tail_ids, right_tail_ids, ancestral_id, species, chr, start, direction FROM genomes WHERE {id_column} IN ({})",
            std::iter::repeat_n("?", ids.len()).collect::<Vec<_>>().join(", ")
        ))?;
        let genes = Self::get_rows(
            query,
            rusqlite::params_from_iter(ids.iter().map(|s| s.as_ref())),
            window,
        )?;
        let species = conn
            .prepare("SELECT DISTINCT species FROM genomes")?
            .query_map([], |row| row.get::<_, String>(0))?
            .collect::<Result<Vec<_>, _>>()?;

        Ok(GeneBook::InMemory { genes, species })
    }

    #[allow(dead_code)]
    pub fn inline(filename: &str, window: usize, id_column: &str) -> Result<Self> {
        let conn = Connection::open(filename).map_err(|e| errors::DataError::FailedToConnect {
            source: e,
            filename: filename.into(),
        })?;
        Ok(GeneBook::Inline {
            conn: Mutex::new(conn),
            window,
            id_column: id_column.to_owned(),
        })
    }

    pub fn get(&self, g: &str) -> Result<Gene> {
        match self {
            GeneBook::InMemory { genes, .. } => genes
                .get(g)
                .cloned()
                .ok_or_else(|| errors::DataError::UnknownId(g.to_owned()).into()),
            GeneBook::Inline {
                conn: conn_mutex,
                window,
                id_column,
            } => {
                let conn = conn_mutex.lock().expect("MUTEX POISONING");
                let mut query = conn.prepare(
                    &format!("SELECT left_tail_ids, right_tail_ids, ancestral_id, species, chr, start, direction FROM genomes WHERE {id_column}=?"),
                )?;
                let (left_str, right_str, ancestral_id, species, chr, pos, strand_str) = query
                    .query_row([g], |r| {
                        rusqlite::Result::Ok((
                            r.get::<_, String>(0)?,
                            r.get::<_, String>(1)?,
                            r.get::<_, i64>(2)?,
                            r.get::<_, String>(3)?,
                            r.get::<_, String>(4)?,
                            r.get::<_, i64>(5)?,
                            r.get::<_, String>(6)?,
                        ))
                    })
                    .with_context(|| "while accessing DB")?;

                let mut left_landscape = Self::parse_landscape(&left_str)?;
                left_landscape.reverse();
                left_landscape.truncate(*window);
                left_landscape.reverse();

                let mut right_landscape = Self::parse_landscape(&right_str)?;
                right_landscape.truncate(*window);

                let strand = strand_str
                    .as_str()
                    .try_into()
                    .with_context(|| format!("invalid strand {strand_str:?} in database"))?;

                Ok(Gene {
                    id: g.to_string(),
                    species,
                    family: ancestral_id
                        .try_into()
                        .expect("SQLite integer should fit in usize"),
                    chr,
                    pos: pos.try_into().expect("SQLite integer should fit in usize"),
                    strand,
                    left_landscape,
                    right_landscape,
                })
            }
        }
    }

    pub fn get_mut(&mut self, g: &str) -> Result<&mut Gene> {
        match self {
            GeneBook::InMemory { genes, .. } => genes
                .get_mut(g)
                .ok_or_else(|| errors::DataError::UnknownId(g.to_owned()).into()),
            GeneBook::Inline { .. } => Err(errors::DataError::ImmutableBook.into()),
        }
    }

    pub fn species(&self) -> Result<Vec<String>> {
        match self {
            GeneBook::InMemory { species, .. } => Ok(species.to_owned()),
            GeneBook::Inline {
                conn: conn_mutex, ..
            } => {
                let conn = conn_mutex.lock().expect("MUTEX POISONING");
                let result = conn
                    .prepare("SELECT DISTINCT species FROM genomes")?
                    .query_map([], |row| row.get::<_, String>(0))?
                    .collect::<std::result::Result<Vec<_>, _>>()
                    .map_err(Into::into);
                result
            }
        }
    }
}
