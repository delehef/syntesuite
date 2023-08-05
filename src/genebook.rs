use anyhow::*;
use log::*;
use rusqlite::Connection;
use std::collections::HashMap;
use std::sync::Mutex;

use crate::{errors, Strand};

pub type FamilyID = usize;

#[allow(dead_code)]
pub enum GeneBook {
    InMemory {
        genes: HashMap<String, Gene>,
        species: Vec<String>,
    },
    Cached {
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
    pub family: FamilyID,
    pub strand: Strand,
}
impl std::fmt::Debug for TailGene {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}/{}", self.family, self.strand)
    }
}
impl std::cmp::PartialEq for TailGene {
    fn eq(&self, other: &Self) -> bool {
        self.family == other.family
    }
}
impl std::cmp::Eq for TailGene {}

#[derive(Clone, Default)]
pub struct Gene {
    pub id: String,
    pub species: String,
    pub family: FamilyID,
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
    fn parse_landscape(landscape: &str) -> Vec<TailGene> {
        fn parse_tailgene(g: &str) -> TailGene {
            let strand = g
                .chars()
                .next()
                .and_then(|c| c.try_into().ok())
                .unwrap_or_default();
            let family_id = g
                .strip_prefix(['+', '-', '.'])
                .unwrap_or(g)
                .parse::<usize>()
                .unwrap();
            TailGene {
                family: family_id,
                strand,
            }
        }

        if landscape.is_empty() {
            Vec::new()
        } else {
            landscape.split('.').map(parse_tailgene).collect::<Vec<_>>()
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
                    r.get::<_, usize>(3)?,  // ancestral id
                    r.get::<_, String>(4)?, // species
                    r.get::<_, String>(5)?, // chr
                    r.get::<_, usize>(6)?,  // position
                    r.get::<_, String>(7)?, // direction
                ))
            })?
            .collect::<Result<Vec<_>, _>>()?;

        Ok(genes
            .into_iter()
            .map(|g| {
                let id = g.0.to_string();
                let mut left_landscape = Self::parse_landscape(&g.1);
                left_landscape.reverse();
                left_landscape.truncate(window);
                left_landscape.reverse();

                let mut right_landscape = Self::parse_landscape(&g.2);
                right_landscape.truncate(window);

                (
                    g.0.clone(),
                    Gene {
                        id,
                        species: g.4,
                        family: g.3,
                        chr: g.5,
                        pos: g.6,
                        strand: g.7.as_str().try_into().unwrap(),
                        left_landscape,
                        right_landscape,
                    },
                )
            })
            .collect())
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
            std::iter::repeat("?").take(ids.len()).collect::<Vec<_>>().join(", ")
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

        Ok(GeneBook::Cached { genes, species })
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
            GeneBook::InMemory { genes, .. } | GeneBook::Cached { genes, .. } => genes
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
                query
                    .query_row([g], |r| {
                        let species = r.get::<_, String>(3)?;

                        let mut left_landscape = Self::parse_landscape(&r.get::<_, String>(0)?);
                        left_landscape.reverse();
                        left_landscape.truncate(*window);
                        left_landscape.reverse();

                        let mut right_landscape = Self::parse_landscape(&r.get::<_, String>(1)?);
                        right_landscape.truncate(*window);

                        let strand = r
                            .get::<_, String>(6)?
                            .chars()
                            .next()
                            .and_then(|c| c.try_into().ok())
                            .unwrap_or_default();

                        rusqlite::Result::Ok(Gene {
                            id: g.to_string(),
                            species,
                            family: r.get::<usize, _>(2)?,
                            chr: r.get::<_, String>(4)?,
                            pos: r.get::<usize, _>(5)?,
                            strand,
                            left_landscape,
                            right_landscape,
                        })
                    })
                    .with_context(|| "while accessing DB")
            }
        }
    }

    pub fn get_mut(&mut self, g: &str) -> Result<&mut Gene> {
        match self {
            GeneBook::InMemory { genes, .. } | GeneBook::Cached { genes, .. } => genes
                .get_mut(g)
                .ok_or_else(|| errors::DataError::UnknownId(g.to_owned()).into()),
            GeneBook::Inline { .. } => Err(errors::DataError::ImmutableBook.into()),
        }
    }

    pub fn species(&self) -> Vec<String> {
        match self {
            GeneBook::InMemory { species, .. } | GeneBook::Cached { species, .. } => {
                species.to_owned()
            }
            GeneBook::Inline {
                conn: conn_mutex, ..
            } => {
                let conn = conn_mutex.lock().expect("MUTEX POISONING");
                let species = conn
                    .prepare("SELECT DISTINCT species FROM genomes")
                    .unwrap()
                    .query_map([], |row| row.get::<_, String>(0))
                    .unwrap()
                    .collect::<Result<Vec<_>, _>>()
                    .unwrap();
                species
            }
        }
    }
}
