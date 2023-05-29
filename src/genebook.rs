use anyhow::*;
use log::*;
use rusqlite::Connection;
use std::collections::HashMap;
use std::sync::Mutex;

use crate::errors;

pub type FamilyID = usize;

#[allow(dead_code)]
pub enum GeneBook {
    InMemory(HashMap<String, Gene>),
    Cached(HashMap<String, Gene>),
    Inline(Mutex<Connection>, usize, String),
}

#[derive(Clone, Default)]
pub struct Gene {
    pub id: String,
    pub species: String,
    pub family: FamilyID,
    pub chr: String,
    pub pos: usize,
    pub left_landscape: Vec<FamilyID>,
    pub right_landscape: Vec<FamilyID>,
}
impl Gene {
    pub fn landscape(&self) -> impl Iterator<Item = FamilyID> + '_ {
        self.left_landscape
            .iter()
            .cloned()
            .chain(vec![self.family].into_iter())
            .chain(self.right_landscape.iter().cloned())
    }
}

impl GeneBook {
    fn parse_landscape(landscape: &str) -> Vec<usize> {
        if landscape.is_empty() {
            Vec::new()
        } else {
            landscape
                .split('.')
                .map(|x| x.parse::<usize>().unwrap())
                .collect::<Vec<_>>()
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
                ))
            })?
            .collect::<Result<Vec<_>, _>>()?;

        Ok(genes
            .into_iter()
            .map(|g| {
                let mut left_landscape = Self::parse_landscape(&g.1);
                left_landscape.reverse();
                left_landscape.truncate(window);
                left_landscape.reverse();

                let mut right_landscape = Self::parse_landscape(&g.2);
                right_landscape.truncate(window);

                (
                    g.0.clone(),
                    Gene {
                        id: g.0.to_string(),
                        species: g.4,
                        family: g.3,
                        chr: g.5,
                        pos: g.6,
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
            "SELECT {id_column}, left_tail_ids, right_tail_ids, ancestral_id, species, chr, start FROM genomes"
        ))?;
        let r = Self::get_rows(query, [], window)?;
        info!("Done.");
        Ok(GeneBook::InMemory(r))
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
            "SELECT {id_column}, left_tail_ids, right_tail_ids, ancestral_id, species, chr, start FROM genomes WHERE {id_column} IN ({})",
            std::iter::repeat("?").take(ids.len()).collect::<Vec<_>>().join(", ")
        ))?;
        let r = Self::get_rows(
            query,
            rusqlite::params_from_iter(ids.iter().map(|s| s.as_ref())),
            window,
        )?;
        Ok(GeneBook::Cached(r))
    }

    #[allow(dead_code)]
    pub fn inline(filename: &str, window: usize, id_column: &str) -> Result<Self> {
        let conn = Connection::open(filename).map_err(|e| errors::DataError::FailedToConnect {
            source: e,
            filename: filename.into(),
        })?;
        Ok(GeneBook::Inline(
            Mutex::new(conn),
            window,
            id_column.to_owned(),
        ))
    }

    pub fn get(&self, g: &str) -> Result<Gene> {
        match self {
            GeneBook::InMemory(book) | GeneBook::Cached(book) => book
                .get(g)
                .cloned()
                .ok_or_else(|| errors::DataError::UnknownId(g.to_owned()).into()),
            GeneBook::Inline(conn_mutex, window, id_column) => {
                let conn = conn_mutex.lock().expect("MUTEX POISONING");
                let mut query = conn.prepare(
                    &format!("SELECT left_tail_ids, right_tail_ids, ancestral_id, species, chr, start FROM genomes WHERE {id_column}=?"),
                )?;
                query
                    .query_row(&[g], |r| {
                        let species = r.get::<_, String>(3)?;

                        let mut left_landscape = Self::parse_landscape(&r.get::<_, String>(0)?);
                        left_landscape.reverse();
                        left_landscape.truncate(*window);
                        left_landscape.reverse();

                        let mut right_landscape = Self::parse_landscape(&r.get::<_, String>(1)?);
                        right_landscape.truncate(*window);

                        rusqlite::Result::Ok(Gene {
                            id: g.to_string(),
                            species,
                            family: r.get::<usize, _>(2)?,
                            chr: r.get::<_, String>(4)?,
                            pos: r.get::<usize, _>(5)?,
                            left_landscape,
                            right_landscape,
                        })
                    })
                    .with_context(|| "while accessing DB")
            }
        }
    }
}
