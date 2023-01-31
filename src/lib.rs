use thiserror::Error;

pub mod dbmaker;
mod errors;
pub mod genebook;
mod gff;

#[derive(Debug, Error)]
pub enum DataError {
    #[error("invalid phase value: {0}")]
    InvalidPhase(String),

    #[error("invalid strand value: {0}")]
    InvalidStrand(String),
}

#[derive(Debug, Copy, Clone)]
pub enum Phase {
    Sync,
    OneShifted,
    TwoShifted,
}
impl TryFrom<&str> for Phase {
    type Error = DataError;

    fn try_from(x: &str) -> Result<Self, Self::Error> {
        match x {
            "0" => Ok(Phase::Sync),
            "1" => Ok(Phase::OneShifted),
            "2" => Ok(Phase::TwoShifted),
            _ => Err(DataError::InvalidPhase(x.to_string())),
        }
    }
}
impl TryFrom<usize> for Phase {
    type Error = DataError;

    fn try_from(x: usize) -> Result<Self, Self::Error> {
        match x {
            0 => Ok(Phase::Sync),
            1 => Ok(Phase::OneShifted),
            2 => Ok(Phase::TwoShifted),
            _ => Err(DataError::InvalidPhase(x.to_string())),
        }
    }
}
impl From<Phase> for usize {
    fn from(p: Phase) -> Self {
        match p {
            Phase::Sync => 0,
            Phase::OneShifted => 1,
            Phase::TwoShifted => 2,
        }
    }
}

#[derive(Debug, Copy, Clone)]
pub enum Strand {
    Direct,
    Reverse,
}
impl TryFrom<&str> for Strand {
    type Error = DataError;

    fn try_from(s: &str) -> Result<Self, Self::Error> {
        match s {
            "+" => Ok(Strand::Direct),
            "-" => Ok(Strand::Reverse),
            _ => Err(DataError::InvalidStrand(s.to_string())),
        }
    }
}
impl From<Strand> for char {
    fn from(s: Strand) -> Self {
        match s {
            Strand::Direct => '+',
            Strand::Reverse => '-',
        }
    }
}
impl From<Strand> for String {
    fn from(s: Strand) -> Self {
        match s {
            Strand::Direct => "+".into(),
            Strand::Reverse => "-".into(),
        }
    }
}
