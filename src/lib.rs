use std::write;

use errors::ParseError;

mod bed;
mod chrom;
pub mod dbmaker;
mod errors;
pub mod genebook;
mod gff;

#[derive(Debug, Copy, Clone)]
pub enum Phase {
    Sync,
    OneShifted,
    TwoShifted,
}
impl TryFrom<&str> for Phase {
    type Error = ParseError;

    fn try_from(x: &str) -> Result<Self, Self::Error> {
        match x {
            "0" => Ok(Phase::Sync),
            "1" => Ok(Phase::OneShifted),
            "2" => Ok(Phase::TwoShifted),
            _ => Err(ParseError::InvalidPhase(x.to_string())),
        }
    }
}
impl TryFrom<usize> for Phase {
    type Error = ParseError;

    fn try_from(x: usize) -> Result<Self, Self::Error> {
        match x {
            0 => Ok(Phase::Sync),
            1 => Ok(Phase::OneShifted),
            2 => Ok(Phase::TwoShifted),
            _ => Err(ParseError::InvalidPhase(x.to_string())),
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

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum Strand {
    Direct,
    Reverse,
    Unknown,
}
impl Strand {
    pub fn reverse(&mut self) {
        match self {
            Strand::Direct => *self = Strand::Reverse,
            Strand::Reverse => *self = Strand::Direct,
            Strand::Unknown => {}
        }
    }
}
impl std::default::Default for Strand {
    fn default() -> Strand {
        Strand::Unknown
    }
}
impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", Into::<char>::into(*self))
    }
}
impl TryFrom<&str> for Strand {
    type Error = ParseError;

    fn try_from(s: &str) -> Result<Self, Self::Error> {
        match s {
            "+" | "1" | "+1" => Ok(Strand::Direct),
            "-" | "-1" => Ok(Strand::Reverse),
            "." => Ok(Strand::Unknown),
            _ => Err(ParseError::InvalidStrand(s.to_string())),
        }
    }
}
impl TryFrom<char> for Strand {
    type Error = ParseError;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c {
            '+' => Ok(Strand::Direct),
            '-' => Ok(Strand::Reverse),
            '.' => Ok(Strand::Unknown),
            _ => Err(ParseError::InvalidStrand(c.to_string())),
        }
    }
}
impl From<Strand> for char {
    fn from(s: Strand) -> Self {
        match s {
            Strand::Direct => '+',
            Strand::Reverse => '-',
            Strand::Unknown => '.',
        }
    }
}
impl From<Strand> for String {
    fn from(s: Strand) -> Self {
        match s {
            Strand::Direct => "+".into(),
            Strand::Reverse => "-".into(),
            Strand::Unknown => "-".into(),
        }
    }
}

enum Record {
    Gff(gff::GffRecord),
    Bed(bed::BedRecord),
    Chrom(chrom::ChromRecord),
}

impl Record {
    fn id(&self) -> Option<&str> {
        match self {
            Record::Bed(r) => r.id(),
            Record::Gff(r) => r.id(),
            Record::Chrom(r) => Some(r.id()),
        }
    }
    fn chr(&self) -> &str {
        match self {
            Record::Gff(r) => r.chr(),
            Record::Bed(r) => r.chr(),
            Record::Chrom(r) => r.chr(),
        }
    }
    fn start(&self) -> usize {
        match self {
            Record::Gff(r) => r.start(),
            Record::Bed(r) => r.start(),
            Record::Chrom(r) => r.start(),
        }
    }
    fn end(&self) -> usize {
        match self {
            Record::Gff(r) => r.end(),
            Record::Bed(r) => r.end(),
            Record::Chrom(r) => r.end(),
        }
    }
    fn strand(&self) -> Strand {
        match self {
            Record::Gff(r) => r.strand().unwrap_or(Strand::Direct),
            Record::Bed(r) => r.strand(),
            Record::Chrom(r) => r.strand(),
        }
    }
    fn is_class(&self, class: &str) -> bool {
        match self {
            Record::Gff(r) => r.class().map(|c| c == class).unwrap_or(false),
            Record::Bed(_) => true,
            Record::Chrom(_) => true,
        }
    }
}

impl From<gff::GffRecord> for Record {
    fn from(r: gff::GffRecord) -> Self {
        Record::Gff(r)
    }
}
impl From<bed::BedRecord> for Record {
    fn from(r: bed::BedRecord) -> Self {
        Record::Bed(r)
    }
}
impl From<chrom::ChromRecord> for Record {
    fn from(r: chrom::ChromRecord) -> Self {
        Record::Chrom(r)
    }
}
