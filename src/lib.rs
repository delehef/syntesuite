use errors::ParseError;

mod bed;
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

#[derive(Debug, Copy, Clone)]
pub enum Strand {
    Direct,
    Reverse,
}
impl TryFrom<&str> for Strand {
    type Error = ParseError;

    fn try_from(s: &str) -> Result<Self, Self::Error> {
        match s {
            "+" => Ok(Strand::Direct),
            "-" => Ok(Strand::Reverse),
            _ => Err(ParseError::InvalidStrand(s.to_string())),
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

enum Record {
    Gff(gff::GffRecord),
    Bed(bed::BedRecord),
}

impl Record {
    fn id(&self) -> Option<&str> {
        match self {
            Record::Bed(r) => r.id(),
            Record::Gff(r) => r.id(),
        }
    }
    fn chr(&self) -> &str {
        match self {
            Record::Gff(r) => r.chr(),
            Record::Bed(r) => r.chr(),
        }
    }
    fn start(&self) -> usize {
        match self {
            Record::Gff(r) => r.start(),
            Record::Bed(r) => r.start(),
        }
    }
    fn end(&self) -> usize {
        match self {
            Record::Gff(r) => r.end(),
            Record::Bed(r) => r.end(),
        }
    }
    fn strand(&self) -> Strand {
        match self {
            Record::Gff(r) => r.strand().unwrap_or(Strand::Direct),
            Record::Bed(r) => r.strand(),
        }
    }
    fn is_class(&self, class: &str) -> bool {
        match self {
            Record::Bed(_) => true,
            Record::Gff(r) => r.class().map(|c| c == class).unwrap_or(false),
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
