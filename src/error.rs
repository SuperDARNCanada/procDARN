use std::fmt;
use std::fmt::Formatter;

#[derive(Debug)]
pub struct BackscatterError {
    pub details: String,
}

impl fmt::Display for BackscatterError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.details)
    }
}

impl BackscatterError {
    pub fn new(details: &str) -> BackscatterError {
        BackscatterError {
            details: details.to_string(),
        }
    }
}
