use thiserror::Error;

#[derive(Debug, Error)]
pub enum ZKError {
    #[error("Invalid field element: {0}")]
    InvalidFieldElement(String),
}
