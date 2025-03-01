use crate::errors::ZKError;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FieldElement {
    pub value: u64,
    pub modulus: u64,
}

impl FieldElement {
    /// Create a new field element with value and modulus.
    pub fn new(value: u64, modulus: u64) -> Result<Self, ZKError> {
        if modulus == 0 {
            return Err(ZKError::InvalidFieldElement(
                "Modulus cannot be zero".into(),
            ));
        }
        Ok(FieldElement { value, modulus })
    }

    /// Add two field elements.
    pub fn add(&self, other: &FieldElement) -> Result<Self, ZKError> {
        if self.modulus != other.modulus {
            return Err(ZKError::InvalidFieldElement(
                "Moduli must be the same for addition".into(),
            ));
        }

        FieldElement::new((self.value + other.value) % self.modulus, self.modulus)
    }

    /// Subtract two field elements.
    pub fn sub(&self, other: &FieldElement) -> Result<Self, ZKError> {
        if self.modulus != other.modulus {
            return Err(ZKError::InvalidFieldElement(
                "Moduli must be the same for subtraction".into(),
            ));
        }

        // Ensure non-negative result by adding the modulus before subtracting.
        let diff = (self.value + self.modulus - other.value) % self.modulus;
        FieldElement::new(diff, self.modulus)
    }

    /// Multiply two field elements.
    pub fn mul(&self, other: &FieldElement) -> Result<FieldElement, ZKError> {
        if self.modulus != other.modulus {
            return Err(ZKError::InvalidFieldElement(
                "Moduli must be the same for multiplication".into(),
            ));
        }
        FieldElement::new((self.value * other.value) % self.modulus, self.modulus)
    }

    /// Find the modular inverse of the field element.
    pub fn inv(&self) -> Result<FieldElement, ZKError> {
        todo!()
    }

    /// Exponentiate the field element by the provided exponent.
    pub fn exp(&self, exponent: u64) -> Result<FieldElement, ZKError> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_addition() {
        let a = FieldElement::new(3, 7).unwrap();
        let b = FieldElement::new(4, 7).unwrap();
        let result = a.add(&b).unwrap();
        assert_eq!(result.value, 0);
    }

    #[test]
    fn test_subtraction() {
        let a = FieldElement::new(3, 7).unwrap();
        let b = FieldElement::new(4, 7).unwrap();
        let result = a.sub(&b).unwrap();
        assert_eq!(result.value, 6);
    }

    #[test]
    fn test_multiplication() {
        let a = FieldElement::new(3, 7).unwrap();
        let b = FieldElement::new(4, 7).unwrap();
        let result = a.mul(&b).unwrap();
        assert_eq!(result.value, 5);
    }
}
