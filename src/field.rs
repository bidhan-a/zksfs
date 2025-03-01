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
                "Modulus cannot be zero.".into(),
            ));
        }
        Ok(FieldElement { value, modulus })
    }

    /// Add two field elements.
    pub fn add(&self, other: &FieldElement) -> Result<Self, ZKError> {
        if self.modulus != other.modulus {
            return Err(ZKError::InvalidFieldElement(
                "Moduli must be the same for addition.".into(),
            ));
        }

        FieldElement::new((self.value + other.value) % self.modulus, self.modulus)
    }

    /// Subtract two field elements.
    pub fn sub(&self, other: &FieldElement) -> Result<Self, ZKError> {
        if self.modulus != other.modulus {
            return Err(ZKError::InvalidFieldElement(
                "Moduli must be the same for subtraction.".into(),
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
                "Moduli must be the same for multiplication.".into(),
            ));
        }
        FieldElement::new((self.value * other.value) % self.modulus, self.modulus)
    }

    /// Find the modular inverse of the field element.
    pub fn inv(&self) -> Result<FieldElement, ZKError> {
        let v = self.value as i64;
        let m = self.modulus as i64;

        let (g, x, _) = Self::eegcd(v, m);
        if g != 1 {
            return Err(ZKError::InvalidFieldElement(
                "Modular inverse does not exist.".into(),
            ));
        }

        // Make sure the inverse is positive.
        let x_pos = ((x % m) + m) % m;
        FieldElement::new(x_pos as u64, self.modulus)
    }

    /// Exponentiate the field element by the provided exponent.
    pub fn exp(&self, exponent: u64) -> Result<FieldElement, ZKError> {
        let mut result = FieldElement::new(1, self.modulus)?;
        let mut base = self.clone();
        let mut exp = exponent;

        while exp > 0 {
            if exp % 2 == 1 {
                result = result.mul(&base)?;
            }
            base = base.mul(&base)?;
            exp /= 2;
        }

        Ok(result)
    }

    fn eegcd(a: i64, b: i64) -> (i64, i64, i64) {
        if a == 0 {
            (b, 0, 1)
        } else {
            let (g, x, y) = Self::eegcd(b % a, a);
            (g, y - (b / a) * x, x)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add() {
        let a = FieldElement::new(3, 7).unwrap();
        let b = FieldElement::new(4, 7).unwrap();
        let result = a.add(&b).unwrap();
        assert_eq!(result.value, 0);
    }

    #[test]
    fn test_sub() {
        let a = FieldElement::new(3, 7).unwrap();
        let b = FieldElement::new(4, 7).unwrap();
        let result = a.sub(&b).unwrap();
        assert_eq!(result.value, 6);
    }

    #[test]
    fn test_mul() {
        let a = FieldElement::new(3, 7).unwrap();
        let b = FieldElement::new(4, 7).unwrap();
        let result = a.mul(&b).unwrap();
        assert_eq!(result.value, 5);
    }

    #[test]
    fn test_inv() {
        let a = FieldElement::new(3, 7).unwrap();
        let a_inv = a.inv().unwrap();
        // a x a_inv should equal to 1 mod 7.
        let one = a.mul(&a_inv).unwrap();
        assert_eq!(one.value, 1);
    }

    #[test]
    fn test_exp() {
        let a = FieldElement::new(3, 7).unwrap();
        let a_exp = a.exp(3).unwrap();
        assert_eq!(a_exp.value, 6);
    }
}
