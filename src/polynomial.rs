use crate::{errors::ZKError, field::FieldElement};

/// Represents a polynomial with coefficients in a finite field.
#[derive(Clone, Debug)]
pub struct Polynomial {
    pub coefficients: Vec<FieldElement>,
}

impl Polynomial {
    /// Creates a new polynomial.
    pub fn new(mut coefficients: Vec<FieldElement>) -> Result<Self, ZKError> {
        if coefficients.is_empty() {
            return Err(ZKError::CircuitError(
                "Polynomial must have at least one coefficient".to_string(),
            ));
        }

        // Verify that all coefficients have the same modulus.
        let modulus = coefficients[0].modulus;
        for coeff in &coefficients {
            if coeff.modulus != modulus {
                return Err(ZKError::PolynomialError(
                    "All coefficients must have the same modulus".to_string(),
                ));
            }
        }

        Ok(Self { coefficients })
    }

    /// Returns the degree of the polynomial.
    pub fn degree(&self) -> usize {
        let mut deg = self.coefficients.len() - 1;
        while deg > 0 && self.coefficients[deg].value == 0 {
            deg -= 1;
        }
        deg
    }

    /// Evaluates the polynomial at the given field element.
    pub fn evaluate(&self, fe: &FieldElement) -> Result<FieldElement, ZKError> {
        if self.coefficients[0].modulus != fe.modulus {
            return Err(ZKError::PolynomialError(
                "Moduli must be the same for evaluation".to_string(),
            ));
        }

        let mut result = FieldElement::new(0, fe.modulus)?;

        // Evaluate from highest degree coefficient downwards.
        for coeff in self.coefficients.iter().rev() {
            result = result.mul(fe)?;
            result = result.add(coeff)?;
        }

        Ok(result)
    }

    /// Adds two polynomials.
    pub fn add(&self, other: &Polynomial) -> Result<Polynomial, ZKError> {
        todo!()
    }

    /// Subtracts two polynomials.
    pub fn sub(&self, other: &Polynomial) -> Result<Polynomial, ZKError> {
        todo!()
    }

    /// Multiplies two polynomials.
    pub fn mul(&self, other: &Polynomial) -> Result<Polynomial, ZKError> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use crate::field::FieldElement;

    use super::Polynomial;

    #[test]
    fn test_evaluate() {
        // Define polynomial: 2 + 4x + 6x^2 mod 97.
        let modulus = 97;
        let coefficients = vec![
            FieldElement::new(2, modulus).unwrap(),
            FieldElement::new(4, modulus).unwrap(),
            FieldElement::new(6, modulus).unwrap(),
        ];
        let polynomial = Polynomial::new(coefficients).unwrap();

        // Evaluate at x = 3.
        let x = FieldElement::new(3, modulus).unwrap();
        let result = polynomial.evaluate(&x).unwrap();
        assert_eq!(result, FieldElement::new(68, modulus).unwrap());
    }
}
