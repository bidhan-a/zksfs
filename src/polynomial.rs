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

        // Evaluate from highest degree coefficient downwards (Horner's method).
        for coeff in self.coefficients.iter().rev() {
            result = result.mul(fe)?;
            result = result.add(coeff)?;
        }

        Ok(result)
    }

    /// Adds two polynomials.
    pub fn add(&self, other: &Polynomial) -> Result<Polynomial, ZKError> {
        if self.coefficients[0].modulus != other.coefficients[0].modulus {
            return Err(ZKError::PolynomialError(
                "Moduli must be the same for addition".to_string(),
            ));
        }

        let max_len = self.coefficients.len().max(other.coefficients.len());
        let modulus = self.coefficients[0].modulus;
        let mut sum = Vec::new();

        for i in 0..max_len {
            let a = self
                .coefficients
                .get(i)
                .cloned()
                .unwrap_or(FieldElement::new(0, modulus)?);
            let b = other
                .coefficients
                .get(i)
                .cloned()
                .unwrap_or(FieldElement::new(0, modulus)?);
            sum.push(a.add(&b)?);
        }

        Polynomial::new(sum)
    }

    /// Subtracts two polynomials.
    pub fn sub(&self, other: &Polynomial) -> Result<Polynomial, ZKError> {
        if self.coefficients[0].modulus != other.coefficients[0].modulus {
            return Err(ZKError::PolynomialError(
                "Moduli must be the same for subtraction".to_string(),
            ));
        }

        let max_len = self.coefficients.len().max(other.coefficients.len());
        let modulus = self.coefficients[0].modulus;
        let mut diff = Vec::new();

        for i in 0..max_len {
            let a = self
                .coefficients
                .get(i)
                .cloned()
                .unwrap_or(FieldElement::new(0, modulus)?);
            let b = other
                .coefficients
                .get(i)
                .cloned()
                .unwrap_or(FieldElement::new(0, modulus)?);
            diff.push(a.sub(&b)?);
        }

        Polynomial::new(diff)
    }

    /// Multiplies two polynomials.
    pub fn mul(&self, other: &Polynomial) -> Result<Polynomial, ZKError> {
        if self.coefficients[0].modulus != other.coefficients[0].modulus {
            return Err(ZKError::PolynomialError(
                "Moduli must be the same for multiplication".to_string(),
            ));
        }

        let n = self.coefficients.len();
        let m = other.coefficients.len();
        let modulus = self.coefficients[0].modulus;
        let mut product = vec![FieldElement::new(0, modulus)?; n + m - 1];

        for i in 0..n {
            for j in 0..m {
                let prod = self.coefficients[i].mul(&other.coefficients[j])?;
                product[i + j] = product[i + j].add(&prod)?;
            }
        }
        Polynomial::new(product)
    }

    /// Performs polynomial long division and returns the quotient and the remainder.
    pub fn div(&self, other: &Polynomial) -> Result<(Polynomial, Polynomial), ZKError> {
        let modulus = self.coefficients[0].modulus;
        if modulus != other.coefficients[0].modulus {
            return Err(ZKError::PolynomialError(
                "Moduli must be the same for division".to_string(),
            ));
        }

        let mut remainder = self.clone();
        let quotient_size = self.degree().saturating_sub(other.degree()) + 1;
        let mut quotient_coefficients = vec![FieldElement::new(0, modulus)?; quotient_size];

        while remainder.degree() >= other.degree()
            && remainder.coefficients.len() > 0
            && remainder.coefficients[remainder.degree()].value != 0
        {
            let deg_diff = remainder.degree() - other.degree();
            let lead_dividend = remainder.coefficients[remainder.degree()].clone();
            let lead_divisor = other.coefficients[other.degree()].clone();
            let factor = lead_dividend.mul(&lead_divisor.inv()?)?;
            // Create a polynomial factor_poly = factor * x^(deg_diff)
            let mut factor_poly_coefficients = vec![FieldElement::new(0, modulus)?; deg_diff];
            factor_poly_coefficients.push(factor.clone());
            let factor_poly = Polynomial::new(factor_poly_coefficients)?;

            quotient_coefficients[deg_diff] = quotient_coefficients[deg_diff].add(&factor)?;
            let subtrahend = factor_poly.mul(other)?;
            remainder = remainder.sub(&subtrahend)?;
        }

        let quotient = Polynomial::new(quotient_coefficients)?;
        Ok((quotient, remainder))
    }

    /// Scales the polynomial by a scalar field element.
    pub fn scale(&self, scalar: &FieldElement) -> Result<Polynomial, ZKError> {
        let scaled_coefficients = self
            .coefficients
            .iter()
            .map(|c| c.mul(scalar).unwrap())
            .collect();
        Polynomial::new(scaled_coefficients)
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

    #[test]
    fn test_add() {
        let modulus = 97;
        // Define polynomial: 1 + 2x mod 97.
        let coefficients1 = vec![
            FieldElement::new(1, modulus).unwrap(),
            FieldElement::new(2, modulus).unwrap(),
        ];
        let polynomial1 = Polynomial::new(coefficients1).unwrap();

        // Define polynomial: 2 + 3x + 4x^2 mod 97.
        let coefficients2 = vec![
            FieldElement::new(2, modulus).unwrap(),
            FieldElement::new(3, modulus).unwrap(),
            FieldElement::new(4, modulus).unwrap(),
        ];
        let polynomial2 = Polynomial::new(coefficients2).unwrap();

        let sum = polynomial1.add(&polynomial2).unwrap();
        assert_eq!(sum.coefficients[0], FieldElement::new(3, modulus).unwrap());
        assert_eq!(sum.coefficients[1], FieldElement::new(5, modulus).unwrap());
        assert_eq!(sum.coefficients[2], FieldElement::new(4, modulus).unwrap());
    }

    #[test]
    fn test_sub() {
        let modulus = 97;
        // Define polynomial: 1 + 2x mod 97.
        let coefficients1 = vec![
            FieldElement::new(1, modulus).unwrap(),
            FieldElement::new(2, modulus).unwrap(),
        ];
        let polynomial1 = Polynomial::new(coefficients1).unwrap();

        // Define polynomial: 2 + 3x + 4x^2 mod 97.
        let coefficients2 = vec![
            FieldElement::new(2, modulus).unwrap(),
            FieldElement::new(3, modulus).unwrap(),
            FieldElement::new(4, modulus).unwrap(),
        ];
        let polynomial2 = Polynomial::new(coefficients2).unwrap();

        let diff = polynomial1.sub(&polynomial2).unwrap();
        assert_eq!(
            diff.coefficients[0],
            FieldElement::new(96, modulus).unwrap()
        );
        assert_eq!(
            diff.coefficients[1],
            FieldElement::new(96, modulus).unwrap()
        );
        assert_eq!(
            diff.coefficients[2],
            FieldElement::new(93, modulus).unwrap()
        );
    }

    #[test]
    fn test_mul() {
        let modulus = 97;
        // Define polynomial: 1 + 2x mod 97.
        let coefficients1 = vec![
            FieldElement::new(1, modulus).unwrap(),
            FieldElement::new(2, modulus).unwrap(),
        ];
        let polynomial1 = Polynomial::new(coefficients1).unwrap();

        // Define polynomial: 2 + 3x + 4x^2 mod 97.
        let coefficients2 = vec![
            FieldElement::new(2, modulus).unwrap(),
            FieldElement::new(3, modulus).unwrap(),
            FieldElement::new(4, modulus).unwrap(),
        ];
        let polynomial2 = Polynomial::new(coefficients2).unwrap();

        // Expected polynomial = 2 + 7x + 10x^2 + 8x^3 mod 97.
        let product = polynomial1.mul(&polynomial2).unwrap();
        assert_eq!(product.coefficients.len(), 4);
        assert_eq!(
            product.coefficients[0],
            FieldElement::new(2, modulus).unwrap()
        );
        assert_eq!(
            product.coefficients[1],
            FieldElement::new(7, modulus).unwrap()
        );
        assert_eq!(
            product.coefficients[2],
            FieldElement::new(10, modulus).unwrap()
        );
        assert_eq!(
            product.coefficients[3],
            FieldElement::new(8, modulus).unwrap()
        );
    }

    #[test]
    fn test_div() {
        let modulus = 97;
        // Define polynomial: 6 + 5x + x^2 mod 97.
        let coefficients1 = vec![
            FieldElement::new(6, modulus).unwrap(),
            FieldElement::new(5, modulus).unwrap(),
            FieldElement::new(1, modulus).unwrap(),
        ];
        let polynomial1 = Polynomial::new(coefficients1).unwrap();

        // Define polynomial: 2 + x mod 97.
        let coefficients2 = vec![
            FieldElement::new(2, modulus).unwrap(),
            FieldElement::new(1, modulus).unwrap(),
        ];
        let polynomial2 = Polynomial::new(coefficients2).unwrap();

        // Expected polynomial = 3 + x
        let (quotient, remainder) = polynomial1.div(&polynomial2).unwrap();

        assert_eq!(
            quotient.coefficients[0],
            FieldElement::new(3, modulus).unwrap()
        );
        assert_eq!(
            quotient.coefficients[1],
            FieldElement::new(1, modulus).unwrap()
        );

        for i in 0..2 {
            assert_eq!(
                remainder.coefficients[i],
                FieldElement::new(0, modulus).unwrap()
            );
        }
    }

    #[test]
    fn test_scale() {
        let modulus = 97;
        // Define polynomial: 1 + 2x + 3x^2 mod 97.
        let coefficients = vec![
            FieldElement::new(1, modulus).unwrap(),
            FieldElement::new(2, modulus).unwrap(),
            FieldElement::new(3, modulus).unwrap(),
        ];
        let polynomial = Polynomial::new(coefficients).unwrap();

        let scaled = polynomial
            .scale(&FieldElement::new(2, modulus).unwrap())
            .unwrap();

        assert_eq!(
            scaled.coefficients[0],
            FieldElement::new(2, modulus).unwrap()
        );
        assert_eq!(
            scaled.coefficients[1],
            FieldElement::new(4, modulus).unwrap()
        );
        assert_eq!(
            scaled.coefficients[2],
            FieldElement::new(6, modulus).unwrap()
        );
    }
}
