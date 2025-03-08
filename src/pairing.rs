use crate::{
    curve::{EllipticCurve, EllipticCurvePoint},
    errors::ZKError,
    field::FieldElement,
};

/// Represents the result of a pairing operation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Pairing {
    pub value: FieldElement,
}

impl Pairing {
    /// Creates a dummy pairing.
    ///
    /// # Parameters
    /// - `curve`: A reference to the elliptic curve.
    /// - `p`: A point from group G1.
    /// - `q`: A point from group G2.
    ///
    /// # Returns
    /// A `Pairing` that:
    ///   - If either point is the point at infinity, returns 1 (the identity in the field).
    ///   - Otherwise, returns the product of the x‑coordinates of `p` and `q` modulo the field's modulus.
    pub fn create(
        curve: &EllipticCurve,
        p: &EllipticCurvePoint,
        q: &EllipticCurvePoint,
    ) -> Result<Self, ZKError> {
        match (p, q) {
            (EllipticCurvePoint::Infinity, _) | (_, EllipticCurvePoint::Infinity) => {
                // If either point is at infinity, the pairing is defined as the identity (1).
                let value = FieldElement::new(1, curve.a.modulus)?;
                Ok(Pairing { value })
            }
            (
                EllipticCurvePoint::Point { x: x1, y: _ },
                EllipticCurvePoint::Point { x: x2, y: _ },
            ) => {
                // Otherwise, multiply the x-coordinates.
                let value = x1.mul(x2)?;
                Ok(Pairing { value })
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve::EllipticCurve;
    use crate::field::FieldElement;

    #[test]
    fn test_pairing() {
        let modulus = 97;
        let curve = EllipticCurve {
            a: FieldElement::new(2, modulus).unwrap(),
            b: FieldElement::new(3, modulus).unwrap(),
        };
        let point_a = EllipticCurvePoint::Point {
            x: FieldElement::new(3, modulus).unwrap(),
            y: FieldElement::new(6, modulus).unwrap(),
        };
        let point_b = EllipticCurvePoint::Point {
            x: FieldElement::new(2, modulus).unwrap(),
            y: FieldElement::new(5, modulus).unwrap(),
        };
        let pairing = Pairing::create(&curve, &point_a, &point_b).unwrap();
        // Dummy pairing multiplies the x-coordinates.
        // For p and q, x = 3, so expected result is 3 * 2 = 6 mod 97.
        assert_eq!(pairing.value, FieldElement::new(6, modulus).unwrap());
    }
}
