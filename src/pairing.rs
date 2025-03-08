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
