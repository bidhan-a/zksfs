use crate::{errors::ZKError, field::FieldElement};

/// Represents an elliptic curve defined by the equation:
/// y^2 = x^3 + ax + b (mod p)
#[derive(Debug, Clone)]
pub struct EllipticCurve {
    pub a: FieldElement,
    pub b: FieldElement,
}

/// Represents a point on the elliptic curve.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EllipticCurvePoint {
    Infinity,
    Point { x: FieldElement, y: FieldElement },
}

impl EllipticCurve {
    /// Check if the given point lies on the elliptic curve.
    pub fn is_on_curve(&self, point: &EllipticCurvePoint) -> Result<bool, ZKError> {
        match point {
            EllipticCurvePoint::Infinity => Ok(true),
            EllipticCurvePoint::Point { x, y } => {
                // Check if the point satisfies the elliptic curve equation.
                let y2 = y.mul(y)?;
                let x3 = x.mul(x)?.mul(x)?;
                let ax = x.mul(&self.a)?;
                let rhs = x3.add(&ax)?.add(&self.b)?;
                Ok(y2 == rhs)
            }
        }
    }

    /// Add two points on the elliptic curve.
    pub fn add_points(
        &self,
        p: &EllipticCurvePoint,
        q: &EllipticCurvePoint,
    ) -> Result<EllipticCurvePoint, ZKError> {
        match (p, q) {
            (EllipticCurvePoint::Infinity, _) => return Ok(q.clone()),
            (_, EllipticCurvePoint::Infinity) => return Ok(p.clone()),
            (
                EllipticCurvePoint::Point { x: x1, y: y1 },
                EllipticCurvePoint::Point { x: x2, y: y2 },
            ) => {
                if x1 == x2 {
                    if y1 == y2 && y1.value != 0 {
                        // Point doubling.

                        // slope(s) = (3x1^2 + a) / 2y1
                        let numerator = FieldElement::new(3, x1.modulus)?
                            .mul(&(x1.mul(x1)?))?
                            .add(&self.a)?;
                        let denominator = FieldElement::new(2, y1.modulus)?.mul(y1)?;
                        let slope = numerator.mul(&denominator.inv()?)?;

                        // x3 = s^2 - 2x1
                        let x3 = slope
                            .mul(&slope)?
                            .sub(&(&FieldElement::new(2, x1.modulus)?.mul(x1)?))?;

                        // y3 = s x (x1 - x3) - y1
                        let y3 = slope.mul(&(x1.sub(&x3))?)?.sub(y1)?;

                        Ok(EllipticCurvePoint::Point { x: x3, y: y3 })
                    } else {
                        // If the points are vertical reflections, their sum is Infinity (identity).
                        Ok(EllipticCurvePoint::Infinity)
                    }
                } else {
                    // Point addition.

                    // slope(s) = (y2 - y1) / (x2 - x1)
                    let numerator = y2.sub(y1)?;
                    let denominator = x2.sub(x1)?;
                    let slope = numerator.mul(&denominator.inv()?)?;

                    // x3 = s^2 - x1 - x2
                    let x3 = slope.mul(&slope)?.sub(x1)?.sub(x2)?;

                    // y3 = s x (x1 - x3) - y1
                    let y3 = slope.mul(&(x1.sub(&x3))?)?.sub(y1)?;

                    Ok(EllipticCurvePoint::Point { x: x3, y: y3 })
                }
            }
        }
    }

    /// Multiply a point with a scalar using the double-and-add algorithm.
    pub fn mul_scalar(
        &self,
        point: &EllipticCurvePoint,
        scalar: u64,
    ) -> Result<EllipticCurvePoint, ZKError> {
        let mut result = EllipticCurvePoint::Infinity;
        let mut addend = point.clone();
        let mut k = scalar;

        while k > 0 {
            if k & 1 == 1 {
                result = self.add_points(&result, &addend)?;
            }
            addend = self.add_points(&addend, &addend)?;
            k >>= 1;
        }

        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn get_test_values() -> (EllipticCurve, EllipticCurvePoint) {
        let modulus = 97;

        // Test curve with the equation y^2 = x^3 + 2x + 3.
        let curve = EllipticCurve {
            a: FieldElement::new(2, modulus).unwrap(),
            b: FieldElement::new(3, modulus).unwrap(),
        };

        // Test point that lies on the curve.
        let x = FieldElement::new(3, modulus).unwrap();
        let y = FieldElement::new(6, modulus).unwrap();
        let point = EllipticCurvePoint::Point { x, y };

        return (curve, point);
    }

    #[test]
    fn test_is_on_curve() {
        let (curve, point) = get_test_values();
        let on_curve = curve.is_on_curve(&point).unwrap();
        assert_eq!(on_curve, true);
    }

    #[test]
    fn test_add_points_identity() {
        let (curve, point) = get_test_values();
        let identity = EllipticCurvePoint::Infinity;
        let result = curve.add_points(&point, &identity).unwrap();
        assert_eq!(result, point);
    }

    #[test]
    fn test_mul_scalar_zero() {
        let (curve, point) = get_test_values();
        let result = curve.mul_scalar(&point, 0).unwrap();
        assert_eq!(result, EllipticCurvePoint::Infinity);
    }

    #[test]
    fn test_mul_scalar() {
        let (curve, point) = get_test_values();
        // P + P
        let double = curve.add_points(&point, &point).unwrap();
        // 2P
        let mul_scalar_result = curve.mul_scalar(&point, 2).unwrap();
        // P + P = 2P
        assert_eq!(double, mul_scalar_result,);
    }
}
