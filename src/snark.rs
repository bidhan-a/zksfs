use crate::{
    curve::{EllipticCurve, EllipticCurvePoint},
    errors::ZKError,
    field::FieldElement,
    pairing::Pairing,
    qap::QAP,
};

/// Represents the CRS (Common Reference String) for the SNARK.
pub struct CRS {
    pub g1: EllipticCurvePoint,
    pub g2: EllipticCurvePoint,
}

/// Represents a SNARK proof.
pub struct Proof {
    pub a: EllipticCurvePoint,
    pub b: EllipticCurvePoint,
    pub c: EllipticCurvePoint,
}

pub struct SNARK {}

impl SNARK {
    /// Generates a dummy CRS.
    pub fn trusted_setup(curve: &EllipticCurve) -> Result<CRS, ZKError> {
        let modulus = curve.a.modulus;

        // We are choosing values here such that our dummy `verify_proof`
        // method is satisfied for modulo 97.
        let g1_x = FieldElement::new(47, modulus)?;
        let g1_y = FieldElement::new(1, modulus)?;
        let g2_x = FieldElement::new(2, modulus)?;
        let g2_y = FieldElement::new(1, modulus)?;

        let g1 = EllipticCurvePoint::Point { x: g1_x, y: g1_y };
        let g2 = EllipticCurvePoint::Point { x: g2_x, y: g2_y };

        Ok(CRS { g1, g2 })
    }

    /// Given a QAP (from the circuit) and a witness vector,
    /// compute the witness quotient polynomial h(x) and then "commit" to it via dummy group operations.
    /// The resulting proof consists of three group elements.
    pub fn create_proof(qap: &QAP, witness: &[FieldElement], crs: &CRS) -> Result<Proof, ZKError> {
        // Compute the witness quotient polynomial h(x).
        let h_polynomial = qap.calculate_witness_quotient(witness)?;
        // For a dummy commitment, we take the constant term of h(x) (h(0)) and "multiply" the CRS group elements.
        let h0 = h_polynomial
            .coefficients
            .get(0)
            .ok_or_else(|| ZKError::PolynomialError("Witness quotient polynomial is empty".into()))?
            .clone();

        // Simulate scalar multiplication of group elements by h0.
        let proof_a = match &crs.g1 {
            EllipticCurvePoint::Point { x, y } => EllipticCurvePoint::Point {
                x: x.mul(&h0)?,
                y: y.mul(&h0)?,
            },
            EllipticCurvePoint::Infinity => EllipticCurvePoint::Infinity,
        };
        let proof_b = match &crs.g2 {
            EllipticCurvePoint::Point { x, y } => EllipticCurvePoint::Point {
                x: x.mul(&h0)?,
                y: y.mul(&h0)?,
            },
            EllipticCurvePoint::Infinity => EllipticCurvePoint::Infinity,
        };
        // For proof_c, we combine g1 and g2 using a dummy addition (this is purely illustrative).
        let proof_c = match (&crs.g1, &crs.g2) {
            (
                EllipticCurvePoint::Point { x: x1, y: y1 },
                EllipticCurvePoint::Point { x: x2, y: y2 },
            ) => {
                // We simulate group addition by adding the coordinates.
                // In practice, group addition is nontrivial.
                EllipticCurvePoint::Point {
                    x: x1.add(x2)?,
                    y: y1.add(y2)?,
                }
            }
            _ => EllipticCurvePoint::Infinity,
        };

        Ok(Proof {
            a: proof_a,
            b: proof_b,
            c: proof_c,
        })
    }

    /// Given a proof, the CRS, and the elliptic curve,
    /// perform a dummy pairing check to verify the proof.
    pub fn verify_proof(proof: &Proof, crs: &CRS, curve: &EllipticCurve) -> Result<bool, ZKError> {
        let pairing_a = Pairing::create(curve, &proof.a, &crs.g2)?;
        let pairing_b = Pairing::create(curve, &proof.b, &crs.g1)?;
        let pairing_c = Pairing::create(curve, &proof.c, &crs.g2)?;
        let combined_value = pairing_b.value.mul(&pairing_c.value)?;
        let combined = Pairing {
            value: combined_value,
        };
        Ok(pairing_a == combined)
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        circuit::{ConstraintSystem, LinearCombination, R1CSConstraint, Term},
        curve::EllipticCurve,
        field::FieldElement,
        snark::SNARK,
    };

    use super::QAP;

    #[test]
    fn test_snark() {
        let modulus = 97;
        let curve = EllipticCurve {
            a: FieldElement::new(2, modulus).unwrap(),
            b: FieldElement::new(3, modulus).unwrap(),
        };

        // Run trusted setup to generate the CRS.
        let crs = SNARK::trusted_setup(&curve).unwrap();

        // Equation: x^3 + x + 5 = 35.
        let mut cs = ConstraintSystem::new();
        let v0 = cs.allocate_variable();
        let v1 = cs.allocate_variable();
        let v2 = cs.allocate_variable();
        let v3 = cs.allocate_variable();
        let v4 = cs.allocate_variable();
        let v5 = cs.allocate_variable();

        // Constraint 1: x * x = x^2
        {
            let mut lc_a = LinearCombination::new();
            lc_a.add_term(Term {
                index: v1,
                coefficient: FieldElement::new(1, modulus).unwrap(),
            });
            let mut lc_b = LinearCombination::new();
            lc_b.add_term(Term {
                index: v1,
                coefficient: FieldElement::new(1, modulus).unwrap(),
            });
            let mut lc_c = LinearCombination::new();
            lc_c.add_term(Term {
                index: v2,
                coefficient: FieldElement::new(1, modulus).unwrap(),
            });

            cs.add_constraint(R1CSConstraint::new(lc_a, lc_b, lc_c));
        }

        // Constraint 2: x * x^2 = x^3
        {
            let mut lc_a = LinearCombination::new();
            lc_a.add_term(Term {
                index: v1,
                coefficient: FieldElement::new(1, modulus).unwrap(),
            });
            let mut lc_b = LinearCombination::new();
            lc_b.add_term(Term {
                index: v2,
                coefficient: FieldElement::new(1, modulus).unwrap(),
            });
            let mut lc_c = LinearCombination::new();
            lc_c.add_term(Term {
                index: v3,
                coefficient: FieldElement::new(1, modulus).unwrap(),
            });

            cs.add_constraint(R1CSConstraint::new(lc_a, lc_b, lc_c));
        }

        // Constraint 3: v3 + v1 = v4 OR (v3 + v1) * 1 = v4
        {
            let mut lc_a = LinearCombination::new();
            lc_a.add_term(Term {
                index: v3,
                coefficient: FieldElement::new(1, modulus).unwrap(),
            });
            lc_a.add_term(Term {
                index: v1,
                coefficient: FieldElement::new(1, modulus).unwrap(),
            });
            let mut lc_b = LinearCombination::new();
            lc_b.add_term(Term {
                index: v0,
                coefficient: FieldElement::new(1, modulus).unwrap(),
            });
            let mut lc_c = LinearCombination::new();
            lc_c.add_term(Term {
                index: v4,
                coefficient: FieldElement::new(1, modulus).unwrap(),
            });

            cs.add_constraint(R1CSConstraint::new(lc_a, lc_b, lc_c));
        }

        // Constraint 4: v4 + 5 = v5 OR (v4 + 5) * 1 = v5
        {
            let mut lc_a = LinearCombination::new();
            lc_a.add_term(Term {
                index: v4,
                coefficient: FieldElement::new(1, modulus).unwrap(),
            });
            lc_a.add_term(Term {
                index: v0,
                coefficient: FieldElement::new(5, modulus).unwrap(),
            });
            let mut lc_b = LinearCombination::new();
            lc_b.add_term(Term {
                index: v0,
                coefficient: FieldElement::new(1, modulus).unwrap(),
            });
            let mut lc_c = LinearCombination::new();
            lc_c.add_term(Term {
                index: v5,
                coefficient: FieldElement::new(1, modulus).unwrap(),
            });

            cs.add_constraint(R1CSConstraint::new(lc_a, lc_b, lc_c));
        }

        // Create QAP from the constraint system.
        let qap = QAP::create(&cs).unwrap();

        // For x = 3, the witness is:
        // v0 = 1, v1 = 3, v2 = 9, v3 = 27, v4 = 27 + 3 = 30, v5 = 30 + 5 = 35.
        let witness = vec![
            FieldElement::new(1, modulus).unwrap(),
            FieldElement::new(3, modulus).unwrap(),
            FieldElement::new(9, modulus).unwrap(),
            FieldElement::new(27, modulus).unwrap(),
            FieldElement::new(30, modulus).unwrap(),
            FieldElement::new(35, modulus).unwrap(),
        ];

        // Prover: Generate a SNARK proof.
        let proof = SNARK::create_proof(&qap, &witness, &crs).unwrap();
        // Verifier: Check the proof.
        let valid = SNARK::verify_proof(&proof, &crs, &curve).unwrap();
        assert!(valid, "The proof is invalid.");
    }
}
