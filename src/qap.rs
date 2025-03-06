use crate::{
    circuit::ConstraintSystem, errors::ZKError, field::FieldElement, polynomial::Polynomial,
};

/// Represents R1CS constraints in QAP form.
pub struct QAP {
    // Interpolated polynomials for a, b, and c constraints.
    pub a_polynomials: Vec<Polynomial>,
    pub b_polynomials: Vec<Polynomial>,
    pub c_polynomials: Vec<Polynomial>,
    // Target polynomial.
    pub target_polynomial: Polynomial,
}

#[derive(Clone, Debug)]
struct Point {
    x: FieldElement,
    y: FieldElement,
}

impl QAP {
    /// Creates a new QAP using the provided R1CS.
    pub fn create(cs: &ConstraintSystem) -> Result<Self, ZKError> {
        let num_constraints = cs.constraints.len();
        if num_constraints == 0 {
            return Err(ZKError::PolynomialError("No constraints available.".into()));
        }

        let num_variables = cs.num_variables;
        let modulus = cs.constraints[0]
            .a
            .terms
            .get(0)
            .ok_or_else(|| ZKError::PolynomialError("Constraint has no terms.".into()))?
            .coefficient
            .modulus;

        // Get evaluation points.
        let evaluation_points: Vec<FieldElement> = (0..num_constraints)
            .map(|i| FieldElement::new((i + 1) as u64, modulus))
            .collect::<Result<_, _>>()?;

        // Construct the target polynomial.
        let mut target_polynomial = Polynomial::new(vec![FieldElement::new(1, modulus)?])?;
        for point in &evaluation_points {
            let factor = Polynomial::new(vec![
                FieldElement::new((modulus - (point.value % modulus)) % modulus, modulus)?,
                FieldElement::new(1, modulus)?,
            ])?;
            target_polynomial = target_polynomial.mul(&factor)?;
        }

        let mut a_polynomials = Vec::with_capacity(num_variables);
        let mut b_polynomials = Vec::with_capacity(num_variables);
        let mut c_polynomials = Vec::with_capacity(num_variables);

        for i in 0..num_variables {
            let mut a_points = Vec::with_capacity(num_constraints);
            let mut b_points = Vec::with_capacity(num_constraints);
            let mut c_points = Vec::with_capacity(num_constraints);

            for (j, constraint) in cs.constraints.iter().enumerate() {
                let r = evaluation_points[j].clone();
                let a_coefficient = constraint
                    .a
                    .terms
                    .iter()
                    .find(|term| term.index == i)
                    .map(|term| term.coefficient.clone())
                    .unwrap_or(FieldElement::new(0, modulus)?);
                let b_coefficient = constraint
                    .b
                    .terms
                    .iter()
                    .find(|term| term.index == i)
                    .map(|term| term.coefficient.clone())
                    .unwrap_or(FieldElement::new(0, modulus)?);
                let c_coefficient = constraint
                    .c
                    .terms
                    .iter()
                    .find(|term| term.index == i)
                    .map(|term| term.coefficient.clone())
                    .unwrap_or(FieldElement::new(0, modulus)?);

                a_points.push(Point {
                    x: r.clone(),
                    y: a_coefficient.clone(),
                });
                b_points.push(Point {
                    x: r.clone(),
                    y: b_coefficient.clone(),
                });
                c_points.push(Point {
                    x: r.clone(),
                    y: c_coefficient.clone(),
                });
            }

            let a_points_interpolated = interpolate_points(&a_points)?;
            let b_points_interpolated = interpolate_points(&b_points)?;
            let c_points_interpolated = interpolate_points(&c_points)?;

            a_polynomials.push(a_points_interpolated);
            b_polynomials.push(b_points_interpolated);
            c_polynomials.push(c_points_interpolated);
        }

        Ok(QAP {
            target_polynomial,
            a_polynomials,
            b_polynomials,
            c_polynomials,
        })
    }
}

// Interpolate points using Lagrange interpolation.
fn interpolate_points(points: &[Point]) -> Result<Polynomial, ZKError> {
    if points.is_empty() {
        return Err(ZKError::PolynomialError("No points to interpolate".into()));
    }

    let modulus = points[0].x.modulus;
    // Start with a zero polynomial.
    let mut result = Polynomial::new(vec![FieldElement::new(0, modulus)?])?;

    for (i, point_outer) in points.iter().enumerate() {
        let mut numerator = Polynomial::new(vec![FieldElement::new(1, modulus)?])?;
        let mut denominator = FieldElement::new(1, modulus)?;

        for (j, point_inner) in points.iter().enumerate() {
            if i == j {
                continue;
            }

            let numerator_factor = Polynomial::new(vec![
                FieldElement::new(
                    (modulus - (point_inner.x.value % modulus)) % modulus,
                    modulus,
                )?,
                FieldElement::new(1, modulus)?,
            ])?;
            numerator = numerator.mul(&numerator_factor)?;
            denominator = denominator.mul(&point_outer.x.sub(&point_inner.x)?)?;
        }

        let denominator_inverse = denominator.inv()?;
        let final_polynomial = numerator.mul(&Polynomial::new(vec![point_outer
            .y
            .mul(&denominator_inverse)?])?)?;
        result = result.add(&final_polynomial)?;
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use crate::{
        circuit::{ConstraintSystem, LinearCombination, R1CSConstraint, Term},
        field::FieldElement,
        polynomial::Polynomial,
    };

    use super::QAP;

    #[test]
    fn test_qap() {
        let modulus = 97;
        let mut cs = ConstraintSystem::new();

        // Let's consider the equation x^3 + x + 5 = 35.

        // We'll flatten this equation first:
        // v0 = 1
        // v1 = x
        // v2 (x^2) = x * x    -> Constraint 1
        // v3 (x^3)= x^2 * x   -> Constraint 2
        // v4 = v3 + v1        -> Constraint 3
        // v5 = v4 + 5         -> Constraint 4

        // Since we have 6 variables, the solution vector will be:
        // [v0, v1, v2, v3, v4, v5]

        // Let's allocate these 6 variables to the constraint system.
        let v0 = cs.allocate_variable();
        let v1 = cs.allocate_variable();
        let v2 = cs.allocate_variable();
        let v3 = cs.allocate_variable();
        let v4 = cs.allocate_variable();
        let v5 = cs.allocate_variable();

        // Now, let's add the constraints.
        // The constraints are represented in the form LC_a X LC_b = LC_c.

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

        // Helper function to check interpolation.
        let check_interpolation = |poly: &Polynomial, expected: &[u64]| {
            for (i, &exp_coeff) in expected.iter().enumerate() {
                let r = FieldElement::new((i + 1) as u64, modulus).unwrap();
                let eval = poly.evaluate(&r).unwrap();
                assert_eq!(
                    eval.value,
                    exp_coeff,
                    "Expected coefficient {} at r = {} but got {}",
                    exp_coeff,
                    i + 1,
                    eval.value
                );
            }
        };

        // Check poly_a for each variable.
        // There are 6 variables (indices 0 to 5).
        let expected_a = vec![
            vec![0, 0, 0, 5], // v0
            vec![1, 1, 1, 0], // v1
            vec![0, 0, 0, 0], // v2
            vec![0, 0, 1, 0], // v3
            vec![0, 0, 0, 1], // v4
            vec![0, 0, 0, 0], // v5
        ];
        for j in 0..6 {
            check_interpolation(&qap.a_polynomials[j], &expected_a[j]);
        }

        // Check poly_b for each variable.
        let expected_b = vec![
            vec![0, 0, 1, 1], // v0
            vec![1, 0, 0, 0], // v1
            vec![0, 1, 0, 0], // v2
            vec![0, 0, 0, 0], // v3
            vec![0, 0, 0, 0], // v4
            vec![0, 0, 0, 0], // v5
        ];
        for j in 0..6 {
            check_interpolation(&qap.b_polynomials[j], &expected_b[j]);
        }

        // Check poly_c for each variable.
        let expected_c = vec![
            vec![0, 0, 0, 0], // v0
            vec![0, 0, 0, 0], // v1
            vec![1, 0, 0, 0], // v2
            vec![0, 1, 0, 0], // v3
            vec![0, 0, 1, 0], // v4
            vec![0, 0, 0, 1], // v5
        ];
        for j in 0..6 {
            check_interpolation(&qap.c_polynomials[j], &expected_c[j]);
        }
    }
}
