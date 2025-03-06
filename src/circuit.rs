use crate::{errors::ZKError, field::FieldElement};

/// Represents a term i.e. a variable with a coefficient at an index.
#[derive(Clone, Debug)]
pub struct Term {
    pub index: usize,
    pub coefficient: FieldElement,
}

/// Represents a linear combination of terms.
#[derive(Clone, Debug)]
pub struct LinearCombination {
    pub terms: Vec<Term>,
}

impl LinearCombination {
    /// Creates a new, empty linear combinaton.
    pub fn new() -> Self {
        LinearCombination { terms: Vec::new() }
    }

    /// Add a term.
    pub fn add_term(&mut self, term: Term) {
        self.terms.push(term);
    }

    /// Evaluates the linear combination given a witness victor.
    /// Each variable's value is taken from the witness by its index.
    pub fn evaluate(&self, witness: &[FieldElement]) -> Result<FieldElement, ZKError> {
        if witness.is_empty() {
            return Err(ZKError::CircuitError("Witness vector is empty.".into()));
        }

        let modulus = witness[0].modulus;
        let mut result = FieldElement::new(0, modulus)?;
        for term in &self.terms {
            if term.index >= witness.len() {
                return Err(ZKError::CircuitError(
                    "Witness index {} is out of bounds.".into(),
                ));
            }
            let term_value = term.coefficient.mul(&witness[term.index])?;
            result = result.add(&term_value)?;
        }

        Ok(result)
    }
}

/// Represents a R1CS constraint which is defined as:
/// (LinearCombination a) x (LinearCombination b) = (LinearCombination c)
#[derive(Clone, Debug)]
pub struct R1CSConstraint {
    pub a: LinearCombination,
    pub b: LinearCombination,
    pub c: LinearCombination,
}

impl R1CSConstraint {
    /// Creates a new R1CSConstraint.
    pub fn new(a: LinearCombination, b: LinearCombination, c: LinearCombination) -> Self {
        R1CSConstraint { a, b, c }
    }
}

/// Stores a set of R1CS constraints and the number of variables.
#[derive(Clone, Debug)]
pub struct ConstraintSystem {
    pub constraints: Vec<R1CSConstraint>,
    pub num_variables: usize,
}

impl ConstraintSystem {
    /// Creates a new, empty constraint system.
    pub fn new() -> Self {
        ConstraintSystem {
            constraints: Vec::new(),
            num_variables: 0,
        }
    }

    /// Adds a new R1CS constraint.
    pub fn add_constraint(&mut self, constraint: R1CSConstraint) {
        self.constraints.push(constraint);
    }

    /// Allocates a new variable and returns its index.
    pub fn allocate_variable(&mut self) -> usize {
        let var_index = self.num_variables;
        self.num_variables += 1;
        var_index
    }

    /// Evaluates the provided witness against all constraints.
    /// For each constraint, it checks that LC a (witness) x LC b (witness) = LC c (witness).
    pub fn evaluate(&self, witness: &[FieldElement]) -> Result<bool, ZKError> {
        for (i, constraint) in self.constraints.iter().enumerate() {
            let a_val = constraint.a.evaluate(witness)?;
            let b_val = constraint.b.evaluate(witness)?;
            let c_val = constraint.c.evaluate(witness)?;
            let product = a_val.mul(&b_val)?;
            if product != c_val {
                return Err(ZKError::CircuitError(format!(
                    "Constraint {} not satisfied: {:?} x {:?} != {:?}",
                    i, a_val, b_val, c_val
                )));
            }
        }

        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::FieldElement;

    #[test]
    fn test_linear_combination() {
        // Create a linear combination: 3v0 + 5v1.
        let modulus = 97;
        let mut lc = LinearCombination::new();
        lc.add_term(Term {
            index: 0,
            coefficient: FieldElement::new(3, modulus).unwrap(),
        });
        lc.add_term(Term {
            index: 1,
            coefficient: FieldElement::new(5, modulus).unwrap(),
        });

        // Witness vector: v0 = 2, v1 = 4.
        let witness = vec![
            FieldElement::new(2, modulus).unwrap(),
            FieldElement::new(4, modulus).unwrap(),
        ];
        let result = lc.evaluate(&witness).unwrap();
        // Expected evaluation: 3 x 2 + 5 x 4 = 26 mod 97.
        assert_eq!(result, FieldElement::new(26, modulus).unwrap());
    }

    #[test]
    fn test_constraint_system() {
        let modulus = 97;
        let mut cs = ConstraintSystem::new();

        // Allocate variables: v0, v1, v2.
        let v0 = cs.allocate_variable(); // index 0
        let v1 = cs.allocate_variable(); // index 1
        let v2 = cs.allocate_variable(); // index 2

        // Enforce the constraint: v0 * v1 = v2.
        // Linear combinations for the constraint:
        // a = 1 * v0, b = 1 * v1, c = 1 * v2.
        let mut a_lc = LinearCombination::new();
        a_lc.add_term(Term {
            index: v0,
            coefficient: FieldElement::new(1, modulus).unwrap(),
        });

        let mut b_lc = LinearCombination::new();
        b_lc.add_term(Term {
            index: v1,
            coefficient: FieldElement::new(1, modulus).unwrap(),
        });

        let mut c_lc = LinearCombination::new();
        c_lc.add_term(Term {
            index: v2,
            coefficient: FieldElement::new(1, modulus).unwrap(),
        });

        let constraint = R1CSConstraint::new(a_lc, b_lc, c_lc);
        cs.add_constraint(constraint);

        // Choose a witness such that v0 * v1 = v2.
        // Let v0 = 3, v1 = 4, then v2 should be 12 mod 97.
        let witness = vec![
            FieldElement::new(3, modulus).unwrap(),
            FieldElement::new(4, modulus).unwrap(),
            FieldElement::new(12, modulus).unwrap(),
        ];
        let result = cs.evaluate(&witness).unwrap();
        assert_eq!(result, true);
    }
}
