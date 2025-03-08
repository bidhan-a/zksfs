# zk-SNARK from scratch

This repository contains a pedagogical implementation of a zk-SNARK system from scratch. The goal of this project is to understand the underlying mathematical concepts and how they translate into real code. 

## Overview of Modules

The implementation is structured into several modular components:

### 1. `field.rs`
- Implements finite field arithmetic, including addition, multiplication, and modular inverses.
- Uses modular arithmetic to ensure all operations are within a prime field.

### 2. `curve.rs`
- Implements an elliptic curve over a finite field.
- Supports basic elliptic curve operations such as point addition and scalar multiplication.

### 3. `circuit.rs`
- Represents arithmetic circuits using R1CS constraints.
- Allows defining computations as a set of constraints on variables.

### 4. `polynomial.rs`
- Implements polynomial arithmetic over finite fields.
- Supports polynomial evaluation, and operations like addition, subtraction, multiplication, and division.

### 5. `qap.rs`
- Converts circuit constraints into a Quadratic Arithmetic Program (QAP).
- Uses Lagrange polynomial interpolation to construct QAP polynomials.

### 6. `pairing.rs`
- Implements a simple bilinear pairing function.
- Used in the zk-SNARK verification step.

### 7. `snark.rs`
- Implements the zk-SNARK protocol including:
  - **Trusted Setup:** Generates Common Reference String (CRS).
  - **Prover:** Constructs a proof given a witness.
  - **Verifier:** Verifies the proof using bilinear pairings.

## Features
- **Fully tested**: All modules have full test coverage.
- **Educational focus**: Code is structured for clarity and learning, not for performance or security.
- **Step-by-step construction**: Each module builds on the previous ones, following the theoretical foundations of zk-SNARKs.

## Disclaimer
⚠️ **This is not a production-grade implementation.** It is for learning purposes only. Security, performance, and cryptographic best practices have not been considered for real-world usage.

## Getting Started
To build:
```sh
cargo build
```

To run tests:
```sh
cargo test
```


## Future Updates
Currently, some components, such as pairing-based cryptography and proof creation and verification, use simplified or dummy implementations. These will be improved in future updates with more robust and mathematically sound implementations.

