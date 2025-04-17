# SystemVerilog CNN Accelerator with Training Support

This project implements a Convolutional Neural Network (CNN) accelerator in **SystemVerilog**, capable of **training and inference** using **INT16 fixed-point format** and the **ReLU activation function**.

---
## Project Highlights
- **Utilized Physical Design Software**
  - Synopsis Design Compiler for preliminary area metrics
  - OpenROAD for PNR, clock tree, and PPA metrics
  - Vivado and ModelSim for simulation

- **Custom RTL CNN Architecture**
  - Includes 4 parallel 3×3 convolution filters
  - ReLU activation stage
  - Flattening layer
  - Fully connected layers with full 8-neuron FC1 and scalar FC2

- **Loop Unrolling and Acceleration**
  - Kernel flattening loop unrolled
  - Parallel processing via MAC units

- **Parameterizable & Modular Design**
  - Modular convolution, activation, FC, and backpropagation blocks
  - External weight/bias management for flexible training and inference

- **Training Support in RTL**
  - Implements mean squared error (MSE) loss
  - Gradient calculation and backpropagation for both FC1 and FC2

- **Optimized MAC Unit**
  - Tree-structured reduction with pipelined accumulation
  - Parameterized for scalable depth (`MAC_DEPTH`)

## Functional Blocks

### Forward Pass
- Input Buffer
- Convolution Unit (MAC array)
- ReLU Activation Unit
- Fully Connected (FC) Layer
- Output Buffer

### Backward Pass (Training)
- Error Computation (Loss Function)
- Backpropagation through FC layer
- Backpropagation through Conv layer
- ReLU Derivative Unit
- Gradient Buffers

### Weight Update
- Learning Rate Control
- Weight Update Logic: `w -= η * ∂L/∂w`
- Bias Update Logic

---

## Control and Pipelining

- FSM for Training Steps:
  - Load Inputs
  - Forward Pass
  - Error Computation
  - Backward Pass
  - Weight Update
- Pipelined MAC Units

---

## Interfaces

- Input/Output Interfaces:
  - Load Input Data and Labels
  - Initialize Weights and Biases
  - Read Output Predictions and Weights
  - Training Control Signals (Start, Reset, Learning Rate)

---

## Fixed-Point Arithmetic

- Format: **INT16 Q8.8**
  - 8 integer bits, 8 fractional bits
- Scaling applied during:
  - Multiplication
  - Activation
  - Gradient updates

---

## Verification and Testing

- System-Level Testbench
  - Both classes included in each epochs
  - Accurate predictions occur after ~18 epochs


## Future Optimizations

- True pipelining
- Clock gating and resource reuse
- Memory tiling for large inputs

---


## Results
### Optimization Part A
- x% reduction in power
- x% reduction in size
- WNS from -2.52 to 0.17

## Resources

- https://github.com/The-OpenROAD-Project/OpenROAD-flow-scripts/blob/master/docs/tutorials/FlowTutorial.md
- https://github.com/zachjs/sv2v