# 🧠 SystemVerilog CNN Accelerator with Training Support (INT16)

This project implements a Convolutional Neural Network (CNN) accelerator in **SystemVerilog**, capable of **training and inference** using **INT16 fixed-point format** and the **ReLU activation function**.

---
## 🚀 Project Highlights: Modular CNN Hardware Accelerator

- **Custom RTL CNN Architecture**  
  - Includes 4 parallel 3×3 convolution filters  
  - ReLU activation stage  
  - Flattening layer  
  - Fully connected layers with full 8-neuron FC1 and scalar FC2  

- **Loop Unrolling and Acceleration**  
  - MAC product and accumulation loops fully unrolled  
  - Kernel flattening loop unrolled  
  - Prepares for future pipelined double-buffering  

- **Parameterizable & Modular Design**  
  - Modular convolution, activation, FC, and backprop blocks  
  - External weight/bias management for flexible training and inference  

- **Training Support in RTL**  
  - Implements mean squared error (MSE) loss  
  - Gradient calculation and backpropagation for both FC1 and FC2  
  - Learning rate is externally tunable  

- **Optimized MAC Unit**  
  - Tree-structured reduction with pipelined accumulation  
  - Parameterized for scalable depth (`MAC_DEPTH`)  

## 🔹 Functional Blocks

### ✅ Forward Pass
- [ ] Input Buffer
- [ ] Convolution Unit (MAC array)
- [ ] ReLU Activation Unit
- [ ] Pooling Layer (optional)
- [ ] Fully Connected (FC) Layer
- [ ] Output Buffer

### ✅ Backward Pass (Training)
- [ ] Error Computation (Loss Function)
- [ ] Backpropagation through FC layer
- [ ] Backpropagation through Conv layer
- [ ] ReLU Derivative Unit
- [ ] Gradient Buffers

### ✅ Weight Update
- [ ] Learning Rate Control
- [ ] Weight Update Logic: `w -= η * ∂L/∂w`
- [ ] Bias Update Logic

---

## 🔁 Control and Pipelining

- [ ] FSM for Training Steps:
  - Load Inputs
  - Forward Pass
  - Error Computation
  - Backward Pass
  - Weight Update
- [ ] Pipelined MAC Units

---

## 💾 Memory Architecture

- [ ] Weight Memory
- [ ] Bias Memory
- [ ] Gradient Memory
- [ ] Intermediate Activation Buffers

---

## 🔌 Interfaces

- [ ] Input/Output Interfaces:
  - Load Input Data and Labels
  - Initialize Weights and Biases
  - Read Output Predictions
  - Training Control Signals (Start, Reset, Learning Rate)
- [ ] Debug Signals

---

## 🧮 Fixed-Point Arithmetic

- Format: **INT16 Q8.8**
  - 8 integer bits, 8 fractional bits
- Scaling applied during:
  - Multiplication
  - Activation
  - Gradient updates

---

## 🧪 Verification and Testing

- [ ] Unit Testbenches for:
  - MAC
  - ReLU
  - Convolution
  - Gradient Update
- [ ] System-Level Testbench
  - Example: XOR or 2×2 image task
- [ ] Compare with Python (NumPy/PyTorch) reference implementation

---

## 🧠 First Milestone

Build a simple CNN with:
- Single Conv Layer → ReLU → FC Layer
- INT16 computations
- Forward and Backward Pass
- Gradient Descent Updates

---

## 📈 Future Optimizations

- Loop unrolling
- MAC pipelining
- Clock gating and resource reuse
- Memory tiling for large inputs

---


## Resources

- https://github.com/The-OpenROAD-Project/OpenROAD-flow-scripts/blob/master/docs/tutorials/FlowTutorial.md
- https://github.com/zachjs/sv2v