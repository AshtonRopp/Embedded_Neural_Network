# 🧠 SystemVerilog CNN Accelerator with Training Support (INT16)

This project implements a Convolutional Neural Network (CNN) accelerator in **SystemVerilog**, capable of **training and inference** using **INT16 fixed-point format** and the **ReLU activation function**.

---

## 🚀 Project Overview

- **Network Type**: Convolutional Neural Network (CNN)
- **Target Use**: Inference + Training
- **Data Format**: INT16 fixed-point (Q8.8)
- **Activation**: ReLU
- **Training Algorithm**: Stochastic Gradient Descent (SGD)
- **Loss Function**: Mean Squared Error (MSE)

---

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

