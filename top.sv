// fixed_point_multiplier.sv
module fixed_point_multiplier (
    input  logic signed [15:0] a,    // Q8.8
    input  logic signed [15:0] b,    // Q8.8
    output logic signed [15:0] result  // Q8.8 result
);
    logic signed [31:0] mult_full;
    
    always_comb begin
        mult_full = a * b;             // Q8.8 × Q8.8 = Q16.16
        result    = mult_full[23:8];   // Convert back to Q8.8 (right shift by 8 bits)
    end
endmodule

// relu_activation.sv
module relu_activation (
    input  logic signed [15:0] in,   // Q8.8
    output logic signed [15:0] out   // Q8.8
);
    always_comb begin
        out = (in[15] == 1) ? 16'sd0 : in; // if negative, output 0
    end
endmodule

// mac_unit.sv
module mac_unit #(
    parameter MAC_DEPTH = 4
) (
    input  logic                  clk,
    input  logic                  rst,
    input  logic                  enable,
    input  logic signed [15:0]    a [MAC_DEPTH],  // Input A vector (Q8.8)
    input  logic signed [15:0]    b [MAC_DEPTH],  // Input B vector (Q8.8)
    output logic signed [15:0]    result          // Accumulated Q8.8 output
);
    logic signed [31:0] products [MAC_DEPTH];
    logic signed [31:0] sum;

    always_comb begin
        for (int i = 0; i < MAC_DEPTH; i++) begin
            products[i] = a[i] * b[i]; // Q8.8 × Q8.8 = Q16.16
        end
        sum = 0;
        for (int i = 0; i < MAC_DEPTH; i++) begin
            sum += products[i];
        end
    end

    // Sequential logic to register output
    always_ff @(posedge clk or posedge rst) begin
        if (rst)
            result <= 16'sd0;
        else if (enable)
            result <= sum[23:8]; // Convert from Q16.16 to Q8.8
    end
endmodule

module conv2d_unit_pipelined #(
    parameter IN_SIZE     = 4,
    parameter KERNEL_SIZE = 3
) (
    input  logic clk,
    input  logic rst,
    input  logic start,
    input  var logic signed [15:0] input_feature [IN_SIZE][IN_SIZE],
    input  var logic signed [15:0] kernel_weights [KERNEL_SIZE][KERNEL_SIZE],
    output logic done,
    output logic signed [15:0] output_feature [IN_SIZE-KERNEL_SIZE+1][IN_SIZE-KERNEL_SIZE+1]
);

    localparam OUT_SIZE  = IN_SIZE - KERNEL_SIZE + 1;
    localparam MAC_DEPTH = KERNEL_SIZE * KERNEL_SIZE;

    logic signed [15:0] a_flat [MAC_DEPTH];
    logic signed [15:0] b_flat [MAC_DEPTH];
    logic signed [15:0] mac_result;
    logic mac_enable;

    // -------------------------------
    // Flatten kernel weights
    // -------------------------------
    int k;
    always_comb begin
        k = 0;
        for (int m = 0; m < KERNEL_SIZE; m++) begin
            for (int n = 0; n < KERNEL_SIZE; n++) begin
                b_flat[k] = kernel_weights[m][n];
                k++;
            end
        end
    end

    // -------------------------------
    // FSM control logic
    // -------------------------------
    typedef enum logic [1:0] {
        IDLE,
        LOAD_PATCH,
        COMPUTE,
        WRITE_RESULT
    } state_t;

    state_t state;
    int x_idx, y_idx;

    always_ff @(posedge clk or posedge rst) begin
        if (rst) begin
            state       <= IDLE;
            done        <= 0;
            x_idx       <= 0;
            y_idx       <= 0;
            mac_enable  <= 0;
        end else begin
            case (state)
                IDLE: begin
                    done <= 0;
                    if (start) begin
                        x_idx <= 0;
                        y_idx <= 0;
                        state <= LOAD_PATCH;
                    end
                end

                LOAD_PATCH: begin
                    for (int m = 0; m < KERNEL_SIZE; m++) begin
                        for (int n = 0; n < KERNEL_SIZE; n++) begin
                            a_flat[m * KERNEL_SIZE + n] = input_feature[x_idx + m][y_idx + n];
                        end
                    end
                    mac_enable <= 1;
                    state <= COMPUTE;
                end

                COMPUTE: begin
                    mac_enable <= 0;  // 1-cycle pulse
                    state <= WRITE_RESULT;
                end

                WRITE_RESULT: begin
                    output_feature[x_idx][y_idx] <= mac_result;

                    if (y_idx < OUT_SIZE - 1) begin
                        y_idx <= y_idx + 1;
                        state <= LOAD_PATCH;
                    end else if (x_idx < OUT_SIZE - 1) begin
                        x_idx <= x_idx + 1;
                        y_idx <= 0;
                        state <= LOAD_PATCH;
                    end else begin
                        state <= IDLE;
                        done <= 1;
                    end
                end
            endcase
        end
    end

    // -------------------------------
    // MAC Unit
    // -------------------------------
    mac_unit #(.MAC_DEPTH(MAC_DEPTH)) mac_inst (
        .clk(clk),
        .rst(rst),
        .enable(mac_enable),
        .a(a_flat),
        .b(b_flat),
        .result(mac_result)
    );

endmodule


// relu_layer_2d.sv
module relu_layer_2d #(
    parameter HEIGHT = 2,
    parameter WIDTH  = 2
) (
    input  logic signed [15:0] in_feature  [HEIGHT][WIDTH], // Q8.8
    output logic signed [15:0] out_feature [HEIGHT][WIDTH]  // Q8.8
);
    always_comb begin
        for (int i = 0; i < HEIGHT; i++) begin
            for (int j = 0; j < WIDTH; j++) begin
                out_feature[i][j] = (in_feature[i][j][15] == 1) ? 16'sd0 : in_feature[i][j];
            end
        end
    end
endmodule

// fc_layer.sv
module fc_layer #(
    parameter INPUT_DIM = 4  // e.g., 4 inputs → 1 output
) (
    input  logic                  clk,
    input  logic                  rst,
    input  logic                  start,
    input  logic signed [15:0]    input_vec [INPUT_DIM],  // Q8.8
    input  logic signed [15:0]    weights   [INPUT_DIM],  // Q8.8
    input  logic signed [15:0]    bias,                   // Q8.8
    output logic                  done,
    output logic signed [15:0]    output_val              // Q8.8
);

    logic mac_enable;
    logic signed [15:0] mac_result;

    typedef enum logic [1:0] {
        IDLE,
        COMPUTE,
        ADD_BIAS,
        DONE
    } state_t;

    state_t state;

    always_ff @(posedge clk or posedge rst) begin
        if (rst) begin
            state <= IDLE;
            mac_enable <= 0;
            done <= 0;
        end else begin
            case (state)
                IDLE: begin
                    done <= 0;
                    if (start) begin
                        mac_enable <= 1;
                        state <= COMPUTE;
                    end
                end

                COMPUTE: begin
                    mac_enable <= 0; // 1-cycle pulse
                    state <= ADD_BIAS;
                end

                ADD_BIAS: begin
                    output_val <= mac_result + bias;
                    state <= DONE;
                end

                DONE: begin
                    done <= 1;
                    state <= IDLE;
                end
            endcase
        end
    end

    mac_unit #(.MAC_DEPTH(INPUT_DIM)) mac_inst (
        .clk(clk),
        .rst(rst),
        .enable(mac_enable),
        .a(input_vec),
        .b(weights),
        .result(mac_result)
    );
endmodule

// loss_gradient.sv
module loss_gradient (
    input  logic signed [15:0] prediction,  // Q8.8
    input  logic signed [15:0] label,       // Q8.8
    output logic signed [15:0] dL_dout      // Q8.8
);
    // ∂L/∂output = output - label (MSE derivative: (ŷ - y))
    assign dL_dout = prediction - label;
endmodule

// fc_backprop.sv
module fc_backprop #(
    parameter INPUT_DIM = 4
) (
    input  logic clk,
    input  logic rst,
    input  logic start,
    input  logic signed [15:0] input_vec   [INPUT_DIM],  // Q8.8
    input  logic signed [15:0] dL_dout,                  // Q8.8
    input  logic signed [15:0] learning_rate,            // Q8.8
    input  logic signed [15:0] weights_in   [INPUT_DIM], // Original weights
    output logic signed [15:0] weights_out  [INPUT_DIM], // Updated weights
    input  logic signed [15:0] bias_in,                  // Original bias
    output logic signed [15:0] bias_out,                 // Updated bias
    output logic signed [15:0] dL_drelu     [INPUT_DIM], // ∂L/∂ReLU input
    output logic done
);
    typedef enum logic [1:0] {
        IDLE,
        COMPUTE,
        DONE
    } state_t;

    state_t state;
    logic [1:0] idx;
    logic update_bias_done;

    logic signed [31:0] grad_mul     [INPUT_DIM];
    logic signed [31:0] update_mul   [INPUT_DIM];
    logic signed [31:0] backprop_mul [INPUT_DIM];
    logic signed [31:0] bias_update;

    always_comb begin
        for (int i = 0; i < INPUT_DIM; i++) begin
            grad_mul[i]     = input_vec[i] * dL_dout;
            update_mul[i]   = learning_rate * grad_mul[i][23:8];
            backprop_mul[i] = weights_in[i] * dL_dout;
        end
        bias_update = learning_rate * dL_dout;
    end

    always_ff @(posedge clk or posedge rst) begin
        if (rst) begin
            state <= IDLE;
            idx   <= 0;
            done  <= 0;
            update_bias_done <= 0;
        end else begin
            case (state)
                IDLE: begin
                    done <= 0;
                    idx  <= 0;
                    update_bias_done <= 0;
                    state <= start ? COMPUTE : IDLE;
                end

                COMPUTE: begin
                    weights_out[idx] <= weights_in[idx] - update_mul[idx][23:8];
                    dL_drelu[idx]    <= backprop_mul[idx][23:8];

                    if (idx == INPUT_DIM - 1) begin
                        bias_out <= bias_in - bias_update[23:8];
                        state <= DONE;
                    end else begin
                        idx <= idx + 1;
                    end
                end

                DONE: begin
                    done <= 1;
                    state <= IDLE;
                end
            endcase
        end
    end
endmodule



// top_forward_pipeline.sv
module top_forward_pipeline (
    input  logic                  clk,
    input  logic                  rst,
    input  logic                  start,
    input  logic signed [15:0]    input_feature  [4][4],  // Q8.8
    input  logic signed [15:0]    kernel_weights  [3][3], // Q8.8
    input  logic signed [15:0]    fc_weights      [3:0],  // Q8.8 (for 2x2 ReLU output → 1 FC output)
    input  logic signed [15:0]    fc_bias,                // Q8.8
    output logic                  done,
    output logic signed [15:0]    output_value            // Q8.8
);

    // Internal connections
    logic                  conv_done;
    logic                  relu_done;
    logic                  fc_done;
    logic signed [15:0]    conv_out  [2][2];  // Output of Conv2D (4x4 input, 3x3 kernel) → 2x2
    logic signed [15:0]    relu_out  [2][2];  // ReLU result
    logic signed [15:0]    relu_flat [3:0];   // Flattened ReLU result to vector for FC

    // ----------------------------------------
    // Convolution Unit (pipelined)
    // ----------------------------------------
    conv2d_unit_pipelined conv_layer (
        .clk(clk),
        .rst(rst),
        .start(start),
        .input_feature(input_feature),
        .kernel_weights(kernel_weights),
        .done(conv_done),
        .output_feature(conv_out)
    );

    // ----------------------------------------
    // ReLU Unit
    // ----------------------------------------
    relu_layer_2d #(.HEIGHT(2), .WIDTH(2)) relu_unit (
        .in_feature(conv_out),
        .out_feature(relu_out)
    );

    // ----------------------------------------
    // Flatten ReLU output for FC input
    // ----------------------------------------
    always_comb begin
        relu_flat[0] = relu_out[0][0];
        relu_flat[1] = relu_out[0][1];
        relu_flat[2] = relu_out[1][0];
        relu_flat[3] = relu_out[1][1];
    end

    // ----------------------------------------
    // Fully Connected Layer
    // ----------------------------------------
    fc_layer #(.INPUT_DIM(4)) fc (
        .clk(clk),
        .rst(rst),
        .start(conv_done),  // Start FC once conv is done
        .input_vec(relu_flat),
        .weights(fc_weights),
        .bias(fc_bias),
        .done(fc_done),
        .output_val(output_value)
    );

    // ----------------------------------------
    // Done signal
    // ----------------------------------------
    assign done = fc_done;

endmodule

module conv_kernel_grad #(
    parameter IN_SIZE = 4,
    parameter KERNEL_SIZE = 3
) (
    input  logic clk,
    input  logic rst,
    input  logic start,
    input  logic signed [15:0] input_feature [IN_SIZE][IN_SIZE],   // Q8.8
    input  logic signed [15:0] dL_dconv_out [IN_SIZE-KERNEL_SIZE+1][IN_SIZE-KERNEL_SIZE+1], // Q8.8
    output logic signed [15:0] dL_dkernel [KERNEL_SIZE][KERNEL_SIZE], // Q8.8
    output logic done
);

    localparam OUT_SIZE = IN_SIZE - KERNEL_SIZE + 1;

    typedef enum logic [1:0] { IDLE, COMPUTE, DONE } state_t;
    state_t state;

    // Accumulators
    logic signed [31:0] acc [KERNEL_SIZE][KERNEL_SIZE];

    // Loop indices
    logic [$clog2(KERNEL_SIZE)-1:0] m, n;
    logic [$clog2(OUT_SIZE)-1:0] i, j;

    always_ff @(posedge clk or posedge rst) begin
        if (rst) begin
            state <= IDLE;
            done  <= 0;
            m <= 0;
            n <= 0;
            i <= 0;
            j <= 0;

            for (int mm = 0; mm < KERNEL_SIZE; mm++) begin
                for (int nn = 0; nn < KERNEL_SIZE; nn++) begin
                    acc[mm][nn] <= 32'sd0;
                end
            end

        end else begin
            case (state)
                IDLE: begin
                    done <= 0;
                    m <= 0;
                    n <= 0;
                    i <= 0;
                    j <= 0;

                    for (int mm = 0; mm < KERNEL_SIZE; mm++) begin
                        for (int nn = 0; nn < KERNEL_SIZE; nn++) begin
                            acc[mm][nn] <= 32'sd0;
                        end
                    end

                    if (start)
                        state <= COMPUTE;
                end

                COMPUTE: begin
                    // Accumulate gradient: acc[m][n] += dout[i][j] * input[i+m][j+n]
                    acc[m][n] <= acc[m][n] + 
                        dL_dconv_out[i][j] * input_feature[i + m][j + n];

                    // Update i, j, n, m counters
                    if (j < OUT_SIZE - 1) begin
                        j <= j + 1;
                    end else begin
                        j <= 0;
                        if (i < OUT_SIZE - 1) begin
                            i <= i + 1;
                        end else begin
                            i <= 0;
                            if (n < KERNEL_SIZE - 1) begin
                                n <= n + 1;
                            end else begin
                                n <= 0;
                                if (m < KERNEL_SIZE - 1) begin
                                    m <= m + 1;
                                end else begin
                                    state <= DONE;
                                end
                            end
                        end
                    end
                end

                DONE: begin
                    for (int mm = 0; mm < KERNEL_SIZE; mm++) begin
                        for (int nn = 0; nn < KERNEL_SIZE; nn++) begin
                            dL_dkernel[mm][nn] <= acc[mm][nn][23:8]; // Q16.16 to Q8.8
                        end
                    end
                    done <= 1;
                    state <= IDLE;
                end
            endcase
        end
    end
endmodule

// conv_weight_update.sv
module conv_weight_update #(
    parameter KERNEL_SIZE = 3
) (
    input  logic clk,
    input  logic rst,
    input  logic start,
    input  logic signed [15:0] dL_dkernel  [KERNEL_SIZE][KERNEL_SIZE], // Q8.8
    input  logic signed [15:0] learning_rate,                          // Q8.8
    input  logic signed [15:0] kernel_in [KERNEL_SIZE][KERNEL_SIZE],   // original
    output logic signed [15:0] kernel_out[KERNEL_SIZE][KERNEL_SIZE],   // updated
    output logic done
);

    typedef enum logic [1:0] { IDLE, UPDATE, DONE } state_t;
    state_t state;

    // Loop indices
    logic [1:0] m, n;

    // Precomputed product of learning rate × dL_dkernel
    logic signed [31:0] update;

    always_ff @(posedge clk or posedge rst) begin
        if (rst) begin
            state <= IDLE;
            done  <= 0;
            m     <= 0;
            n     <= 0;
        end else begin
            case (state)
                IDLE: begin
                    done  <= 0;
                    m     <= 0;
                    n     <= 0;
                    state <= start ? UPDATE : IDLE;
                end

                UPDATE: begin
                    update = learning_rate * dL_dkernel[m][n];
                    kernel_out[m][n] <= kernel_in[m][n] - update[23:8];

                    if (n < KERNEL_SIZE - 1) begin
                        n <= n + 1;
                    end else if (m < KERNEL_SIZE - 1) begin
                        n <= 0;
                        m <= m + 1;
                    end else begin
                        state <= DONE;
                    end
                end

                DONE: begin
                    done  <= 1;
                    state <= IDLE;
                end
            endcase
        end
    end
endmodule


// top_training_pipeline.sv
module top_training_pipeline (
    input  logic clk,
    input  logic rst,
    input  logic start,

    // Input feature map
    input  logic signed [15:0] input_feature [4][4],       // Q8.8

    // Kernel (Conv layer)
    input  logic signed [15:0] kernel_in [3][3],            // Q8.8
    output logic signed [15:0] kernel_out [3][3],           // Q8.8

    // Fully Connected (FC) weights and bias
    input  logic signed [15:0] fc_weights_in [3:0],         // Q8.8
    output logic signed [15:0] fc_weights_out [3:0],        // Q8.8
    input  logic signed [15:0] fc_bias_in,                  // Q8.8
    output logic signed [15:0] fc_bias_out,                 // Q8.8

    // Training label and learning rate
    input  logic signed [15:0] label,                       // Q8.8
    input  logic signed [15:0] learning_rate,               // Q8.8

    // Completion flag
    output logic done
);

    // Intermediate signals
    logic signed [15:0] conv_out  [2][2];
    logic signed [15:0] relu_out  [2][2];
    logic signed [15:0] relu_flat [3:0];
    logic signed [15:0] fc_output;
    logic signed [15:0] loss_grad;
    logic signed [15:0] dL_drelu [3:0];
    logic signed [15:0] dL_dconv_out [2][2];
    logic signed [15:0] dL_dkernel [3][3];

    // Done flags
    logic conv_done, fc_done, fc_bp_done, conv_bp_done, conv_upd_done;

    // Updated weights/bias (local, forwarded to top-level outputs)
    logic signed [15:0] updated_fc_weights [3:0];
    logic signed [15:0] updated_fc_bias;
    logic signed [15:0] updated_kernel_weights [3][3];

    // -------------------------------
    // Forward Path
    // -------------------------------
    conv2d_unit_pipelined conv (
        .clk(clk),
        .rst(rst),
        .start(start),
        .input_feature(input_feature),
        .kernel_weights(kernel_in),
        .done(conv_done),
        .output_feature(conv_out)
    );

    relu_layer_2d #(.HEIGHT(2), .WIDTH(2)) relu (
        .in_feature(conv_out),
        .out_feature(relu_out)
    );

    always_comb begin
        relu_flat[0] = relu_out[0][0];
        relu_flat[1] = relu_out[0][1];
        relu_flat[2] = relu_out[1][0];
        relu_flat[3] = relu_out[1][1];
    end

    fc_layer #(.INPUT_DIM(4)) fc (
        .clk(clk),
        .rst(rst),
        .start(conv_done),
        .input_vec(relu_flat),
        .weights(fc_weights_in),
        .bias(fc_bias_in),
        .done(fc_done),
        .output_val(fc_output)
    );

    // -------------------------------
    // Loss Gradient
    // -------------------------------
    loss_gradient loss (
        .prediction(fc_output),
        .label(label),
        .dL_dout(loss_grad)
    );

    // -------------------------------
    // FC Backprop + Weight Update
    // -------------------------------
    fc_backprop #(.INPUT_DIM(4)) fc_bp (
        .clk(clk),
        .rst(rst),
        .start(fc_done),
        .input_vec(relu_flat),
        .dL_dout(loss_grad),
        .learning_rate(learning_rate),
        .weights_in(fc_weights_in),
        .weights_out(updated_fc_weights),
        .bias_in(fc_bias_in),
        .bias_out(updated_fc_bias),
        .dL_drelu(dL_drelu),
        .done(fc_bp_done)
    );

    // -------------------------------
    // ReLU Derivative Backprop
    // -------------------------------
    always_comb begin
        for (int i = 0; i < 4; i++) begin
            dL_dconv_out[i / 2][i % 2] = (relu_flat[i][15] == 1) ? 16'sd0 : dL_drelu[i];
        end
    end

    // -------------------------------
    // Conv Kernel Gradient
    // -------------------------------
    conv_kernel_grad conv_bp (
        .clk(clk),
        .rst(rst),
        .start(fc_bp_done),
        .input_feature(input_feature),
        .dL_dconv_out(dL_dconv_out),
        .dL_dkernel(dL_dkernel),
        .done(conv_bp_done)
    );

    // -------------------------------
    // Conv Weight Update
    // -------------------------------
    conv_weight_update conv_upd (
        .clk(clk),
        .rst(rst),
        .start(conv_bp_done),
        .dL_dkernel(dL_dkernel),
        .learning_rate(learning_rate),
        .kernel_in(kernel_in),
        .kernel_out(updated_kernel_weights),
        .done(conv_upd_done)
    );

    // -------------------------------
    // Final Output Assignments
    // -------------------------------
    assign fc_weights_out = updated_fc_weights;
    assign fc_bias_out    = updated_fc_bias;
    assign kernel_out     = updated_kernel_weights;
    assign done           = conv_upd_done;

endmodule
