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
    parameter INPUT_DIM = 4,
    parameter LOG_INPUT_DIM = 2
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
    logic [3:0] idx; // 4-bit counter
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


// cnn_top_modular.sv
// Modular CNN with external weight inputs
module cnn_top_modular (
    input  logic clk,
    input  logic rst,
    input  logic start,

    input  logic signed [15:0] input_image [4][4],   // Q8.8
    input  logic signed [15:0] label,                // Q8.8
    input  logic signed [15:0] learning_rate,        // Q8.8

    input  logic signed [15:0] conv_weights [4][3][3],
    input  logic signed [15:0] fc1_weights  [8][16],
    input  logic signed [15:0] fc1_bias     [8],
    input  logic signed [15:0] fc2_weights  [8],
    input  logic signed [15:0] fc2_bias,

    output logic signed [15:0] output_value,         // Q8.8
    output logic done
);

    localparam CONV = 4;
    localparam CONV_OUT = 2;
    localparam FLAT_SIZE = CONV * CONV_OUT * CONV_OUT;
    localparam FC1 = 8;

    logic signed [15:0] conv_output [CONV][CONV_OUT][CONV_OUT];
    logic conv_done [CONV];
    logic conv_all_done;

    logic signed [15:0] relu_output [CONV][CONV_OUT][CONV_OUT];
    logic signed [15:0] flat_relu [FLAT_SIZE];

    logic fc1_done [FC1];
    logic signed [15:0] fc1_out [FC1];
    logic all_fc1_done;

    logic fc2_done;
    logic signed [15:0] fc2_out;

    logic signed [15:0] loss_grad;
    logic fc2_bp_done;
    logic signed [15:0] dL_drelu [FC1];
    logic signed [15:0] updated_fc2_weights [FC1];
    logic signed [15:0] updated_fc2_bias;
    logic signed [15:0] updated_fc1_weights [FC1][FLAT_SIZE];
    logic signed [15:0] updated_fc1_bias [FC1];
    logic fc1_bp_done [FC1];

    always_comb begin
        conv_all_done = 1;
        for (int i = 0; i < CONV; i++) begin
            if (!conv_done[i])
                conv_all_done = 0;
        end

        all_fc1_done = 1;
        for (int i = 0; i < FC1; i++) begin
            if (!fc1_done[i])
                all_fc1_done = 0;
        end

        done = 1;
        for (int i = 0; i < FC1; i++) begin
            if (!fc1_bp_done[i])
                done = 0;
        end
    end

    assign output_value = fc2_out;

    genvar f;
    generate
        for (f = 0; f < CONV; f++) begin : CONV_INST
            conv2d_unit_pipelined conv_inst (
                .clk(clk),
                .rst(rst),
                .start(start),
                .input_feature(input_image),
                .kernel_weights(conv_weights[f]),
                .done(conv_done[f]),
                .output_feature(conv_output[f])
            );
        end
    endgenerate

    generate
        for (f = 0; f < CONV; f++) begin : RELU_INST
            relu_layer_2d #(.HEIGHT(CONV_OUT), .WIDTH(CONV_OUT)) relu (
                .in_feature(conv_output[f]),
                .out_feature(relu_output[f])
            );
        end
    endgenerate

    always_comb begin
        for (int f = 0; f < CONV; f++) begin
            for (int i = 0; i < CONV_OUT; i++) begin
                for (int j = 0; j < CONV_OUT; j++) begin
                    flat_relu[f * 4 + i * 2 + j] = relu_output[f][i][j];
                end
            end
        end
    end

    genvar n;
    generate
        for (n = 0; n < FC1; n++) begin : FC1_INST
            fc_layer #(.INPUT_DIM(FLAT_SIZE)) fc1 (
                .clk(clk),
                .rst(rst),
                .start(conv_all_done),
                .input_vec(flat_relu),
                .weights(fc1_weights[n]),
                .bias(fc1_bias[n]),
                .done(fc1_done[n]),
                .output_val(fc1_out[n])
            );
        end
    endgenerate

    fc_layer #(.INPUT_DIM(FC1)) fc2 (
        .clk(clk),
        .rst(rst),
        .start(all_fc1_done),
        .input_vec(fc1_out),
        .weights(fc2_weights),
        .bias(fc2_bias),
        .done(fc2_done),
        .output_val(fc2_out)
    );

    loss_gradient loss_inst (
        .prediction(fc2_out),
        .label(label),
        .dL_dout(loss_grad)
    );

    fc_backprop #(.INPUT_DIM(FC1)) fc2_bp (
        .clk(clk),
        .rst(rst),
        .start(fc2_done),
        .input_vec(fc1_out),
        .dL_dout(loss_grad),
        .learning_rate(learning_rate),
        .weights_in(fc2_weights),
        .weights_out(updated_fc2_weights),
        .bias_in(fc2_bias),
        .bias_out(updated_fc2_bias),
        .dL_drelu(dL_drelu),
        .done(fc2_bp_done)
    );

    generate
        for (n = 0; n < FC1; n++) begin : FC1_BP
            fc_backprop #(.INPUT_DIM(FLAT_SIZE)) fc1_bp (
                .clk(clk),
                .rst(rst),
                .start(fc2_bp_done),
                .input_vec(flat_relu),
                .dL_dout(dL_drelu[n]),
                .learning_rate(learning_rate),
                .weights_in(fc1_weights[n]),
                .weights_out(updated_fc1_weights[n]),
                .bias_in(fc1_bias[n]),
                .bias_out(updated_fc1_bias[n]),
                .dL_drelu(), // unused
                .done(fc1_bp_done[n])
            );
        end
    endgenerate

endmodule

