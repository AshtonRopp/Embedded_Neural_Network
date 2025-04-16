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


module mac_unit #(
    parameter MAC_DEPTH = 4
) (
    input  logic                  clk,
    input  logic                  rst,
    input  logic signed [15:0]    a [MAC_DEPTH],  // Q8.8
    input  logic signed [15:0]    b [MAC_DEPTH],  // Q8.8
    output logic signed [15:0]    result          // Q8.8
);

    // Q16.16 products
    logic signed [31:0] products [MAC_DEPTH];

    // Step 1: Compute products
    always_comb begin
        for (int i = 0; i < MAC_DEPTH; i++) begin
            products[i] = a[i] * b[i];
        end
    end

    // Step 2: Tree adder (works for powers of 2)
    localparam LEVELS = $clog2(MAC_DEPTH);
    logic signed [31:0] sums [LEVELS-1:0][MAC_DEPTH-1:0];
    generate
        // First level = N/2, second = N/(2*2), ... , last = N/N
        for (genvar i = 0; i < LEVELS; i ++) begin : ADDER_SUMS
            if (i == 0) begin
                for (genvar j = 0; j < MAC_DEPTH/2; j++) begin : MULT_SUMS
                    assign sums[i][j] = products[2*j] + products[2*j+1];
                end
            end
            else begin
                for (genvar j = 0; j < MAC_DEPTH >> i; j++) begin : INIT_ADDERS
                    assign sums[i][j] = sums[i-1][2*j] + sums[i-1][2*j+1];
                end
            end
        end
    endgenerate

    assign result = (sums[LEVELS-1][0]) >>> 8; // Q16.16 → Q8.8

endmodule


module mac_unit_9 (
    input  logic                  clk,
    input  logic                  rst,
    input  logic signed [15:0]    a [9],  // Q8.8
    input  logic signed [15:0]    b [9],  // Q8.8
    output logic signed [15:0]    result  // Q8.8
);
    localparam MAC_DEPTH = 9;

    // Q16.16 products
    logic signed [31:0] products [MAC_DEPTH-1:0];

    // Step 1: Compute products
    always_comb begin
        for (int i = 0; i < MAC_DEPTH; i++) begin
            products[i] = a[i] * b[i];
        end
    end

    // Step 2: Tree adder (works for powers of 2)
    localparam LEVELS = $clog2(MAC_DEPTH);
    logic signed [31:0] sums [MAC_DEPTH-2:0];
    always_comb begin
        sums[0] = products[0] + products[1];
        sums[1] = products[2] + products[3];
        sums[2] = products[4] + products[5];
        sums[3] = products[6] + products[7];
        sums[4] = sums[0] + sums[1];
        sums[5] = sums[2] + sums[3];
        sums[6] = sums[4] + sums[5];
        sums[7] = sums[6] + products[MAC_DEPTH-1];
    end

    assign result = (sums[7]) >>> 8;

endmodule


module conv2d_unit_pipelined #(
    parameter IN_SIZE     = 4,
    parameter KERNEL_SIZE = 3
) (
    input  logic clk,
    input  logic rst,
    input  logic start,
    input  logic signed [15:0] input_feature [IN_SIZE][IN_SIZE],
    input  logic signed [15:0] kernel_weights [KERNEL_SIZE][KERNEL_SIZE],
    output logic done,
    output logic signed [15:0] output_feature [IN_SIZE-KERNEL_SIZE+1][IN_SIZE-KERNEL_SIZE+1]
);

    localparam OUT_SIZE  = IN_SIZE - KERNEL_SIZE + 1;
    localparam MAC_DEPTH = KERNEL_SIZE * KERNEL_SIZE;

    // Flatten kernel weights once
    logic signed [15:0] b_flat [MAC_DEPTH];
    genvar m, n;
    generate
        for (m = 0; m < KERNEL_SIZE; m++) begin : FLATTEN_M
            for (n = 0; n < KERNEL_SIZE; n++) begin : FLATTEN_N
                localparam int IDX = m * KERNEL_SIZE + n;
                assign b_flat[IDX] = kernel_weights[m][n];
            end
        end
    endgenerate

    // One MAC unit per output pixel
    logic signed [15:0] mac_inputs_a [OUT_SIZE][OUT_SIZE][MAC_DEPTH];
    logic signed [15:0] mac_outputs [OUT_SIZE][OUT_SIZE];

    genvar i, j, k;
    generate
        for (i = 0; i < OUT_SIZE; i++) begin : ROW
            for (j = 0; j < OUT_SIZE; j++) begin : COL
                for (k = 0; k < MAC_DEPTH; k++) begin : PATCH
                    localparam int mi = k / KERNEL_SIZE;
                    localparam int ni = k % KERNEL_SIZE;
                    always_comb begin
                        mac_inputs_a[i][j][k] = input_feature[i + mi][j + ni];
                    end
                end

                mac_unit_9 mac_inst (
                    .clk(clk),
                    .rst(rst),
                    .a(mac_inputs_a[i][j]),
                    .b(b_flat),
                    .result(mac_outputs[i][j])
                );
            end
        end
    endgenerate

    // Assign results
    always_ff @(posedge clk) begin
        if (rst) begin
            done <= 0;
        end else if (start) begin
            for (int i = 0; i < OUT_SIZE; i++) begin
                for (int j = 0; j < OUT_SIZE; j++) begin
                    output_feature[i][j] <= mac_outputs[i][j];
                end
            end
            done <= 1;
        end else begin
            done <= 0;
        end
    end

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

    logic signed [15:0] mac_result;

    always_ff @(posedge clk) begin
        if (rst) begin
            output_val <= 32'b0;
            done <= 0;
        end
        else begin
            if (start) begin
                output_val <= mac_result + bias;
                done <= 1;
            end
            else begin
                output_val <= mac_result + bias;
                done <= 0;
            end
        end
    end

    mac_unit #(.MAC_DEPTH(INPUT_DIM)) mac_inst (
        .clk(clk),
        .rst(rst),
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

    logic signed [31:0] grad_mul     [INPUT_DIM];
    logic signed [31:0] update_mul   [INPUT_DIM];
    logic signed [31:0] backprop_mul [INPUT_DIM];
    logic signed [31:0] bias_update;

    always_comb begin
        for (int i = 0; i < INPUT_DIM; i++) begin
            grad_mul[i]     = input_vec[i] * dL_dout;
            update_mul[i]   = learning_rate * (grad_mul[i] >>> 8);
            backprop_mul[i] = weights_in[i] * dL_dout;
        end
        bias_update = learning_rate * dL_dout;
    end

    always_ff @(posedge clk) begin
        if (rst) begin
            for (int i = 0; i < INPUT_DIM; i++) begin
                weights_out[i] <= 16'd0;
                dL_drelu[i]    <= 16'd0;
                bias_out <= 16'd0;
            end
            done <= 0;
        end
        else begin
            if (start) begin
                for (int i = 0; i < INPUT_DIM; i++) begin
                    weights_out[i] <= weights_in[i] - update_mul[i][23:8];
                    dL_drelu[i]    <= backprop_mul[i][23:8];
                    bias_out <= bias_in - bias_update[23:8];
                end
                done <= 1;
            end
            else begin
                for (int i = 0; i < INPUT_DIM; i++) begin
                    weights_out[i] <= 16'd0;
                    dL_drelu[i]    <= 16'd0;
                    bias_out <= 16'd0;
                end
                done <= 0;
            end
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
