`timescale 1ns/1ps

module tb_training_pipeline;

    logic clk = 0;
    logic rst;
    logic start_forward, start_train;

    // Input
    logic signed [15:0] input_feature [4][4];

    // Kernel weights: original and trained
    logic signed [15:0] kernel_in   [3][3];
    logic signed [15:0] kernel_out  [3][3];

    // FC weights & bias: original and trained
    logic signed [15:0] fc_weights_in   [3:0];
    logic signed [15:0] fc_weights_out  [3:0];
    logic signed [15:0] fc_bias_in;
    logic signed [15:0] fc_bias_out;

    // Active weights used for forward pass
    logic signed [15:0] kernel_active   [3][3];
    logic signed [15:0] fc_weights_active [3:0];
    logic signed [15:0] fc_bias_active;

    // Training params
    logic signed [15:0] label;
    logic signed [15:0] learning_rate;

    // Forward output
    logic signed [15:0] forward_output;
    logic forward_done, train_done;

    // Clock generation
    always #5 clk = ~clk;

    // ----------------------------------------
    // Single Forward Pipeline Instance
    // ----------------------------------------
    top_forward_pipeline forward (
        .clk(clk),
        .rst(rst),
        .start(start_forward),
        .input_feature(input_feature),
        .kernel_weights(kernel_active),
        .fc_weights(fc_weights_active),
        .fc_bias(fc_bias_active),
        .done(forward_done),
        .output_value(forward_output)
    );

    // ----------------------------------------
    // Training Pipeline
    // ----------------------------------------
    top_training_pipeline train (
        .clk(clk),
        .rst(rst),
        .start(start_train),
        .input_feature(input_feature),

        .kernel_in(kernel_in),
        .kernel_out(kernel_out),

        .fc_weights_in(fc_weights_in),
        .fc_weights_out(fc_weights_out),

        .fc_bias_in(fc_bias_in),
        .fc_bias_out(fc_bias_out),

        .label(label),
        .learning_rate(learning_rate),
        .done(train_done)
    );

    // ----------------------------------------
    // Test Procedure
    // ----------------------------------------
    initial begin
        // Reset
        rst = 1;
        start_forward = 0;
        start_train = 0;
        #20 rst = 0;
        #10;

        // Input image = all 1s
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                input_feature[i][j] = 16'sd256;

        // Kernel = 0.25
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                kernel_in[i][j] = 16'sd64;

        // FC weights = 0.5
        for (int i = 0; i < 4; i++)
            fc_weights_in[i] = 16'sd128;
        fc_bias_in = 16'sd0;

        label = 16'sd512;         // 2.0 in Q8.8
        learning_rate = 16'sd26;  // ≈ 0.1 in Q8.8

        // -------------------------------------
        // Forward BEFORE Training
        // -------------------------------------
        $display("\n=== Forward Pass (Before Training) ===");
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                kernel_active[i][j] = kernel_in[i][j];
        for (int i = 0; i < 4; i++)
            fc_weights_active[i] = fc_weights_in[i];
        fc_bias_active = fc_bias_in;

        start_forward = 1;
        #10 start_forward = 0;
        wait(forward_done);
        $display("Output Before Training: %0d (≈ %0.2f)", 
                 forward_output, $itor(forward_output)/256.0);

        // -------------------------------------
        // Train
        // -------------------------------------
        $display("\n=== Training Step ===");
        start_train = 1;
        #10 start_train = 0;
        wait(train_done);

        $display("\nUpdated FC Weights:");
        for (int i = 0; i < 4; i++)
            $display("  fc_weights_out[%0d] = %0d", i, fc_weights_out[i]);
        $display("Updated FC Bias: %0d", fc_bias_out);

        $display("Updated Kernel:");
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                $display("  kernel_out[%0d][%0d] = %0d", i, j, kernel_out[i][j]);

        // -------------------------------------
        // Forward AFTER Training
        // -------------------------------------
        $display("\n=== Forward Pass (After Training) ===");
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                kernel_active[i][j] = kernel_out[i][j];
        for (int i = 0; i < 4; i++)
            fc_weights_active[i] = fc_weights_out[i];
        fc_bias_active = fc_bias_out;

        start_forward = 1;
        #10 start_forward = 0;
        wait(forward_done);
        $display("Output After Training: %0d (≈ %0.2f)", 
                 forward_output, $itor(forward_output)/256.0);

        $finish;
    end

endmodule
