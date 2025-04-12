`timescale 1ns/1ps

module testbench;

    logic clk = 0;
    logic rst;
    logic start_forward, start_train;

    // Input
    logic signed [15:0] input_feature [4][4];

    // Kernel
    logic signed [15:0] kernel_in  [3][3];
    logic signed [15:0] kernel_out [3][3];

    // FC weights & bias
    logic signed [15:0] fc_weights_in  [3:0];
    logic signed [15:0] fc_weights_out [3:0];
    logic signed [15:0] fc_bias_in;
    logic signed [15:0] fc_bias_out;

    // Training parameters
    logic signed [15:0] label;
    logic signed [15:0] learning_rate;

    // Outputs
    logic signed [15:0] output_before;
    logic signed [15:0] output_after;

    logic forward_done_before, forward_done_after, train_done;

    // Clock generation
    always #5 clk = ~clk;

    // --------------------------------
    // Forward path before training
    // --------------------------------
    top_forward_pipeline forward_before (
        .clk(clk),
        .rst(rst),
        .start(start_forward),
        .input_feature(input_feature),
        .kernel_weights(kernel_in),
        .fc_weights(fc_weights_in),
        .fc_bias(fc_bias_in),
        .done(forward_done_before),
        .output_value(output_before)
    );

    // --------------------------------
    // Forward path after training
    // --------------------------------
    top_forward_pipeline forward_after (
        .clk(clk),
        .rst(rst),
        .start(start_forward),
        .input_feature(input_feature),
        .kernel_weights(kernel_out),
        .fc_weights(fc_weights_out),
        .fc_bias(fc_bias_out),
        .done(forward_done_after),
        .output_value(output_after)
    );

    // --------------------------------
    // Training pipeline
    // --------------------------------
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

    // --------------------------------
    // Test procedure
    // --------------------------------
    initial begin
        // Reset
        rst = 1;
        start_forward = 0;
        start_train = 0;
        #20;
        rst = 0;
        #10;

        // Input = all 1s (Q8.8 = 256)
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                input_feature[i][j] = 16'sd256;

        // Kernel weights = 0.25 (Q8.8 = 64)
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                kernel_in[i][j] = 16'sd64;

        // FC weights = 0.5 (Q8.8 = 128)
        for (int i = 0; i < 4; i++)
            fc_weights_in[i] = 16'sd128;

        fc_bias_in = 16'sd0;
        label = 16'sd512;           // 2.0 in Q8.8
        learning_rate = 16'sd26;    // ~0.1 in Q8.8

        // ----------------------------------
        // Inference BEFORE Training
        // ----------------------------------
        $display("\n=== Forward Pass (Before Training) ===");
        start_forward = 1;
        #20;
        start_forward = 0;
        wait(forward_done_before);
        $display("Output Before Training: %0d (≈ %0.2f)", 
                 output_before, $itor(output_before)/256.0);

        // ----------------------------------
        // Run One Training Step
        // ----------------------------------
        $display("\n=== Training Step ===");
        start_train = 1;
        #20;
        start_train = 0;
        wait(train_done);

        // Show updated weights
        $display("\nUpdated FC Weights:");
        for (int i = 0; i < 4; i++)
            $display("  fc_weights_out[%0d] = %0d", i, fc_weights_out[i]);
        $display("Updated FC Bias: %0d", fc_bias_out);

        $display("Updated Kernel:");
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                $display("  kernel_out[%0d][%0d] = %0d", i, j, kernel_out[i][j]);

        // ----------------------------------
        // Inference AFTER Training
        // ----------------------------------
        $display("\n=== Forward Pass (After Training) ===");
        start_forward = 1;
        #20;
        start_forward = 0;
        wait(forward_done_after);
        $display("Output After Training: %0d (≈ %0.2f)", 
                 output_after, $itor(output_after)/256.0);

        $finish;
    end

endmodule
