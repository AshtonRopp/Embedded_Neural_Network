
`timescale 1ns / 1ps

module cnn_top_modular_tb;

    logic clk, rst, start, done;
    logic signed [15:0] input_image [4][4];
    logic signed [15:0] label;
    logic signed [15:0] learning_rate;
    logic signed [15:0] output_value;

    // Training data set
    logic signed [15:0] images   [2][4][4];
    logic signed [15:0] labels   [2];

    logic signed [15:0] err;
    logic signed [15:0] abs_error;
    assign abs_error = (err < 0) ? -err : err;

    cnn_top_modular uut (
        .clk(clk),
        .rst(rst),
        .start(start),
        .input_image(input_image),
        .label(label),
        .output_value(output_value),
        .done(done)
    );

    // Clock generation
    always #5 clk = ~clk;

    initial begin
        $display("\n=== CNN Modular Testbench: Multi-Epoch Training ===");

        clk = 0;
        rst = 1;
        start = 0;
        learning_rate = 16'sd2; // Q8.8 = 0.007812


        // Input #0: All 1s → label 1.0
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                images[0][i][j] = 16'sd256;
        labels[0] = 16'sd256;

        // Input #1: All 0s → label 0.0
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                images[1][i][j] = 16'sd0;
        labels[1] = 16'sd0;

        #20 rst = 0;

        for (int epoch = 0; epoch < 20; epoch++) begin
            $display("\n=== Epoch %0d ===", epoch);

            for (int sample = 0; sample < 2; sample++) begin
                // Load input and label
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 4; j++)
                        input_image[i][j] = images[sample][i][j];
                label = labels[sample];

                // Start
                #10;
                start = 1;
                #20;
                start = 0;

                // Wait for complete
                wait (done);

                // Print loss
                err = output_value - label;
                $display("Prediction: %0d (~= %.2f) | Label: %0d | Error: %0d",output_value, output_value / 256.0, label, abs_error);
            end
        end

        $display("\n=== Training Complete ===");
        $finish;
    end

endmodule
