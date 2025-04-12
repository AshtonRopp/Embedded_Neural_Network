
`timescale 1ns / 1ps

module cnn_top_modular_tb;

    logic clk, rst, start, done;
    logic signed [15:0] input_image [4][4];
    logic signed [15:0] label;
    logic signed [15:0] learning_rate;
    logic signed [15:0] output_value;

    // External weight/bias storage
    logic signed [15:0] conv_weights [4][3][3];
    logic signed [15:0] fc1_weights [8][16];
    logic signed [15:0] fc1_bias [8];
    logic signed [15:0] fc2_weights [8];
    logic signed [15:0] fc2_bias;

    // Training data set
    logic signed [15:0] images   [2][4][4];
    logic signed [15:0] labels   [2];

    logic signed [15:0] err;

    cnn_top_modular uut (
        .clk(clk),
        .rst(rst),
        .start(start),
        .input_image(input_image),
        .label(label),
        .learning_rate(learning_rate),
        .output_value(output_value),
        .done(done),
        .conv_weights(conv_weights),
        .fc1_weights(fc1_weights),
        .fc1_bias(fc1_bias),
        .fc2_weights(fc2_weights),
        .fc2_bias(fc2_bias)
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

        // Initialize weights small
        for (int f = 0; f < 4; f++)
            for (int m = 0; m < 3; m++)
                for (int n = 0; n < 3; n++)
                    conv_weights[f][m][n] = 16'sd16; // 0.0625

        for (int o = 0; o < 8; o++) begin
            for (int i = 0; i < 16; i++)
                fc1_weights[o][i] = 16'sd8; // 0.03125
            fc1_bias[o] = 0;
        end

        for (int i = 0; i < 8; i++)
            fc2_weights[i] = 16'sd16; // 0.0625
        fc2_bias = 0;

        #20 rst = 0;

        for (int epoch = 0; epoch < 100; epoch++) begin
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

                // Copy updated FC2 weights/bias
                for (int i = 0; i < 8; i++)
                    fc2_weights[i] = uut.updated_fc2_weights[i];
                fc2_bias = uut.updated_fc2_bias;

                // Copy updated FC1 weights/bias
                for (int o = 0; o < 8; o++) begin
                    for (int i = 0; i < 16; i++)
                        fc1_weights[o][i] = uut.updated_fc1_weights[o][i];
                    fc1_bias[o] = uut.updated_fc1_bias[o];
                end

                // Debug
                $display("--- Sample %0d ---", sample);
                for (int i = 0; i < 8; i++) begin
                    $display("fc1_out[%0d] = %0d", i, uut.fc1_out[i]);
                end
                $display("fc2_out = %0d (≈ %0.2f)", uut.fc2_out, uut.fc2_out / 256.0);

                // Print loss

                err = output_value - label;
                $display("Prediction: %0d (≈ %.2f) | Label: %0d | Error: %0d",output_value, output_value / 256.0, label, err);
            end
        end

        $display("\n=== Training Complete ===");
        $finish;
    end

endmodule
