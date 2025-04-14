module fixed_point_multiplier (
	a,
	b,
	result
);
	reg _sv2v_0;
	input wire signed [15:0] a;
	input wire signed [15:0] b;
	output reg signed [15:0] result;
	reg signed [31:0] mult_full;
	always @(*) begin
		if (_sv2v_0)
			;
		mult_full = a * b;
		result = mult_full[23:8];
	end
	initial _sv2v_0 = 0;
endmodule
module relu_activation (
	in,
	out
);
	reg _sv2v_0;
	input wire signed [15:0] in;
	output reg signed [15:0] out;
	always @(*) begin
		if (_sv2v_0)
			;
		out = (in[15] == 1 ? 16'sd0 : in);
	end
	initial _sv2v_0 = 0;
endmodule
module mac_unit (
	clk,
	rst,
	enable,
	a,
	b,
	result
);
	reg _sv2v_0;
	parameter MAC_DEPTH = 4;
	input wire clk;
	input wire rst;
	input wire enable;
	input wire signed [(MAC_DEPTH * 16) - 1:0] a;
	input wire signed [(MAC_DEPTH * 16) - 1:0] b;
	output reg signed [15:0] result;
	reg signed [31:0] products [0:MAC_DEPTH - 1];
	always @(*) begin
		if (_sv2v_0)
			;
		begin : sv2v_autoblock_1
			reg signed [31:0] i;
			for (i = 0; i < MAC_DEPTH; i = i + 1)
				products[i] = a[((MAC_DEPTH - 1) - i) * 16+:16] * b[((MAC_DEPTH - 1) - i) * 16+:16];
		end
	end
	localparam LEVELS = $clog2(MAC_DEPTH);
	wire signed [31:0] sums [LEVELS - 1:0][MAC_DEPTH - 1:0];
	genvar _gv_i_1;
	generate
		for (_gv_i_1 = 0; _gv_i_1 < LEVELS; _gv_i_1 = _gv_i_1 + 1) begin : ADDER_SUMS
			localparam i = _gv_i_1;
			if (i == 0) begin : genblk1
				genvar _gv_j_1;
				for (_gv_j_1 = 0; _gv_j_1 < (MAC_DEPTH / 2); _gv_j_1 = _gv_j_1 + 1) begin : MULT_SUMS
					localparam j = _gv_j_1;
					assign sums[i][j] = products[2 * j] + products[(2 * j) + 1];
				end
			end
			else begin : genblk1
				genvar _gv_j_2;
				for (_gv_j_2 = 0; _gv_j_2 < (MAC_DEPTH >> i); _gv_j_2 = _gv_j_2 + 1) begin : INIT_ADDERS
					localparam j = _gv_j_2;
					assign sums[i][j] = sums[i - 1][2 * j] + sums[i - 1][(2 * j) + 1];
				end
			end
		end
	endgenerate
	wire signed [31:0] sum;
	assign sum = sums[LEVELS - 1][0];
	always @(posedge clk or posedge rst)
		if (rst)
			result <= 16'sd0;
		else if (enable)
			result <= sum >>> 8;
	initial _sv2v_0 = 0;
endmodule
module mac_unit_9 (
	clk,
	rst,
	enable,
	a,
	b,
	result
);
	reg _sv2v_0;
	input wire clk;
	input wire rst;
	input wire enable;
	input wire signed [143:0] a;
	input wire signed [143:0] b;
	output reg signed [15:0] result;
	localparam MAC_DEPTH = 9;
	reg signed [31:0] products [8:0];
	always @(*) begin
		if (_sv2v_0)
			;
		begin : sv2v_autoblock_1
			reg signed [31:0] i;
			for (i = 0; i < MAC_DEPTH; i = i + 1)
				products[i] = a[(8 - i) * 16+:16] * b[(8 - i) * 16+:16];
		end
	end
	localparam LEVELS = 4;
	reg signed [31:0] sums [7:0];
	always @(*) begin
		if (_sv2v_0)
			;
		sums[0] = products[0] + products[1];
		sums[1] = products[2] + products[3];
		sums[2] = products[4] + products[5];
		sums[3] = products[6] + products[7];
		sums[4] = sums[0] + sums[1];
		sums[5] = sums[2] + sums[3];
		sums[6] = sums[4] + sums[5];
		sums[7] = sums[6] + products[8];
	end
	wire signed [31:0] sum;
	assign sum = sums[7];
	always @(posedge clk or posedge rst)
		if (rst)
			result <= 16'sd0;
		else if (enable)
			result <= sum >>> 8;
	initial _sv2v_0 = 0;
endmodule
module conv2d_unit_pipelined (
	clk,
	rst,
	start,
	input_feature,
	kernel_weights,
	done,
	output_feature
);
	parameter IN_SIZE = 4;
	parameter KERNEL_SIZE = 3;
	input wire clk;
	input wire rst;
	input wire start;
	input wire signed [((IN_SIZE * IN_SIZE) * 16) - 1:0] input_feature;
	input wire signed [((KERNEL_SIZE * KERNEL_SIZE) * 16) - 1:0] kernel_weights;
	output reg done;
	output reg signed [((0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? ((1 - ((IN_SIZE - KERNEL_SIZE) + 0)) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) + ((((IN_SIZE - KERNEL_SIZE) + 0) + (((IN_SIZE - KERNEL_SIZE) + 0) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0)))) - 1) : ((1 - ((IN_SIZE - KERNEL_SIZE) + 0)) * ((IN_SIZE - KERNEL_SIZE) + 1)) + ((((IN_SIZE - KERNEL_SIZE) + 0) * ((IN_SIZE - KERNEL_SIZE) + 1)) - 1)) : (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (((IN_SIZE - KERNEL_SIZE) + 1) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) + ((IN_SIZE - KERNEL_SIZE) - 1) : (((IN_SIZE - KERNEL_SIZE) + 1) * ((IN_SIZE - KERNEL_SIZE) + 1)) - 1)) >= (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? ((IN_SIZE - KERNEL_SIZE) + 0) + (((IN_SIZE - KERNEL_SIZE) + 0) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) : ((IN_SIZE - KERNEL_SIZE) + 0) * ((IN_SIZE - KERNEL_SIZE) + 1)) : (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (IN_SIZE - KERNEL_SIZE) + 0 : 0)) ? ((((0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? ((1 - ((IN_SIZE - KERNEL_SIZE) + 0)) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) + ((((IN_SIZE - KERNEL_SIZE) + 0) + (((IN_SIZE - KERNEL_SIZE) + 0) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0)))) - 1) : ((1 - ((IN_SIZE - KERNEL_SIZE) + 0)) * ((IN_SIZE - KERNEL_SIZE) + 1)) + ((((IN_SIZE - KERNEL_SIZE) + 0) * ((IN_SIZE - KERNEL_SIZE) + 1)) - 1)) : (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (((IN_SIZE - KERNEL_SIZE) + 1) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) + ((IN_SIZE - KERNEL_SIZE) - 1) : (((IN_SIZE - KERNEL_SIZE) + 1) * ((IN_SIZE - KERNEL_SIZE) + 1)) - 1)) - (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? ((IN_SIZE - KERNEL_SIZE) + 0) + (((IN_SIZE - KERNEL_SIZE) + 0) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) : ((IN_SIZE - KERNEL_SIZE) + 0) * ((IN_SIZE - KERNEL_SIZE) + 1)) : (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (IN_SIZE - KERNEL_SIZE) + 0 : 0))) + 1) * 16) + (((0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? ((IN_SIZE - KERNEL_SIZE) + 0) + (((IN_SIZE - KERNEL_SIZE) + 0) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) : ((IN_SIZE - KERNEL_SIZE) + 0) * ((IN_SIZE - KERNEL_SIZE) + 1)) : (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (IN_SIZE - KERNEL_SIZE) + 0 : 0)) * 16) - 1) : ((((0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? ((IN_SIZE - KERNEL_SIZE) + 0) + (((IN_SIZE - KERNEL_SIZE) + 0) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) : ((IN_SIZE - KERNEL_SIZE) + 0) * ((IN_SIZE - KERNEL_SIZE) + 1)) : (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (IN_SIZE - KERNEL_SIZE) + 0 : 0)) - (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? ((1 - ((IN_SIZE - KERNEL_SIZE) + 0)) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) + ((((IN_SIZE - KERNEL_SIZE) + 0) + (((IN_SIZE - KERNEL_SIZE) + 0) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0)))) - 1) : ((1 - ((IN_SIZE - KERNEL_SIZE) + 0)) * ((IN_SIZE - KERNEL_SIZE) + 1)) + ((((IN_SIZE - KERNEL_SIZE) + 0) * ((IN_SIZE - KERNEL_SIZE) + 1)) - 1)) : (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (((IN_SIZE - KERNEL_SIZE) + 1) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) + ((IN_SIZE - KERNEL_SIZE) - 1) : (((IN_SIZE - KERNEL_SIZE) + 1) * ((IN_SIZE - KERNEL_SIZE) + 1)) - 1))) + 1) * 16) + (((0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? ((1 - ((IN_SIZE - KERNEL_SIZE) + 0)) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) + ((((IN_SIZE - KERNEL_SIZE) + 0) + (((IN_SIZE - KERNEL_SIZE) + 0) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0)))) - 1) : ((1 - ((IN_SIZE - KERNEL_SIZE) + 0)) * ((IN_SIZE - KERNEL_SIZE) + 1)) + ((((IN_SIZE - KERNEL_SIZE) + 0) * ((IN_SIZE - KERNEL_SIZE) + 1)) - 1)) : (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (((IN_SIZE - KERNEL_SIZE) + 1) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) + ((IN_SIZE - KERNEL_SIZE) - 1) : (((IN_SIZE - KERNEL_SIZE) + 1) * ((IN_SIZE - KERNEL_SIZE) + 1)) - 1)) * 16) - 1)):((0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? ((1 - ((IN_SIZE - KERNEL_SIZE) + 0)) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) + ((((IN_SIZE - KERNEL_SIZE) + 0) + (((IN_SIZE - KERNEL_SIZE) + 0) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0)))) - 1) : ((1 - ((IN_SIZE - KERNEL_SIZE) + 0)) * ((IN_SIZE - KERNEL_SIZE) + 1)) + ((((IN_SIZE - KERNEL_SIZE) + 0) * ((IN_SIZE - KERNEL_SIZE) + 1)) - 1)) : (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (((IN_SIZE - KERNEL_SIZE) + 1) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) + ((IN_SIZE - KERNEL_SIZE) - 1) : (((IN_SIZE - KERNEL_SIZE) + 1) * ((IN_SIZE - KERNEL_SIZE) + 1)) - 1)) >= (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? ((IN_SIZE - KERNEL_SIZE) + 0) + (((IN_SIZE - KERNEL_SIZE) + 0) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) : ((IN_SIZE - KERNEL_SIZE) + 0) * ((IN_SIZE - KERNEL_SIZE) + 1)) : (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (IN_SIZE - KERNEL_SIZE) + 0 : 0)) ? (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? ((IN_SIZE - KERNEL_SIZE) + 0) + (((IN_SIZE - KERNEL_SIZE) + 0) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) : ((IN_SIZE - KERNEL_SIZE) + 0) * ((IN_SIZE - KERNEL_SIZE) + 1)) : (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (IN_SIZE - KERNEL_SIZE) + 0 : 0)) * 16 : (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? ((1 - ((IN_SIZE - KERNEL_SIZE) + 0)) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) + ((((IN_SIZE - KERNEL_SIZE) + 0) + (((IN_SIZE - KERNEL_SIZE) + 0) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0)))) - 1) : ((1 - ((IN_SIZE - KERNEL_SIZE) + 0)) * ((IN_SIZE - KERNEL_SIZE) + 1)) + ((((IN_SIZE - KERNEL_SIZE) + 0) * ((IN_SIZE - KERNEL_SIZE) + 1)) - 1)) : (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (((IN_SIZE - KERNEL_SIZE) + 1) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) + ((IN_SIZE - KERNEL_SIZE) - 1) : (((IN_SIZE - KERNEL_SIZE) + 1) * ((IN_SIZE - KERNEL_SIZE) + 1)) - 1)) * 16)] output_feature;
	localparam OUT_SIZE = (IN_SIZE - KERNEL_SIZE) + 1;
	localparam MAC_DEPTH = KERNEL_SIZE * KERNEL_SIZE;
	reg signed [(MAC_DEPTH * 16) - 1:0] a_flat;
	wire signed [(MAC_DEPTH * 16) - 1:0] b_flat;
	wire signed [15:0] mac_result;
	reg mac_enable;
	localparam KERNEL_AREA = KERNEL_SIZE * KERNEL_SIZE;
	genvar _gv_m_1;
	genvar _gv_n_1;
	generate
		for (_gv_m_1 = 0; _gv_m_1 < KERNEL_SIZE; _gv_m_1 = _gv_m_1 + 1) begin : FLATTEN_M
			localparam m = _gv_m_1;
			for (_gv_n_1 = 0; _gv_n_1 < KERNEL_SIZE; _gv_n_1 = _gv_n_1 + 1) begin : FLATTEN_N
				localparam n = _gv_n_1;
				localparam signed [31:0] IDX = (m * KERNEL_SIZE) + n;
				assign b_flat[((MAC_DEPTH - 1) - IDX) * 16+:16] = kernel_weights[((((KERNEL_SIZE - 1) - m) * KERNEL_SIZE) + ((KERNEL_SIZE - 1) - n)) * 16+:16];
			end
		end
	endgenerate
	reg [1:0] state;
	reg signed [31:0] x_idx;
	reg signed [31:0] y_idx;
	always @(posedge clk or posedge rst)
		if (rst) begin
			state <= 2'd0;
			done <= 0;
			x_idx <= 0;
			y_idx <= 0;
			mac_enable <= 0;
		end
		else
			case (state)
				2'd0: begin
					done <= 0;
					if (start) begin
						x_idx <= 0;
						y_idx <= 0;
						state <= 2'd1;
					end
				end
				2'd1: begin
					begin : sv2v_autoblock_1
						reg signed [31:0] m;
						for (m = 0; m < KERNEL_SIZE; m = m + 1)
							begin : sv2v_autoblock_2
								reg signed [31:0] n;
								for (n = 0; n < KERNEL_SIZE; n = n + 1)
									a_flat[((MAC_DEPTH - 1) - ((m * KERNEL_SIZE) + n)) * 16+:16] = input_feature[((((IN_SIZE - 1) - (x_idx + m)) * IN_SIZE) + ((IN_SIZE - 1) - (y_idx + n))) * 16+:16];
							end
					end
					mac_enable <= 1;
					state <= 2'd2;
				end
				2'd2: begin
					mac_enable <= 0;
					state <= 2'd3;
				end
				2'd3: begin
					output_feature[((0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? ((1 - ((IN_SIZE - KERNEL_SIZE) + 0)) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) + ((((IN_SIZE - KERNEL_SIZE) + 0) + (((IN_SIZE - KERNEL_SIZE) + 0) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0)))) - 1) : ((1 - ((IN_SIZE - KERNEL_SIZE) + 0)) * ((IN_SIZE - KERNEL_SIZE) + 1)) + ((((IN_SIZE - KERNEL_SIZE) + 0) * ((IN_SIZE - KERNEL_SIZE) + 1)) - 1)) : (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (((IN_SIZE - KERNEL_SIZE) + 1) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) + ((IN_SIZE - KERNEL_SIZE) - 1) : (((IN_SIZE - KERNEL_SIZE) + 1) * ((IN_SIZE - KERNEL_SIZE) + 1)) - 1)) >= (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? ((IN_SIZE - KERNEL_SIZE) + 0) + (((IN_SIZE - KERNEL_SIZE) + 0) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) : ((IN_SIZE - KERNEL_SIZE) + 0) * ((IN_SIZE - KERNEL_SIZE) + 1)) : (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (IN_SIZE - KERNEL_SIZE) + 0 : 0)) ? ((0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? x_idx : ((IN_SIZE - KERNEL_SIZE) + 0) - x_idx) * (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? 1 - ((IN_SIZE - KERNEL_SIZE) + 0) : (IN_SIZE - KERNEL_SIZE) + 1)) + (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? y_idx : ((IN_SIZE - KERNEL_SIZE) + 0) - y_idx) : (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? ((IN_SIZE - KERNEL_SIZE) + 0) + (((IN_SIZE - KERNEL_SIZE) + 0) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) : ((IN_SIZE - KERNEL_SIZE) + 0) * ((IN_SIZE - KERNEL_SIZE) + 1)) : (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (IN_SIZE - KERNEL_SIZE) + 0 : 0)) - ((((0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? x_idx : ((IN_SIZE - KERNEL_SIZE) + 0) - x_idx) * (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? 1 - ((IN_SIZE - KERNEL_SIZE) + 0) : (IN_SIZE - KERNEL_SIZE) + 1)) + (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? y_idx : ((IN_SIZE - KERNEL_SIZE) + 0) - y_idx)) - (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? ((1 - ((IN_SIZE - KERNEL_SIZE) + 0)) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) + ((((IN_SIZE - KERNEL_SIZE) + 0) + (((IN_SIZE - KERNEL_SIZE) + 0) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0)))) - 1) : ((1 - ((IN_SIZE - KERNEL_SIZE) + 0)) * ((IN_SIZE - KERNEL_SIZE) + 1)) + ((((IN_SIZE - KERNEL_SIZE) + 0) * ((IN_SIZE - KERNEL_SIZE) + 1)) - 1)) : (0 >= ((IN_SIZE - KERNEL_SIZE) + 0) ? (((IN_SIZE - KERNEL_SIZE) + 1) * (1 - ((IN_SIZE - KERNEL_SIZE) + 0))) + ((IN_SIZE - KERNEL_SIZE) - 1) : (((IN_SIZE - KERNEL_SIZE) + 1) * ((IN_SIZE - KERNEL_SIZE) + 1)) - 1)))) * 16+:16] <= mac_result;
					if (y_idx < (OUT_SIZE - 1)) begin
						y_idx <= y_idx + 1;
						state <= 2'd1;
					end
					else if (x_idx < (OUT_SIZE - 1)) begin
						x_idx <= x_idx + 1;
						y_idx <= 0;
						state <= 2'd1;
					end
					else begin
						state <= 2'd0;
						done <= 1;
					end
				end
			endcase
	mac_unit_9 mac_inst(
		.clk(clk),
		.rst(rst),
		.enable(mac_enable),
		.a(a_flat),
		.b(b_flat),
		.result(mac_result)
	);
endmodule
module relu_layer_2d (
	in_feature,
	out_feature
);
	reg _sv2v_0;
	parameter HEIGHT = 2;
	parameter WIDTH = 2;
	input wire signed [((HEIGHT * WIDTH) * 16) - 1:0] in_feature;
	output reg signed [((HEIGHT * WIDTH) * 16) - 1:0] out_feature;
	always @(*) begin
		if (_sv2v_0)
			;
		begin : sv2v_autoblock_1
			reg signed [31:0] i;
			for (i = 0; i < HEIGHT; i = i + 1)
				begin : sv2v_autoblock_2
					reg signed [31:0] j;
					for (j = 0; j < WIDTH; j = j + 1)
						out_feature[((((HEIGHT - 1) - i) * WIDTH) + ((WIDTH - 1) - j)) * 16+:16] = (in_feature[(((((HEIGHT - 1) - i) * WIDTH) + ((WIDTH - 1) - j)) * 16) + 15] == 1 ? 16'sd0 : in_feature[((((HEIGHT - 1) - i) * WIDTH) + ((WIDTH - 1) - j)) * 16+:16]);
				end
		end
	end
	initial _sv2v_0 = 0;
endmodule
module fc_layer (
	clk,
	rst,
	start,
	input_vec,
	weights,
	bias,
	done,
	output_val
);
	parameter INPUT_DIM = 4;
	input wire clk;
	input wire rst;
	input wire start;
	input wire signed [(INPUT_DIM * 16) - 1:0] input_vec;
	input wire signed [(INPUT_DIM * 16) - 1:0] weights;
	input wire signed [15:0] bias;
	output reg done;
	output reg signed [15:0] output_val;
	reg mac_enable;
	wire signed [15:0] mac_result;
	reg [1:0] state;
	always @(posedge clk or posedge rst)
		if (rst) begin
			state <= 2'd0;
			mac_enable <= 0;
			done <= 0;
		end
		else
			case (state)
				2'd0: begin
					done <= 0;
					if (start) begin
						mac_enable <= 1;
						state <= 2'd1;
					end
				end
				2'd1: begin
					mac_enable <= 0;
					state <= 2'd2;
				end
				2'd2: begin
					output_val <= mac_result + bias;
					state <= 2'd3;
				end
				2'd3: begin
					done <= 1;
					state <= 2'd0;
				end
			endcase
	mac_unit #(.MAC_DEPTH(INPUT_DIM)) mac_inst(
		.clk(clk),
		.rst(rst),
		.enable(mac_enable),
		.a(input_vec),
		.b(weights),
		.result(mac_result)
	);
endmodule
module loss_gradient (
	prediction,
	label,
	dL_dout
);
	input wire signed [15:0] prediction;
	input wire signed [15:0] label;
	output wire signed [15:0] dL_dout;
	assign dL_dout = prediction - label;
endmodule
module fc_backprop (
	clk,
	rst,
	start,
	input_vec,
	dL_dout,
	learning_rate,
	weights_in,
	weights_out,
	bias_in,
	bias_out,
	dL_drelu,
	done
);
	reg _sv2v_0;
	parameter INPUT_DIM = 4;
	parameter LOG_INPUT_DIM = 2;
	input wire clk;
	input wire rst;
	input wire start;
	input wire signed [(INPUT_DIM * 16) - 1:0] input_vec;
	input wire signed [15:0] dL_dout;
	input wire signed [15:0] learning_rate;
	input wire signed [(INPUT_DIM * 16) - 1:0] weights_in;
	output reg signed [(INPUT_DIM * 16) - 1:0] weights_out;
	input wire signed [15:0] bias_in;
	output reg signed [15:0] bias_out;
	output reg signed [(INPUT_DIM * 16) - 1:0] dL_drelu;
	output reg done;
	reg [1:0] state;
	reg [3:0] idx;
	reg update_bias_done;
	reg signed [31:0] grad_mul [0:INPUT_DIM - 1];
	reg signed [31:0] update_mul [0:INPUT_DIM - 1];
	reg signed [31:0] backprop_mul [0:INPUT_DIM - 1];
	reg signed [31:0] bias_update;
	always @(*) begin
		if (_sv2v_0)
			;
		begin : sv2v_autoblock_1
			reg signed [31:0] i;
			for (i = 0; i < INPUT_DIM; i = i + 1)
				begin
					grad_mul[i] = input_vec[((INPUT_DIM - 1) - i) * 16+:16] * dL_dout;
					update_mul[i] = learning_rate * (grad_mul[i] >>> 8);
					backprop_mul[i] = weights_in[((INPUT_DIM - 1) - i) * 16+:16] * dL_dout;
				end
		end
		bias_update = learning_rate * dL_dout;
	end
	always @(posedge clk or posedge rst)
		if (rst) begin
			state <= 2'd0;
			idx <= 0;
			done <= 0;
			update_bias_done <= 0;
		end
		else
			case (state)
				2'd0: begin
					done <= 0;
					idx <= 0;
					update_bias_done <= 0;
					state <= (start ? 2'd1 : 2'd0);
				end
				2'd1: begin
					weights_out[((INPUT_DIM - 1) - idx) * 16+:16] <= weights_in[((INPUT_DIM - 1) - idx) * 16+:16] - update_mul[idx][23:8];
					dL_drelu[((INPUT_DIM - 1) - idx) * 16+:16] <= backprop_mul[idx][23:8];
					if (idx == (INPUT_DIM - 1)) begin
						bias_out <= bias_in - bias_update[23:8];
						state <= 2'd2;
					end
					else
						idx <= idx + 1;
				end
				2'd2: begin
					done <= 1;
					state <= 2'd0;
				end
			endcase
	initial _sv2v_0 = 0;
endmodule
module cnn_top_modular (
	clk,
	rst,
	start,
	input_image,
	label,
	learning_rate,
	conv_weights,
	fc1_weights,
	fc1_bias,
	fc2_weights,
	fc2_bias,
	output_value,
	done
);
	reg _sv2v_0;
	input wire clk;
	input wire rst;
	input wire start;
	input wire signed [255:0] input_image;
	input wire signed [15:0] label;
	input wire signed [15:0] learning_rate;
	input wire signed [575:0] conv_weights;
	input wire signed [2047:0] fc1_weights;
	input wire signed [127:0] fc1_bias;
	input wire signed [127:0] fc2_weights;
	input wire signed [15:0] fc2_bias;
	output wire signed [15:0] output_value;
	output reg done;
	localparam CONV = 4;
	localparam CONV_OUT = 2;
	localparam FLAT_SIZE = 16;
	localparam FC1 = 8;
	wire signed [63:0] conv_output [0:3];
	wire conv_done [0:3];
	reg conv_all_done;
	wire signed [63:0] relu_output [0:3];
	reg signed [255:0] flat_relu;
	wire fc1_done [0:7];
	wire signed [127:0] fc1_out;
	reg all_fc1_done;
	wire fc2_done;
	wire signed [15:0] fc2_out;
	wire signed [15:0] loss_grad;
	wire fc2_bp_done;
	wire signed [127:0] dL_drelu;
	wire signed [127:0] updated_fc2_weights;
	wire signed [15:0] updated_fc2_bias;
	wire signed [255:0] updated_fc1_weights [0:7];
	wire signed [15:0] updated_fc1_bias [0:7];
	wire fc1_bp_done [0:7];
	always @(*) begin
		if (_sv2v_0)
			;
		conv_all_done = 1;
		begin : sv2v_autoblock_1
			reg signed [31:0] i;
			for (i = 0; i < CONV; i = i + 1)
				if (!conv_done[i])
					conv_all_done = 0;
		end
		all_fc1_done = 1;
		begin : sv2v_autoblock_2
			reg signed [31:0] i;
			for (i = 0; i < FC1; i = i + 1)
				if (!fc1_done[i])
					all_fc1_done = 0;
		end
		done = 1;
		begin : sv2v_autoblock_3
			reg signed [31:0] i;
			for (i = 0; i < FC1; i = i + 1)
				if (!fc1_bp_done[i])
					done = 0;
		end
	end
	assign output_value = fc2_out;
	genvar _gv_f_1;
	generate
		for (_gv_f_1 = 0; _gv_f_1 < CONV; _gv_f_1 = _gv_f_1 + 1) begin : CONV_INST
			localparam f = _gv_f_1;
			conv2d_unit_pipelined conv_inst(
				.clk(clk),
				.rst(rst),
				.start(start),
				.input_feature(input_image),
				.kernel_weights(conv_weights[16 * (3 * ((3 - f) * 3))+:144]),
				.done(conv_done[f]),
				.output_feature(conv_output[f])
			);
		end
		for (_gv_f_1 = 0; _gv_f_1 < CONV; _gv_f_1 = _gv_f_1 + 1) begin : RELU_INST
			localparam f = _gv_f_1;
			relu_layer_2d #(
				.HEIGHT(CONV_OUT),
				.WIDTH(CONV_OUT)
			) relu(
				.in_feature(conv_output[f]),
				.out_feature(relu_output[f])
			);
		end
	endgenerate
	always @(*) begin
		if (_sv2v_0)
			;
		begin : sv2v_autoblock_4
			reg signed [31:0] f;
			for (f = 0; f < CONV; f = f + 1)
				begin : sv2v_autoblock_5
					reg signed [31:0] i;
					for (i = 0; i < CONV_OUT; i = i + 1)
						begin : sv2v_autoblock_6
							reg signed [31:0] j;
							for (j = 0; j < CONV_OUT; j = j + 1)
								flat_relu[(15 - (((f * 4) + (i * 2)) + j)) * 16+:16] = relu_output[f][(((1 - i) * 2) + (1 - j)) * 16+:16];
						end
				end
		end
	end
	genvar _gv_n_2;
	generate
		for (_gv_n_2 = 0; _gv_n_2 < FC1; _gv_n_2 = _gv_n_2 + 1) begin : FC1_INST
			localparam n = _gv_n_2;
			fc_layer #(.INPUT_DIM(FLAT_SIZE)) fc1(
				.clk(clk),
				.rst(rst),
				.start(conv_all_done),
				.input_vec(flat_relu),
				.weights(fc1_weights[16 * ((7 - n) * 16)+:256]),
				.bias(fc1_bias[(7 - n) * 16+:16]),
				.done(fc1_done[n]),
				.output_val(fc1_out[(7 - n) * 16+:16])
			);
		end
	endgenerate
	fc_layer #(.INPUT_DIM(FC1)) fc2(
		.clk(clk),
		.rst(rst),
		.start(all_fc1_done),
		.input_vec(fc1_out),
		.weights(fc2_weights),
		.bias(fc2_bias),
		.done(fc2_done),
		.output_val(fc2_out)
	);
	loss_gradient loss_inst(
		.prediction(fc2_out),
		.label(label),
		.dL_dout(loss_grad)
	);
	fc_backprop #(.INPUT_DIM(FC1)) fc2_bp(
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
		for (_gv_n_2 = 0; _gv_n_2 < FC1; _gv_n_2 = _gv_n_2 + 1) begin : FC1_BP
			localparam n = _gv_n_2;
			fc_backprop #(.INPUT_DIM(FLAT_SIZE)) fc1_bp(
				.clk(clk),
				.rst(rst),
				.start(fc2_bp_done),
				.input_vec(flat_relu),
				.dL_dout(dL_drelu[(7 - n) * 16+:16]),
				.learning_rate(learning_rate),
				.weights_in(fc1_weights[16 * ((7 - n) * 16)+:256]),
				.weights_out(updated_fc1_weights[n]),
				.bias_in(fc1_bias[(7 - n) * 16+:16]),
				.bias_out(updated_fc1_bias[n]),
				.dL_drelu(),
				.done(fc1_bp_done[n])
			);
		end
	endgenerate
	initial _sv2v_0 = 0;
endmodule
