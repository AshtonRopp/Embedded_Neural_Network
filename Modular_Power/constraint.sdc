# Define the clock
create_clock -name clk -period 10.0 [get_ports clk]

# Input delay
set_input_delay 2.5 -clock clk [get_ports {conv_weights fc1_bias fc1_weights fc2_bias fc2_weights input_image label learning_rate rst start}]

# Output delay
set_output_delay 2.5 -clock clk [get_ports {done output_value}]