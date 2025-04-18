# Define the clock
create_clock -name clk -period 10.0 [get_ports clk]

# Input delay
set_input_delay 2.5 -clock clk [get_ports {input_image label rst start}]

# Output delay
set_output_delay 2.5 -clock clk [get_ports {done output_value}]