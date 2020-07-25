#include "generator.h"

class multi_generator
{
private:
	uint64_t n;
	uint64_t output_n;
	vector<bool> gamma;
public:
	vector<function> function_to_reg;
	function general_function;
	multi_generator(uint64_t  n, uint64_t  N, vector<vector<bool>> input, vector<bool> input_general) {
		this->n = n;
		output_n = N;
		for (int k = 0; k < input.size(); ++k)
			function_to_reg.push_back(function(n, input[k]));
		general_function = function(N, input_general);
	}
	void generate_sequence(uint64_t first_state) {
		uint64_t cut = (static_cast<uint64_t>(1) << n) - static_cast<uint64_t>(1);
		vector<uint64_t > states;
		states.resize(function_to_reg.size());
		for (int k = 0; k < function_to_reg.size(); ++k)
			states[k] = (first_state & (cut << (k * n))) >> (k * n);
		uint64_t general_argument = 0;
		uint64_t N = (static_cast<uint64_t>(1) << n);
		ofstream out("C:\\Users\\0\\source\\repos\\generator\\out_file.txt", std::ios::app);
		vector<uint64_t> x; x.resize(8);
		int c = 0;
		int print = 8000;
		int k = print;
		for (int k0 = 0; k0 < N; ++k0) {
			for (int k1 = 0; k1 < N; ++k1) {
				for (int k2 = 0; k2 < N; ++k2) {
					for (int k3 = 0; k3 < N; ++k3) {
						for (int k4 = 0; k4 < N; ++k4) {
							for (int k5 = 0; k5 < N; ++k5) {
								for (int k6 = 0; k6 < N; ++k6) {
									for (int k7 = 0; k7 < N; ++k7) {
										if (function_to_reg[0].values[states[0]]) general_argument ^= 1;
										if (function_to_reg[1].values[states[1]]) general_argument ^= (1 << 1);
										if (function_to_reg[2].values[states[2]]) general_argument ^= (1 << 2);
										if (function_to_reg[3].values[states[3]]) general_argument ^= (1 << 3);
										if (function_to_reg[4].values[states[4]]) general_argument ^= (1 << 4);
										if (function_to_reg[5].values[states[5]]) general_argument ^= (1 << 5);
										if (function_to_reg[6].values[states[6]]) general_argument ^= (1 << 6);
										if (function_to_reg[7].values[states[7]]) general_argument ^= (1 << 7);
										//if (general_function.values[general_argument]) gamma.push_back(true);
										//else gamma.push_back(false);
										if (general_function.return_value(general_argument)) out << 1;
										else out << 0;
										switch (c) {
										case 1:
											x[0] = function_to_reg[0].return_value(states[0]);
											states[0] = states[0] << 1;
											states[0] = (states[0] + x[0]) & cut;
										case 7:
											x[1] = function_to_reg[1].return_value(states[1]);
											states[1] = states[1] << 1;
											states[1] = (states[1] + x[1]) & cut;
										case 0:
											x[2] = function_to_reg[2].return_value(states[2]);
											states[2] = states[2] << 1;
											states[2] = (states[2] + x[2]) & cut;
										case 5:
											x[3] = function_to_reg[3].return_value(states[3]);
											states[3] = states[3] << 1;
											states[3] = (states[3] + x[3]) & cut;
										case 2:
											x[4] = function_to_reg[4].return_value(states[4]);
											states[4] = states[4] << 1;
											states[4] = (states[4] + x[4]) & cut;
										case 4:
											x[5] = function_to_reg[5].return_value(states[5]);
											states[5] = states[5] << 1;
											states[5] = (states[5] + x[5]) & cut;
										case 3:
											x[6] = function_to_reg[6].return_value(states[6]);
											states[6] = states[6] << 1;
											states[6] = (states[6] + x[6]) & cut;
										case 6:
											x[7] = function_to_reg[7].return_value(states[7]);
											states[7] = states[7] << 1;
											states[7] = (states[7] + x[7]) & cut;
										}
										if (!k) {
											out << endl;
											k = print;
										}
										++c;
										--k;
										c = c % 8;
										general_argument = 0;
									}
								}
							}
						}
					}
				}
			}
		}
		out.close();
	}
};
