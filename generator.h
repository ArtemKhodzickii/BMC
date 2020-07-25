#pragma once
#include "function.h"
using namespace std;

vector<bool> podvector(vector<bool> source, int index_1, int index_2) {
	vector<bool> k(index_2 - index_1 + 1);
	for (int i = 0; i < index_2 - index_1 + 1; i++)
		k[i] = source[index_1 + i];
	return k;
}

class generator
{
private:
	uint64_t precycle;       //k0 с которого начинаетс€ цикл
	uint64_t initial_state;  //начальное состо€ние
	uint64_t n;              //кол-во УпеременныхФ регистра
	uint64_t current_state;  //текущее состо€ние
	vector<bool> gamma;
public:
	uint64_t length_of_cycle;//длина полученного цикла
	function f1;             //изначальное Ц это анф
	function f2;             //функци€
	generator(uint64_t N, uint64_t Initial_state, uint64_t value_of_function) :n(N), initial_state(Initial_state) {
		f1 = function(n, value_of_function);
		f2 = function(n, value_of_function);
		current_state = initial_state;
		precycle = 0;
		length_of_cycle = 0;
	}
	generator(uint64_t N, uint64_t Initial_state, const vector<bool>& value_of_function) :n(N), initial_state(Initial_state) {
		f1 = function(n, value_of_function);
		f2 = function(n, value_of_function);
		current_state = initial_state;
		precycle = 0;
		length_of_cycle = 0;
	}
	void shift_true_table(uint64_t steps) {
		uint64_t x;
		uint64_t cut = (static_cast<uint64_t>(1) << n) - static_cast<uint64_t>(1);
		for (int j = 0; j < steps; j++) {
			gamma.push_back(static_cast<bool>((current_state >> (n - 1)) & 1));
			//gamma += (current_state >> (n - 1)) & 1;
			//gamma = gamma << 1;
			x = f2.return_value(current_state);
			current_state = current_state << 1;
			current_state = (current_state + x) & cut;
		}
	}
	void transformation() {
		uint64_t blocks = 1;
		for (int k = (1 << (n - 1)); k >= 1; k /= 2) {
			blocks *= 2;
			for (int i = 0; i < blocks - 1; i += 2)
				for (int j = 0; j < k; ++j)
					if (f2.values[(i + 1) * k + j] != f2.values[i * k + j])
						f2.values[(i + 1) * k + j] = true;
					else
						f2.values[(i + 1) * k + j] = false;
		}
	}
	void print_gamma(std::ofstream& out) {
		//for (int j = static_cast<uint64_t>(5) *(static_cast<uint64_t>(1) <<n); j > 0; --j) out << (int)((gamma >> j) & 1);
		for (int j = static_cast<uint64_t>(5) * (static_cast<uint64_t>(1) << n); j > 0; --j) out << (static_cast<int>(gamma[j]));
	}
	void print_gamma(string& out) {
		//for (int j = static_cast<uint64_t>(5) * (static_cast<uint64_t>(1) << n); j > 0; --j) out+=to_string((gamma >> j) & 1);
		for (int j = 0; j < gamma.size(); ++j) out += to_string(static_cast<int>(gamma[j]));
	}
	uint64_t translate(uint64_t index1, uint64_t index2)
	{
		uint64_t ret = static_cast<uint64_t>(0);
		for (int j = gamma.size() - index2; j < gamma.size() - index1; ++j) {
			if (gamma[j]) ret += 1;
			ret = ret << 1;
		}
		return ret;
	}
	void cycle_search(){
		int size = gamma.size();
		uint64_t N = static_cast<uint64_t>(1) << n;
		bool flag = false;
		for (auto e = 1; e <= N; ++e) {
			for (auto k = 2; k < (int)((3*N) / (e)); ++k){
				if (podvector(gamma, size - (e - 1) - 1, size - 1) != (podvector(gamma, (size - k*(e)), size - k * (e)+(e - 1)))){
					flag = false;
					break;
				}
				flag = true;
			}
			if (flag) { length_of_cycle = e; break; }
		}
	}
	uint64_t get_lenght() { return length_of_cycle; }
	uint64_t get_n()      { return n; }
};