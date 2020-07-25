#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>

using namespace std;

class function
{
private:
	uint64_t n;
	vector<int>  deg_of_function;
	uint64_t value;
public:
	vector<bool> values;
	function() :n(0), value(0) {}
	uint64_t return_value(uint64_t index) { return static_cast<uint64_t>(values[index]); }
	function(uint64_t N, const uint64_t value) :value(0) {
		n = N;
		values.resize(static_cast<uint64_t>(1) << n);
		deg_of_function.resize(static_cast<uint64_t>(1) << n);
		int c = 0;
		for (uint64_t j = 0; j < (static_cast<uint64_t>(1) << n); ++j)
		{
			if ((value >> j) & 1) values[values.size() - 1 - j] = true;
			for (uint64_t k = 0; k < n; ++k)
				if ((j >> k) & 1) ++c;
			deg_of_function[j] = c;
			c = 0;
		}
	}
	function(uint64_t N, const vector<bool>& value) :value(0) {
		n = N;
		values.resize(value.size());
		deg_of_function.resize(value.size());
		int c = 0;
		for (uint64_t j = 0; j < value.size(); ++j) {
			if (value[j]) values[j] = true;
			for (uint64_t k = 0; k < n; ++k)
				if ((j >> k) & 1) ++c;
			deg_of_function[j] = c;
			c = 0;
		}
	}
	float equilibrium() {
		uint64_t one = 0;
		for (uint64_t j = 0; j < (static_cast<uint64_t>(1) << n); ++j)
			if (values[j]) ++one;
		return static_cast<float>(one) / static_cast<float>(static_cast<uint64_t>(1) << n);
	}
	void print_table(std::ofstream& out) {
		string m;
		for (uint64_t j = 0; j < (static_cast<uint64_t>(1) << n); ++j) {
			for (uint64_t i = 0; i < n; ++i)
				m += (((j >> i) & 1) + '0');
			reverse(m.begin(), m.end());
			out << m;
			m = "";
			out << " " << (values[j]) << endl;
		}
	}
	void print_values(std::ofstream& out) {
		for (uint64_t j = 0; j < (static_cast<uint64_t>(1) << n); ++j)
			out << (values[j]);
		out << endl;
	}
	void print_values(string& out) {
		for (uint64_t j = 0; j < (static_cast<uint64_t>(1) << n); ++j)
			out += to_string(values[j]);
		out += '\n';
	}
	void print_anf(std::ofstream& out) {
		if (values[0])
			out << "1+";
		for (uint64_t j = 1; j < static_cast<uint64_t>(1) << n; ++j) {
			if (values[j]) {
				for (uint64_t k = 0; k < n; ++k)
					if ((j >> k) & 1) {
						out << 'X';
						out << (n - (k));
					}
				out << '+';
			}
		}
	}
	void print_anf(string& out) {
		if (values[0])
			out += "1+";
		for (uint64_t j = 1; j < static_cast<uint64_t>(1) << n; ++j) {
			if (values[j]) {
				for (uint64_t k = 0; k < n; ++k)
					if ((j >> k) & 1) {
						out += 'X';
						out += to_string(k + 1);
					}
				out += '+';
			}
		}
	}
	int return_degree() {
		int deg = 0;
		for (uint64_t j = 0; j < (static_cast<uint64_t>(1) << n); ++j)
			if (values[j])
				if (deg < deg_of_function[j]) deg = deg_of_function[j];
		return deg;
	}
	void set_values(uint64_t val) {
		value = val;
		int c = 0;
		for (uint64_t j = 0; j < (uint64_t(1) << n); ++j)
		{
			if ((val >> j) & uint64_t(1)) values[j] = true;
			for (int k = 0; k < (uint64_t(1) << n); ++k)
				if ((j >> k) & uint64_t(1)) ++c;
			deg_of_function[j] = c;
			c = 0;
		}
	}
	void set_n(uint64_t num) {
		n = num;
		values.clear();
		deg_of_function.clear();
		values.resize(uint64_t(1) << n);
		deg_of_function.resize(uint64_t(1) << n);
	}
	void transformation() {
		uint64_t blocks = 1;
		for (int k = (1 << (n - 1)); k >= 1; k /= 2) {
			blocks *= 2;
			for (int i = 0; i < blocks - 1; i += 2)
				for (int j = 0; j < k; ++j)
					if (values[(i + 1) * k + j] != values[i * k + j])
						values[(i + 1) * k + j] = true;
					else
						values[(i + 1) * k + j] = false;
		}
	}
};
