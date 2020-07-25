#include "generator.h"
#include "multi_generator.h"
#include <sstream>
#include <set>

void transformation(vector<bool>& values, int n) {
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

int coeff_nonlinear(vector<bool>& values, int n) {
	int blocks = 1;
	int tmp;
	vector<int> value; value.resize(values.size());
	for (int j = 0; j < values.size(); ++j)
		if (values[j]) value[j] = -1;
		else value[j] = 1;
	for (int k = (1 << (n - 1)); k >= 1; k /= 2) {
		blocks *= 2;
		for (int i = 0; i < blocks - 1; i += 2)
			for (int j = 0; j < k; ++j){
				tmp = value[i * k + j];
				value[i * k + j]       = tmp - value[(i + 1) * k + j];
				value[(i + 1) * k + j] = tmp + value[(i + 1) * k + j];
			}
	}
	tmp = value[0];
	for (int j = 1; j < values.size(); ++j)
		if (abs(value[j]) > tmp) tmp = abs(value[j]);
	return tmp;
}

bool is_equilibrium(int n, uint64_t value_of_function) {
	uint64_t sum = 0;
	for (int k = 0; k < (1 << n); ++k) {
		sum += value_of_function & 1;
		value_of_function = value_of_function >> 1;
	}
	if (sum != ((1 << n) / 2)) return false;
	else return true;
}

float is_equilibrium(vector<bool> source) {
	int k = 0;
	for (int i = 0; i < source.size(); ++i)
		if (source[i] == true) ++k;
	return static_cast<float>(k) / static_cast<float>(source.size());
}

int weight(vector<bool>& values) {
	int weight = 0;
	for (int k = 0; k < values.size(); ++k)
		if (values[k]) ++weight;
	return weight;
}

bool is_BMC(generator& gen, int n) {
	//deg
	if (gen.f1.return_degree() < 2) return false;
	gen.transformation();
	//anf
	if(!gen.f1.values[0]) return false;
	if (weight(gen.f1.values) % 2) return false;
	//balance
	if (is_equilibrium(gen.f2.values) != 0.5) return false;
	//dual
	for (int k = 0; k < ((1 << n) / 2); ++k)
		if (gen.f2.values[k] == gen.f2.values[k + (1 << (n - 1))]) return false;
	return true;
}

void increment(vector<bool>& func, int N) {
	int k = func.size() - N - 2;
	while (func[k] == true) --k;
	func[k] = true;
	for (int j = k + 1; j < func.size() - N - 1; ++j) func[j] = false;
}

void enumeration(uint64_t n) {
	uint64_t N = (static_cast<uint64_t>(1) << n);
	uint64_t N_ = static_cast<uint64_t>(1) << (n - 1);
	string teta; teta = "C:\\Users\\0\\source\\repos\\generator\\length_of_cycle_";
	int index_of_reg = 0;
	vector<string> out_file; out_file.resize(N);
	for (int k = 0; k < N; ++k) {
		out_file[k] = teta;
		out_file[k] += to_string(k + 1);
		out_file[k] += ".txt";
	}
	teta = "";
	vector<bool> func; func.resize(N); 
	func[N_-1] = true; func[N_] = true; func[0] = true;
	while (weight(func) != (N_+1)) {
		increment(func, N_);
		for (uint64_t value_of_register = 0; value_of_register < N; ++value_of_register) {
			generator gen(n, value_of_register, func);
			if (!is_BMC(gen, n)) break;
			gen.shift_true_table(static_cast<uint64_t>(5) * N);
			gen.cycle_search();
			if (index_of_reg < (gen.length_of_cycle - 1)) index_of_reg = gen.length_of_cycle - 1;
			if (!value_of_register) {
				teta = "function in anf: ";
				gen.f1.print_anf(teta);
				teta += " and true table: ";
				gen.f2.print_values(teta);
			}
			//gen.print_gamma(teta);
			//teta += '\n';
		}
		ofstream out(out_file[index_of_reg], std::ios::app);
		out << teta;
		out.close();
		teta = "";
		index_of_reg = 0;
	}
}

vector<string> find(ifstream& inp, int constanta) {
	vector<string> out;
	string teta;
	int count = 0;
	while (!inp.eof()) {
		getline(inp, teta);
		if (!(count % constanta)) out.push_back(teta);
		++count;
	}
	out.pop_back();
	return out;
}

void find_anf(vector<string>& vect) {
	string str = "function in anf: ";
	int j = str.size();
	str = "";
	for (int index = 0; index < vect.size(); ++index) {
		int t = j;
		while ((vect[index])[t] != ' ') {
			str += (vect[index])[t];
			++t;
		}
		vect[index] = str;
		str = "";
	}
}

void find_tt(vector<string>& vect, int n) {
	string str = "";
	for (int index = 0; index < vect.size(); ++index) {
		reverse(vect[index].begin(), vect[index].begin() + vect[index].size());
		for (int t = 0; t < (1 << n); ++t)
			str += (vect[index])[t];
		reverse(str.begin(), str.begin() + str.size());
		vect[index] = str;
		str = "";
	}
}

void string_to_function_vect_bool(string& data, vector<bool>& vect) {
	for (int j = 0; j < data.size(); ++j)
		if (data[j] == '1') vect.push_back(true);
		else vect.push_back(false);
}

void function_vect_bool_to_string(string& data, vector<bool>& vect) {
	for (int j = 0; j < vect.size(); ++j)
		if (vect[j]) data[j] = '1';
		else data[j] = '0';
}

uint64_t string_to_function_uint(string data) {
	uint64_t vect = 0;
	for (uint64_t j = 0; j < data.size(); ++j) {
		vect = vect << 1;
		if (data[j] == '1') vect ^= 1;
	}
	return vect;
}

void reg_convert_to_vect(std::ofstream& out, uint64_t value_of_register, int n) {
	for (int j = (n - 1); j > -1; --j)
		out << ((value_of_register >> j) & 1);
}

void enumeration_for_function(vector<uint64_t>& incoming, int n, uint64_t control_value_of_the_length_cycles) {
	int lenght_of_reg = 0;
	uint64_t N = (static_cast<uint64_t>(1) << n);
	ofstream out("C:\\Users\\0\\source\\repos\\generator\\out_file.txt", std::ios::app);
	string teta;
	for (auto value_of_function : incoming) {
		for (uint64_t value_of_register = 0; value_of_register < N; ++value_of_register) {
			generator gen(n, value_of_register, value_of_function);
			gen.f1.transformation();
			if (!value_of_register) {
				out << "function in anf: ";
				gen.f1.print_anf(out);
				out << " and true table: ";
				gen.f2.print_values(out);
			}
			gen.shift_true_table(static_cast<uint64_t>(5) * N);
			gen.cycle_search();
			if (lenght_of_reg < gen.length_of_cycle) lenght_of_reg = gen.length_of_cycle;
			if (control_value_of_the_length_cycles != lenght_of_reg) {
				out << "register: ";
				reg_convert_to_vect(out, value_of_register, n);
				out << " lenght of cycle " << lenght_of_reg << ' ';
				gen.print_gamma(teta);
				out << teta << '\n';
				teta = "";
			}
		}
		lenght_of_reg = 0;
	}
	out.close();
}

void vector_to_value_print_in_file(vector<string>& vector) {
	ofstream out("C:\\Users\\0\\source\\repos\\generator\\out_file.txt", std::ios::app);
	uint64_t out_value;
	for (int index = 0; index < vector.size(); ++index) {
		string to_print;
		for (int h = 0; h < 4; ++h)
			to_print += vector[index][h];
		out_value = string_to_function_uint(to_print);
		out << (int)out_value << ' ';
		to_print.clear();
		for (int k = 1; k < vector[index].size() / 4; ++k) {
			for (int h = 0; h < 4; ++h)
				to_print += vector[index][h+k*4];
			out_value = string_to_function_uint(to_print);
			out << (int)out_value << ' ';
			to_print.clear();
		}
		out << '\n';
	}
	out.close();
}

int main(){
	int n = 8; uint64_t N = (static_cast<uint64_t>(1) << n);
	//uint64_t control_value_of_the_length_cycles = 16;//проверить значение цикла длины control_value_of_the_length_cycles
	//ifstream ifs("D:\\functions\\5_values\\32.txt"); 
	//ofstream out("C:\\Users\\0\\source\\repos\\generator\\out_file.txt", std::ios::app);
	//string str = "function in anf :";
	//vector<string> vector1, vector2, vector_to_out;//вектора для поиска анф, истинности и анф таблицы
	//vector1 = find(ifs, N+1); vector2 = vector1;
	//vector1 = find(ifs, 1); vector2 = vector1;
	//find_anf(vector1); find_tt(vector2, n);
	//for (auto& d : vector2) cout << d << endl;
	//for (auto& d : vector1) cout << d << endl;
	//vector<bool> vect, vect_to_nonlinear;
	/*
	for (int j = 0; j < vector2.size(); ++j) {
		vect.clear();
		string_to_function_vect_bool(vector2[j], vect);
		out << "ANF" << vector1[j] << " TT " << vector2[j];
		transformation(vect, n);
		function_vect_bool_to_string(vector2[j], vect);
		out <<" ANF TABLE "<< vector2[j];
		vector_to_out.push_back(vector2[j]);
		out << endl;
	}
	*/

	//vector_to_value_print_in_file(vector_to_out);
	/*
	set <int> imm_index, nonlinear_index_12, nonlinear_index_20;
	
	for (int j = 0; j < vector2.size(); ++j) {
		string_to_function_vect_bool(vector2[j], vect_to_nonlinear);
		//out << coeff_nonlinear(vect_to_nonlinear, n) << endl;
		int t = coeff_nonlinear(vect_to_nonlinear, n);
		if (t == 20) nonlinear_index_20.insert(j);
		else nonlinear_index_12.insert(j);
		vect_to_nonlinear.clear();
	}
	ifstream ifs2("D:\\functions\\5_values\\32_imm.txt");
	vector<int> index;
	vector<string> ooout;
	string teta;
	while (!ifs2.eof()) {
		getline(ifs2, teta);
		ooout.push_back(teta);
	}
	index.resize(ooout.size());
	for (int k = 0; k < ooout.size(); ++k) {
		istringstream iss(ooout[k], istringstream::in);
		iss >> index[k];
		imm_index.insert(index[k]);
	}
	ifs2.close();

	out << "IMM BEZ 12" << endl;
	for (const auto& element : imm_index)
		if (nonlinear_index_12.find(element) == nonlinear_index_12.end())
			out << element << endl;
	out << "12 BEZ IMM" << endl;
	for (const auto& element : nonlinear_index_12)
		if (imm_index.find(element) == imm_index.end())
			out << element << endl;
	out << "IMM BEZ 20" << endl;
	for (const auto& element : imm_index)
		if (nonlinear_index_20.find(element) == nonlinear_index_20.end())
			out << element << endl;
	out << "20 BEZ IMM" << endl;
	for (const auto& element : nonlinear_index_20)
		if (imm_index.find(element) == imm_index.end())
			out << element << endl;
	*/

	//for (auto d : vect) cout << (int)d;
	//uint64_t incoming=0;
	//vector<uint64_t> inc;
	//for (int j = 0; j < vector2.size(); ++j) {
	//	string_to_function_uint(vector2[j], incoming);
	//	inc.push_back(incoming);
	//}
	//cout << (int)incoming;
	//enumeration_for_function(inc, n, control_value_of_the_length_cycles);
	//enumeration(n);
	/*
	ifstream ifs2("D:\\functions\\5_values\\32_imm.txt");
	vector<int> index;
	vector<string> ooout;
	string teta;
	while (!ifs2.eof()) {
		getline(ifs2, teta);
		ooout.push_back(teta);
	}
	index.resize(ooout.size());

	for (int k = 0; k < ooout.size(); ++k) {
		istringstream iss(ooout[k], istringstream::in);
		iss >> index[k];
	}
	vector<string> anf_5, tt_5; anf_5.resize(ooout.size()); tt_5.resize(ooout.size());
	int indeX=0;

	for (int k = 0; k < 2039; ++k)
		if (k == index[indeX]) {
			anf_5[indeX] = vector1[k];
			tt_5[indeX] = vector2[k];
			++indeX;
		}
	for (int k = 0; k < anf_5.size(); ++k)
		out << "anf is " << anf_5[k] << " and TTable is " << tt_5[k] << endl;
	ifs2.close();
	*/
	//enumeration(n);
	//ifs.close();
	//out.close();
	vector<string> func; func.resize(9);
	func[0] = "1111111111110000110111100101110100000000000011110010000110100010";
	func[1] = "1101100010111100000011111111111100100111010000111111000000000000";
	func[2] = "1111111111110001111101101111001100000000000011100000100100001100";
	func[3] = "1100110101111001110111111111111100110010100001100010000000000000";
	func[4] = "1111111111111101010101001001000100000000000000101010101101101110";
	func[5] = "1000100000001111110111111111111101110111111100000010000000000000";
	func[6] = "1111111111111001000111110011100100000000000001101110000011000110";
	func[7] = "1011111000100001000111111111111101000001110111101110000000000000";
	func[8] = "1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111110001001100000011000011100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000011101100111111001111000";
	vector<vector<bool>> vect; vect.resize(8);
	vector<bool> vv;
	for (int k = 0; k < 8; ++k)
		string_to_function_vect_bool(func[k], vect[k]);
	string_to_function_vect_bool(func[8], vv);
	multi_generator seq(6, 8, vect, vv);
	seq.generate_sequence(0);
	return 0;
}