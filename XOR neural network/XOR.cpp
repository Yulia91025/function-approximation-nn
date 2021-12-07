#include <iostream>
#include <math.h>
#include <vector>
using namespace std;
const int num_in_neurons = 2; // количество входных нейронов
const int num_hid_neurons = 2; // количество скрытых нейронов
const int num_out_neurons = 1; // количество выходных нейронов
const double E = 0.7; // скорость обучения
const double A = 0.3; // момент
const double accuracy = 0.002;
double activation_function(double in) // функция активации
{
	return (1 / (1 + exp(-in))); // сигмоид
}
struct Neuron // структура для скрытых нейронов и нейронов выхода
{
	double input;
	double output;
	double delta = 0; // переменная для ввычисления новых значений синапсов, связанных с этим нейроном 
	Neuron()
	{
		input = 0;
		output = 0;
	}
	Neuron(double in)
	{
		input = in;
		output = activation_function(in);
	}
};
struct Input_Neuron : public Neuron  // структура для входных нейронов
{
	Input_Neuron()
	{
		output = 0;
	}
	Input_Neuron(double in)
	{
		output = in;
	}
};
struct Synapse // структура для синапсов
{
	Neuron* Out; // указатель на нейрон, из которого синапс выходит 
	Neuron* In; // указатель на нейрон, в который синапс входит 
	double value; // значение синапса
	Synapse()
	{
		In = new Neuron();
		Out = new Neuron();
		value = 0;
	}
	Synapse(Neuron &Out_n, Neuron &In_n, double val)
	{
		Out = &Out_n;
		In = &In_n;
		value = val;
	}
};
class Neural_network
{
	Input_Neuron Input_Neurons[num_in_neurons];
	Neuron Hidden_Neurons[num_hid_neurons];
	Neuron Output_Neurons[num_out_neurons];
	const int num_of_syn = num_in_neurons * num_hid_neurons + num_hid_neurons * num_out_neurons;
	Synapse* Synapses = new Synapse[num_of_syn];
	int counter_syn = 0;
	double Ideal[num_out_neurons];
	double Error;
public:
	Neural_network(double* input_values, double* val_synapses, double* answer)
	{
		for (int i = 0; i < num_out_neurons; ++i) // присваиваем идеальным значениям значения ответов
		{
			Ideal[i] = answer[i]; 
		}
		for (int i = 0; i < num_in_neurons; ++i)
		{
			Input_Neurons[i] = Input_Neuron(input_values[i]);  // инициализируем i-й нейрон входного слоя
			for (int j = 0; j < num_hid_neurons; ++j)
			{
				Synapses[counter_syn] = Synapse(Input_Neurons[i], Hidden_Neurons[j], val_synapses[counter_syn]); // инициализируем синапсы для каждого входного нейрона
				++counter_syn;
			}
		}
		for (int i = 0; i < num_hid_neurons; ++i)
		{
			hid_neu_init(i); // инициализируем i-й нейрон скрытого слоя
			for (int j = 0; j < num_out_neurons; ++j)
			{
				Synapses[counter_syn] = Synapse(Hidden_Neurons[i], Output_Neurons[j], val_synapses[counter_syn]); // инициализируем исходящие синапсы для каждого скрытого нейрона
				++counter_syn;
			}
		}
		for (int i = 0; i < num_out_neurons; ++i)
		{
			out_neu_init(i); // инициализируем i-й нейрон выходного слоя
			cout << "Result " << i + 1 << " output neuron: " << Output_Neurons[i].output << endl;
		}
		error_rate(); // вычисляем ошибку
		cout << "-------------------------------------------------------------------------------------" << endl;
		backpropagation();
	}
	void hid_neu_init(int i) // функция инициализации i-ого нейрона скрытого слоя
	{
		double in_val = 0; // переменная значения входных данных для i-ого нейрона скрытого слоя
		for (int j = 0; j < num_in_neurons; ++j)
		{
			in_val += Input_Neurons[j].output * Synapses[i + j * num_in_neurons].value;
		}
		Hidden_Neurons[i] = Neuron(in_val);
	}
	void out_neu_init(int i) // функция инициализации i-ого нейрона выходного слоя
	{
		double in_val = 0; // переменная значения входных данных для i-ого нейрона выходного слоя
		for (int j = 0; j < num_hid_neurons; ++j)
		{
			in_val += Hidden_Neurons[j].output * Synapses[i + j + num_in_neurons * num_hid_neurons].value;
		}
		Output_Neurons[i] = Neuron(in_val);
	}
	void error_rate() // функция вычисления ошибки
	{
		for (int i = 0; i < num_out_neurons; ++i)
		{
			Error += (Ideal[i] - Output_Neurons[i].output) * (Ideal[i] - Output_Neurons[i].output); // метод среднего квадрата ошибки (складываем квадраты разностей идеального значения и полученного, затем делим на количество выходных нейронов)
		}
		Error /= num_out_neurons; // делим полученную сумму квадратов на количество выходных нейронов
		cout << "Error rate: " << Error << endl;
	}
	double GRAD(Synapse &w) // вычисление градиента для синапса w
	{
		return ((*w.Out).output * (*w.In).delta);
	}
	void backpropagation() // Метод обратного распространения ошибки
	{
		int changing_synapses_counter = 0; // показывает сколько раз мы меняли все синапсы 
		vector<double> syn_change_val(100000); // вектор величин, на которые необходимо будет изменить синапсы
		while (Error > accuracy)
		{
			Error = 0; // зануляем ошибку для новой итерации
			for (int i = changing_synapses_counter * num_of_syn; i < (changing_synapses_counter + 1)* num_of_syn ; ++i)
			{
				syn_change_val[i] = 0;
			}
			for (int i = 0; i < num_out_neurons; ++i) 
			{
				Output_Neurons[i].delta = (Ideal[i] - Output_Neurons[i].output) * (1 - Output_Neurons[i].output) * Output_Neurons[i].output; // находим дельты для выходных нейронов
				for (int j = 0; j < num_hid_neurons; ++j) // находим величины, на которые необходимо изменить синапсы, связанные с выходными нейронами
				{
					if (changing_synapses_counter > 0)
						syn_change_val[i + j + num_in_neurons * num_hid_neurons + changing_synapses_counter * num_of_syn] = E * GRAD(Synapses[i + j + num_in_neurons * num_hid_neurons]) + A * syn_change_val[i + j + num_in_neurons * num_hid_neurons + (changing_synapses_counter - 1) * num_of_syn];
					else syn_change_val[i + j + num_in_neurons * num_hid_neurons + changing_synapses_counter * num_of_syn] = E * GRAD(Synapses[i + j + num_in_neurons * num_hid_neurons]);
				}
			}
			for (int i = 0; i < num_hid_neurons; ++i) 
			{
				for (int j = 0; j < num_out_neurons; ++j) // находим дельты для скрытых нейронов
				{
					Hidden_Neurons[i].delta += (1 - Hidden_Neurons[i].output) * Hidden_Neurons[i].output * Synapses[i + j + num_in_neurons * num_hid_neurons].value * Output_Neurons[j].delta;
				}
				for (int j = 0; j < num_in_neurons; ++j) // находим величины, на которые необходимо изменить синапсы, входящие в скрытые нейроны
				{
					if (changing_synapses_counter > 0)
						syn_change_val[i + j * num_in_neurons + changing_synapses_counter * num_of_syn] = E * GRAD(Synapses[i + j * num_in_neurons]) + A * syn_change_val[i + j * num_in_neurons + (changing_synapses_counter - 1) * num_of_syn];
					else syn_change_val[i + j * num_in_neurons + changing_synapses_counter * num_of_syn] = E * GRAD(Synapses[i + j * num_in_neurons]);
				}
			}
			for (int i = 0; i < num_of_syn; ++i) // меняем значения синапсов
			{
				Synapses[i].value += syn_change_val[i + changing_synapses_counter * num_of_syn];
			}
			for (int i = 0; i < num_hid_neurons; ++i)
			{
				hid_neu_init(i); // пересчитываем значения i-ого нейрона скрытого слоя
			}
			for (int i = 0; i < num_out_neurons; ++i)
			{
				out_neu_init(i); // пересчитываем значения i-ого нейрона выходного слоя
				cout << "Result of " << i + 1 << " output neuron: " << Output_Neurons[i].output << endl;
			}
			error_rate(); // вычисляем ошибку
			++changing_synapses_counter;
			cout << changing_synapses_counter << " times synapse values changed" << endl;
			cout << "-------------------------------------------------------------------------------------" << endl;
		}
	}
};
int main()
{
	double* input_values = new double[num_in_neurons]; // массив входных значений
	for (int i = 0; i < num_in_neurons; ++i) // считываем входные значения
	{
		cout << "Enter " << i + 1 << " value:" << endl;
		cin >> input_values[i];
	}
	int num_of_synapses = num_in_neurons * num_hid_neurons + num_hid_neurons * num_out_neurons; // количество синапсов
	double* val_synapses = new double[num_of_synapses]; // массив синапсов
	for (int i = 0; i < num_of_synapses; ++i) // считываем синапсы
	{
		cout << "Enter value of " << i + 1 << " synapse:" << endl;
		cin >> val_synapses[i];
	}
	double answer[num_out_neurons];
	for (int i = 0; i < num_out_neurons; ++i) // считываем ответы, которые необходимо получить 
	{
		cout << "Enter answer: " << endl;
		cin >> answer[i];
	}
	Neural_network Net(input_values, val_synapses, answer); // создаём нейросеть 
	return 0;
};