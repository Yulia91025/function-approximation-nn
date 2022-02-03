#include <iostream>
#include <vector>
#include <time.h>
#include "neurons.cpp"
#include "timer.cpp"
#include <fstream>
#include <GL/freeglut.h>
using namespace std;

// static ����� - ��� ������� ��� ���������� ���������� 

static string output_file = "results.txt"; // ��� �����, � ������� ������� ���������� 

static int num_in_neurons = 2; // ���������� ������� ��������
static int num_hid_neurons = 2; // ���������� ������� ��������
static int num_out_neurons = 1; // ���������� �������� ��������
static int num_of_experiments = 1; // ���������� ���������� ������������� 
static double accuracy = 0.001; // �������� ������

struct DATA
{
	double** input_values; // ������ ���������� �� ������� � �������� �������
	double** answers; // ������ � ����������� �� ������� � ���������� ������������ ����������, ������� � ���������� ���� ����� �������
	int max_d = 1; // ��������� �� ������ � ������������ ������������ ������������
	DATA(ifstream& file)
	{
		cout << "Enter accuracy" << endl;
		cin >> accuracy;
		file >> num_in_neurons; // ��������� ���������� "���������" �������� 
		++num_in_neurons; //����������� �� 1 ���������� ������� ��������, ����� ��� ���� ������ ��������
		file >> num_out_neurons; // ��������� ���������� ����������, ������� ���� ����� ������� 
		file >> num_of_experiments; // ��������� ���������� �������������, ������� ���� ��������� (�.�. � ������� ��� �������� ��������)

		input_values = new double* [num_of_experiments];
		answers = new double* [num_of_experiments];
		for (int n = 0; n < num_of_experiments; ++n)
		{
			input_values[n] = new double[num_in_neurons];
			answers[n] = new double[num_out_neurons];
			for (int i = 0; i < num_in_neurons - 1; ++i) // ��������� ������� ��������
			{
				file >> input_values[n][i];
				int d = dig(input_values[n][i]);
				if (max_d < d) max_d = d; // ���� ����������� i-��� �������� �������� ������ ������������ �����������, ������ ���� ����������� 
			}
			input_values[n][num_in_neurons - 1] = 1; // ����� �������� ��� ������� ��������
			for (int i = 0; i < num_out_neurons; ++i) // ��������� ������, ������� ���������� �������� 
			{
				file >> answers[n][i];
				int d = dig(answers[n][i]);
				if (max_d < d) max_d = d; // ���� ����������� i-��� ������ ������ ������������ �����������, ������ ���� ����������� 
			}
		}
		corr_val();
		num_hid_neurons = 2 * num_in_neurons; // ����� ��������� �������� �������� ���� � 2 ���� ������ �������� �������� ����
		//cout << "MAX D " << max_d << endl;
	}
	int dig(double num) // ���������� 10 � ������� ������� ����� num
	{
		double n = abs(num);
		int counter = 0;
		while (n > 1)
		{
			n /= 10;
			++counter;
		}
		return pow(10, counter);
	}
	void corr_val() // ������������ �������� 
	{
		for (int n = 0; n < num_of_experiments; ++n)
		{
			for (int i = 0; i < num_in_neurons; ++i)
			{                                         //
				input_values[n][i] /= max_d;          // 
			}                                         //        ������ ��� �������� �� ������ 1
			for (int i = 0; i < num_out_neurons; ++i) //  (�.�. ������� ��������� ��������� �������� �� 0 �� 1)
			{                                         // 
				answers[n][i] /= max_d;               //
			}
		}
	}
	double res(int n, double output_value) // ������� � �������� ����������� ��������� n-��� ������������
	{
		return output_value * max_d;
	}
};

class Neural_network
{
	Input_Neuron* Input_Neurons = new Input_Neuron[num_in_neurons];
	Neuron* Hidden_Neurons = new Neuron[num_hid_neurons];
	Neuron* Output_Neurons = new Neuron[num_out_neurons];
	const int num_of_syn = num_in_neurons * num_hid_neurons + num_hid_neurons * num_out_neurons;
	Synapse* Synapses = new Synapse[num_of_syn];
	int counter_syn = 0;
	double* Ideal = new double[num_out_neurons];
	double* Error = new double[num_of_experiments];
	int curr_experiment = 0; // ����������, ������������ � ����� ������������� �������� 
	double E = 0.7; // �������� ��������
	double A = 0.3; // ������
public:
	Neural_network(DATA data)
	{
		ofstream fout(output_file, ofstream::app); // ����� � ����
		fout << "Accuracy is " << accuracy << endl; // ���������� �������� � ���� 
		Timer t; //������� ������� �������
		int num_of_synapses = num_in_neurons * num_hid_neurons + num_hid_neurons * num_out_neurons; // ���������� ��������
		double* val_synapses = new double[num_of_synapses]; // ������ ��������
		cout << "Do you already have synapses values?" << endl;
		string ans_1 = "";
		cin >> ans_1;
		if (ans_1 == "Yes" || ans_1 == "YES" || ans_1 == "yes") import_syn(val_synapses);
		else for (int i = 0; i < num_of_synapses; ++i) // ����� ��������� �������� �������� ��������
		{
			val_synapses[i] = (double)(rand() % 1000 + 1) * (rand() % 10 - 5) / 10000;
		}
		cout << "E = " << E << ";  A = " << A << endl;
		fout << "E = " << E << ";  A = " << A << endl;

		for (int i = 0; i < num_out_neurons; ++i) // ����������� ��������� ��������� �������� ������� ������ �������
		{
			Ideal[i] = data.answers[0][i];
		}
		for (int i = 0; i < num_in_neurons; ++i)
		{
			Input_Neurons[i] = Input_Neuron(data.input_values[0][i]);  // �������������� i-� ������ �������� ���� ������ ������� ��������
			for (int j = 0; j < num_hid_neurons; ++j)
			{
				Synapses[counter_syn] = Synapse(Input_Neurons[i], Hidden_Neurons[j], val_synapses[counter_syn]); // �������������� ������� ��� ������� �������� �������
				++counter_syn;
			}
		}
		for (int i = 0; i < num_hid_neurons; ++i)
		{
			hid_neu_init(i); // �������������� i-� ������ �������� ����
			for (int j = 0; j < num_out_neurons; ++j)
			{
				Synapses[counter_syn] = Synapse(Hidden_Neurons[i], Output_Neurons[j], val_synapses[counter_syn]); // �������������� ��������� ������� ��� ������� �������� �������
				++counter_syn;
			}
		}
		for (int i = 0; i < num_out_neurons; ++i)
		{
			out_neu_init(i); // �������������� i-� ������ ��������� ����
		}
		for (int n = 0; n < num_of_experiments; ++n)
		{
			Error[n] = 0; // �������� �������� ������
			curr_experiment = n;
			error_rate(); // ��������� ������ ��� �������� ������������ 
			change_input_values(data.input_values[(n + 1) % num_of_experiments], data.answers[((n + 1) % num_of_experiments)]); // ������ ����������� 
		}
		cout << "-------------------------------------------------------------------------------------" << endl;
		backpropagation(data); // ��������� ����� ��������� ��������������� ������
		cout << "Neural network worked for " << t.elapsed() << " seconds" << endl;
		fout << "Neural network worked for " << t.elapsed() << " seconds" << endl;
		fout << "-------------------------------------------------------------------------------------" << endl;
		for_graph(data);
		bool flag = 1;
		while (flag)
		{
			guess(data.max_d);
			cout << "Do you want to finish?" << endl;
			string ans = "";
			cin >> ans;
			if (ans == "Yes" || ans == "YES" || ans == "yes") flag = 0;
		}
		cout << "Do you want to save synapses values?" << endl;
		string ans_2 = "";
		cin >> ans_2;
		if (ans_2 == "Yes" || ans_2 == "YES" || ans_2 == "yes") save_syn();
	}
	void hid_neu_init(int i) // ������� ������������� i-��� ������� �������� ����
	{
		double in_val = 0; // ���������� �������� ������� ������ ��� i-��� ������� �������� ����
		for (int j = 0; j < num_in_neurons; ++j)
		{
			in_val += Input_Neurons[j].output * Synapses[i + j * num_in_neurons].value;
		}
		Hidden_Neurons[i] = Neuron(in_val);
		//cout << "Hidden neuon " << i + 1 << "   " << in_val << endl;
	}
	void out_neu_init(int i) // ������� ������������� i-��� ������� ��������� ����
	{
		double in_val = 0; // ���������� �������� ������� ������ ��� i-��� ������� ��������� ����
		for (int j = 0; j < num_hid_neurons; ++j)
		{
			in_val += Hidden_Neurons[j].output * Synapses[i + j + num_in_neurons * num_hid_neurons].value;
		}
		Output_Neurons[i] = Neuron(in_val);
	}
	void error_rate() // ������� ���������� ������
	{
		cout << "Current number " << curr_experiment << endl;
		Error[curr_experiment] = 0; // �������� ������ ��� ����� ��������
		for (int i = 0; i < num_out_neurons; ++i)
		{
			Error[curr_experiment] += (Ideal[i] - Output_Neurons[i].output) * (Ideal[i] - Output_Neurons[i].output); // ����� �������� �������� ������ (���������� �������� ��������� ���������� �������� � �����������, ����� ����� �� ���������� �������� ��������)
		}
		Error[curr_experiment] /= num_out_neurons; // ����� ���������� ����� ��������� �� ���������� �������� ��������
		cout << "Error rate: " << Error[curr_experiment] << endl;
	}
	double GRAD(Synapse& w) // ���������� ��������� ��� ������� w
	{
		return ((*w.Out).output * (*w.In).delta);
	}
	void backpropagation(DATA data) // ����� ��������� ��������������� ������
	{
		int changing_synapses_counter = 0; // ���������� ������� ��� �� ������ ��� ������� 
		double* syn_change_val = new double[num_of_syn]; // ������ �������, �� ������� ���������� ����� �������� �������
		for (int i = 0; i < num_of_syn; ++i)
		{
			syn_change_val[i] = 0;
		}
		double sum_err = num_of_experiments; // ���������� ����� ���� ������ ������������ (1 * num_of_experiments)
		while (sum_err > accuracy) // ����� �� ������ �������� ���� ������ ������ �������� �������� 
		{
			for (int i = 0; i < num_out_neurons; ++i)
			{
				//Output_Neurons[i].delta = (Ideal[i] - Output_Neurons[i].output) * (1 - Output_Neurons[i].output) * Output_Neurons[i].output; // ������� ������ ��� �������� �������� ��� ��������
				Output_Neurons[i].delta = (Ideal[i] - Output_Neurons[i].output) * (1 - Output_Neurons[i].output * Output_Neurons[i].output); // ������� ������ ��� �������� �������� ��� ���������������� ��������
				for (int j = 0; j < num_hid_neurons; ++j) // ������� ��������, �� ������� ���������� �������� �������, ��������� � ��������� ���������
				{
					if (changing_synapses_counter > 0)
						syn_change_val[i + j + num_in_neurons * num_hid_neurons] = E * GRAD(Synapses[i + j + num_in_neurons * num_hid_neurons]) + A * syn_change_val[i + j + num_in_neurons * num_hid_neurons];
					else syn_change_val[i + j + num_in_neurons * num_hid_neurons] = E * GRAD(Synapses[i + j + num_in_neurons * num_hid_neurons]);
				}
			}
			for (int i = 0; i < num_hid_neurons; ++i)
			{
				for (int j = 0; j < num_out_neurons; ++j) // ������� ������ ��� ������� ��������
				{
					//Hidden_Neurons[i].delta += (1 - Hidden_Neurons[i].output) * Hidden_Neurons[i].output * Synapses[i + j + num_in_neurons * num_hid_neurons].value * Output_Neurons[j].delta; // ��� �������� 
					Hidden_Neurons[i].delta += (1 - Hidden_Neurons[i].output * Hidden_Neurons[i].output) * Synapses[i + j + num_in_neurons * num_hid_neurons].value * Output_Neurons[j].delta; // ��� ���������������� �������� 
				}
				for (int j = 0; j < num_in_neurons; ++j) // ������� ��������, �� ������� ���������� �������� �������, �������� � ������� �������
				{
					if (changing_synapses_counter > 0)
						syn_change_val[i + j * num_in_neurons] = E * GRAD(Synapses[i + j * num_in_neurons]) + A * syn_change_val[i + j * num_in_neurons];
					else syn_change_val[i + j * num_in_neurons] = E * GRAD(Synapses[i + j * num_in_neurons]);
				}
			}
			for (int i = 0; i < num_of_syn; ++i) // ������ �������� ��������
			{
				Synapses[i].value += syn_change_val[i];
				//cout << "Synapse " << i + 1 << "value: " << Synapses[i].value << endl;
			}
			for (int i = 0; i < num_hid_neurons; ++i)
			{
				hid_neu_init(i); // ������������� �������� i-��� ������� �������� ����
			}
			for (int i = 0; i < num_out_neurons; ++i)
			{
				out_neu_init(i); // ������������� �������� i-��� ������� ��������� ����
				cout << "Result of " << i + 1 << " output neuron: " << Output_Neurons[i].output << endl;
			}
			error_rate(); // ��������� ������
			++changing_synapses_counter;
			cout << changing_synapses_counter << " times synapse values changed" << endl;
			cout << "-------------------------------------------------------------------------------------" << endl;
			if (num_of_experiments > 1)
			{
				change_input_values(data.input_values[changing_synapses_counter % num_of_experiments], data.answers[changing_synapses_counter % num_of_experiments]);
				curr_experiment = changing_synapses_counter % num_of_experiments; // ������ ������� ����� ������������ 
			}
			sum_err = 0;
			for (int n = 0; n < num_of_experiments; ++n)
			{
				sum_err += Error[n];
			}
			cout << "Sum of errors " << sum_err << endl;
		}
	}
	void change_input_values(double* input_values, double* answers) // ������ ������� ��������
	{
		for (int i = 0; i < num_out_neurons; ++i) // ����������� ��������� ��������� �������� �������
		{
			Ideal[i] = answers[i];
		}
		for (int i = 0; i < num_in_neurons; ++i)
		{
			Input_Neurons[i] = Input_Neuron(input_values[i]);  // �������������� i-� ������ �������� ����
		}
		for (int i = 0; i < num_hid_neurons; ++i)
		{
			hid_neu_init(i); // �������������� i-� ������ �������� ����
		}
		for (int i = 0; i < num_out_neurons; ++i)
		{
			out_neu_init(i); // �������������� i-� ������ ��������� ����
			cout << "Result " << i + 1 << " output neuron: " << Output_Neurons[i].output << endl;

		}
		cout << "-------------------------------------------------------------------------------------" << endl;
	}
	void guess(int max_d) // �����, ����������� ����� 
	{
		double* input_values = new double[num_in_neurons]; // ������ ������� �������� 
		for (int i = 0; i < num_in_neurons - 1; ++i)
		{
			cout << "Enter " << i + 1 << " value: ";
			cin >> input_values[i];
			Input_Neurons[i] = Input_Neuron(input_values[i] / max_d);  // �������������� i-� ������ �������� ����
		}
		for (int i = 0; i < num_hid_neurons; ++i)
		{
			hid_neu_init(i); // �������������� i-� ������ �������� ����
		}
		for (int i = 0; i < num_out_neurons; ++i)
		{
			out_neu_init(i); // �������������� i-� ������ ��������� ����
			cout << "Answer " << i + 1 << " is " << Output_Neurons[i].output * max_d << endl;

		}
		cout << "-------------------------------------------------------------------------------------" << endl;
	}
	void save_syn() // �����, ����������� �������� �������� 
	{
		string filename;
		cout << "Enter filename where synapses values will be saved:" << endl;
		cin >> filename;
		filename += ".txt";
		ofstream fout;
		fout.open(filename);
		for (int i = 0; i < num_of_syn; ++i)
		{
			fout << Synapses[i].value << endl;
		}
		fout.close();
	}
	void import_syn(double* val_synapses)// ��������� �������� �������� �� �����
	{
		string filename;
		cout << "Enter filename with synapses values:" << endl;
		cin >> filename;
		filename += ".txt";
		ifstream file(filename);
		// ���� ������� ������ �� �������� ���������� �������� � �����
		for (int i = 0; i < num_of_syn; ++i)
		{
			file >> val_synapses[i];
		}
		file.close();
	}
	double** for_graph(DATA data) // ����� ��� �������� ��������� ��� ������� 
	{
		double** coord = new double* [num_of_experiments]; // ������ ��������� 
		for (int n = 0; n < num_of_experiments; ++n) // ����������� �� ���� ������������� 
		{
			coord[n] = new double[num_in_neurons + num_out_neurons - 1]; // �������� �������, ������ ��� �� ����� ��������� ���������� ������� �������� ��� ������ ������� �������� 
			for (int i = 0; i < num_in_neurons - 1; ++i)
			{
				coord[n][i] = data.input_values[n][i] * data.max_d;
				Input_Neurons[i] = Input_Neuron(data.input_values[n][i]);  // �������������� i-� ������ �������� ����
			}
			Input_Neurons[num_in_neurons - 1] = Input_Neuron(data.input_values[n][num_in_neurons - 1]);  // �������������� ������� ������ �������� ����
			for (int i = 0; i < num_hid_neurons; ++i)
			{
				hid_neu_init(i); // �������������� i-� ������ �������� ����
			}
			for (int i = 0; i < num_out_neurons; ++i)
			{
				out_neu_init(i); // �������������� i-� ������ ��������� ����
				coord[n][num_in_neurons + i - 1] = Output_Neurons[i].output * data.max_d;
			}
			cout << "Look at this results " << coord[n][0] << "               " << coord[n][num_in_neurons - 1] << endl;
		}
		return coord; // ���������� ������ �������� � ������������ 
	}
};

