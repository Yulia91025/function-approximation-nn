#include <iostream>
#include <vector>
#include <time.h>
#include "neurons.cpp"
#include "timer.cpp"
#include <fstream>
#include <GL/freeglut.h>
using namespace std;

// static здесь - это костыль для раздельной компиляции 

static string output_file = "results.txt"; // имя файла, в который выводим результаты 

static int num_in_neurons = 2; // количество входных нейронов
static int num_hid_neurons = 2; // количество скрытых нейронов
static int num_out_neurons = 1; // количество выходных нейронов
static int num_of_experiments = 1; // количество проведённых экспериментов 
static double accuracy = 0.001; // точность ответа

struct DATA
{
	double** input_values; // массив указателей на массивы с входными данными
	double** answers; // массив с указателями на массивы с известными результатами переменных, которые в дальнейшем надо будет угадать
	int max_d = 1; // указатель на массив с максимальной разрядностью эксперимента
	DATA(ifstream& file)
	{
		cout << "Enter accuracy" << endl;
		cin >> accuracy;
		file >> num_in_neurons; // считываем количество "известных" значений 
		++num_in_neurons; //увеличиваем на 1 количество входных нейронов, чтобы был один нейрон смещения
		file >> num_out_neurons; // считываем количество переменных, которые надо будет считать 
		file >> num_of_experiments; // считываем количество экспериментов, которые были проведены (т.е. в которых все значения известны)

		input_values = new double* [num_of_experiments];
		answers = new double* [num_of_experiments];
		for (int n = 0; n < num_of_experiments; ++n)
		{
			input_values[n] = new double[num_in_neurons];
			answers[n] = new double[num_out_neurons];
			for (int i = 0; i < num_in_neurons - 1; ++i) // считываем входные значения
			{
				file >> input_values[n][i];
				int d = dig(input_values[n][i]);
				if (max_d < d) max_d = d; // если разрядность i-ого входного значения больше максимальной разрядности, меняем макс разрядность 
			}
			input_values[n][num_in_neurons - 1] = 1; // задаём значение для нейрона смещения
			for (int i = 0; i < num_out_neurons; ++i) // считываем ответы, которые необходимо получить 
			{
				file >> answers[n][i];
				int d = dig(answers[n][i]);
				if (max_d < d) max_d = d; // если разрядность i-ого ответа больше максимальной разрядности, меняем макс разрядность 
			}
		}
		corr_val();
		num_hid_neurons = 2 * num_in_neurons; // задаём количетво нейронов скрытого слоя в 2 раза больше нейронов входного слоя
		//cout << "MAX D " << max_d << endl;
	}
	int dig(double num) // возвращает 10 в степени разряда числа num
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
	void corr_val() // корректируем значения 
	{
		for (int n = 0; n < num_of_experiments; ++n)
		{
			for (int i = 0; i < num_in_neurons; ++i)
			{                                         //
				input_values[n][i] /= max_d;          // 
			}                                         //        делаем все значения не больше 1
			for (int i = 0; i < num_out_neurons; ++i) //  (т.к. функция активации принимает значения от 0 до 1)
			{                                         // 
				answers[n][i] /= max_d;               //
			}
		}
	}
	double res(int n, double output_value) // выводит в исходной разрядности результат n-ого эксперимента
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
	int curr_experiment = 0; // переменная, показывающая с каким экспериментом работаем 
	double E = 0.7; // скорость обучения
	double A = 0.3; // момент
public:
	Neural_network(DATA data)
	{
		ofstream fout(output_file, ofstream::app); // вывод в файл
		fout << "Accuracy is " << accuracy << endl; // записываем точность в файл 
		Timer t; //заводим счётчик времени
		int num_of_synapses = num_in_neurons * num_hid_neurons + num_hid_neurons * num_out_neurons; // количество синапсов
		double* val_synapses = new double[num_of_synapses]; // массив синапсов
		cout << "Do you already have synapses values?" << endl;
		string ans_1 = "";
		cin >> ans_1;
		if (ans_1 == "Yes" || ans_1 == "YES" || ans_1 == "yes") import_syn(val_synapses);
		else for (int i = 0; i < num_of_synapses; ++i) // задаём начальные значения синапсов рандомно
		{
			val_synapses[i] = (double)(rand() % 1000 + 1) * (rand() % 10 - 5) / 10000;
		}
		cout << "E = " << E << ";  A = " << A << endl;
		fout << "E = " << E << ";  A = " << A << endl;

		for (int i = 0; i < num_out_neurons; ++i) // присваиваем идеальным значениям значения первого набора ответов
		{
			Ideal[i] = data.answers[0][i];
		}
		for (int i = 0; i < num_in_neurons; ++i)
		{
			Input_Neurons[i] = Input_Neuron(data.input_values[0][i]);  // инициализируем i-й нейрон входного слоя первым набором значений
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
		}
		for (int n = 0; n < num_of_experiments; ++n)
		{
			Error[n] = 0; // зануляем значения ошибок
			curr_experiment = n;
			error_rate(); // вычисляем ошибку для текущего эксперимента 
			change_input_values(data.input_values[(n + 1) % num_of_experiments], data.answers[((n + 1) % num_of_experiments)]); // меняем экмперимент 
		}
		cout << "-------------------------------------------------------------------------------------" << endl;
		backpropagation(data); // запускаем метод обратного распространения ошибки
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
	void hid_neu_init(int i) // функция инициализации i-ого нейрона скрытого слоя
	{
		double in_val = 0; // переменная значения входных данных для i-ого нейрона скрытого слоя
		for (int j = 0; j < num_in_neurons; ++j)
		{
			in_val += Input_Neurons[j].output * Synapses[i + j * num_in_neurons].value;
		}
		Hidden_Neurons[i] = Neuron(in_val);
		//cout << "Hidden neuon " << i + 1 << "   " << in_val << endl;
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
		cout << "Current number " << curr_experiment << endl;
		Error[curr_experiment] = 0; // зануляем ошибку для новой итерации
		for (int i = 0; i < num_out_neurons; ++i)
		{
			Error[curr_experiment] += (Ideal[i] - Output_Neurons[i].output) * (Ideal[i] - Output_Neurons[i].output); // метод среднего квадрата ошибки (складываем квадраты разностей идеального значения и полученного, затем делим на количество выходных нейронов)
		}
		Error[curr_experiment] /= num_out_neurons; // делим полученную сумму квадратов на количество выходных нейронов
		cout << "Error rate: " << Error[curr_experiment] << endl;
	}
	double GRAD(Synapse& w) // вычисление градиента для синапса w
	{
		return ((*w.Out).output * (*w.In).delta);
	}
	void backpropagation(DATA data) // Метод обратного распространения ошибки
	{
		int changing_synapses_counter = 0; // показывает сколько раз мы меняли все синапсы 
		double* syn_change_val = new double[num_of_syn]; // вектор величин, на которые необходимо будет изменить синапсы
		for (int i = 0; i < num_of_syn; ++i)
		{
			syn_change_val[i] = 0;
		}
		double sum_err = num_of_experiments; // изначально сумма всех ошибок максимальная (1 * num_of_experiments)
		while (sum_err > accuracy) // Метод ОР ошибки работает пока ошибка больше заданной точности 
		{
			for (int i = 0; i < num_out_neurons; ++i)
			{
				//Output_Neurons[i].delta = (Ideal[i] - Output_Neurons[i].output) * (1 - Output_Neurons[i].output) * Output_Neurons[i].output; // находим дельты для выходных нейронов для сигмоида
				Output_Neurons[i].delta = (Ideal[i] - Output_Neurons[i].output) * (1 - Output_Neurons[i].output * Output_Neurons[i].output); // находим дельты для выходных нейронов для гиперболического тангенса
				for (int j = 0; j < num_hid_neurons; ++j) // находим величины, на которые необходимо изменить синапсы, связанные с выходными нейронами
				{
					if (changing_synapses_counter > 0)
						syn_change_val[i + j + num_in_neurons * num_hid_neurons] = E * GRAD(Synapses[i + j + num_in_neurons * num_hid_neurons]) + A * syn_change_val[i + j + num_in_neurons * num_hid_neurons];
					else syn_change_val[i + j + num_in_neurons * num_hid_neurons] = E * GRAD(Synapses[i + j + num_in_neurons * num_hid_neurons]);
				}
			}
			for (int i = 0; i < num_hid_neurons; ++i)
			{
				for (int j = 0; j < num_out_neurons; ++j) // находим дельты для скрытых нейронов
				{
					//Hidden_Neurons[i].delta += (1 - Hidden_Neurons[i].output) * Hidden_Neurons[i].output * Synapses[i + j + num_in_neurons * num_hid_neurons].value * Output_Neurons[j].delta; // для сигмоида 
					Hidden_Neurons[i].delta += (1 - Hidden_Neurons[i].output * Hidden_Neurons[i].output) * Synapses[i + j + num_in_neurons * num_hid_neurons].value * Output_Neurons[j].delta; // для гиперболического тангенса 
				}
				for (int j = 0; j < num_in_neurons; ++j) // находим величины, на которые необходимо изменить синапсы, входящие в скрытые нейроны
				{
					if (changing_synapses_counter > 0)
						syn_change_val[i + j * num_in_neurons] = E * GRAD(Synapses[i + j * num_in_neurons]) + A * syn_change_val[i + j * num_in_neurons];
					else syn_change_val[i + j * num_in_neurons] = E * GRAD(Synapses[i + j * num_in_neurons]);
				}
			}
			for (int i = 0; i < num_of_syn; ++i) // меняем значения синапсов
			{
				Synapses[i].value += syn_change_val[i];
				//cout << "Synapse " << i + 1 << "value: " << Synapses[i].value << endl;
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
			if (num_of_experiments > 1)
			{
				change_input_values(data.input_values[changing_synapses_counter % num_of_experiments], data.answers[changing_synapses_counter % num_of_experiments]);
				curr_experiment = changing_synapses_counter % num_of_experiments; // меняем текущий номер эксперимента 
			}
			sum_err = 0;
			for (int n = 0; n < num_of_experiments; ++n)
			{
				sum_err += Error[n];
			}
			cout << "Sum of errors " << sum_err << endl;
		}
	}
	void change_input_values(double* input_values, double* answers) // меняем входные значения
	{
		for (int i = 0; i < num_out_neurons; ++i) // присваиваем идеальным значениям значения ответов
		{
			Ideal[i] = answers[i];
		}
		for (int i = 0; i < num_in_neurons; ++i)
		{
			Input_Neurons[i] = Input_Neuron(input_values[i]);  // инициализируем i-й нейрон входного слоя
		}
		for (int i = 0; i < num_hid_neurons; ++i)
		{
			hid_neu_init(i); // инициализируем i-й нейрон скрытого слоя
		}
		for (int i = 0; i < num_out_neurons; ++i)
		{
			out_neu_init(i); // инициализируем i-й нейрон выходного слоя
			cout << "Result " << i + 1 << " output neuron: " << Output_Neurons[i].output << endl;

		}
		cout << "-------------------------------------------------------------------------------------" << endl;
	}
	void guess(int max_d) // метод, угадывающий ответ 
	{
		double* input_values = new double[num_in_neurons]; // массив входных значений 
		for (int i = 0; i < num_in_neurons - 1; ++i)
		{
			cout << "Enter " << i + 1 << " value: ";
			cin >> input_values[i];
			Input_Neurons[i] = Input_Neuron(input_values[i] / max_d);  // инициализируем i-й нейрон входного слоя
		}
		for (int i = 0; i < num_hid_neurons; ++i)
		{
			hid_neu_init(i); // инициализируем i-й нейрон скрытого слоя
		}
		for (int i = 0; i < num_out_neurons; ++i)
		{
			out_neu_init(i); // инициализируем i-й нейрон выходного слоя
			cout << "Answer " << i + 1 << " is " << Output_Neurons[i].output * max_d << endl;

		}
		cout << "-------------------------------------------------------------------------------------" << endl;
	}
	void save_syn() // метод, сохраняющий значения синапсов 
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
	void import_syn(double* val_synapses)// считывает значения синапсов из файла
	{
		string filename;
		cout << "Enter filename with synapses values:" << endl;
		cin >> filename;
		filename += ".txt";
		ifstream file(filename);
		// надо сделать защиту от меньшего кольчества синапсов в файле
		for (int i = 0; i < num_of_syn; ++i)
		{
			file >> val_synapses[i];
		}
		file.close();
	}
	double** for_graph(DATA data) // метод для создания координат для графика 
	{
		double** coord = new double* [num_of_experiments]; // массив координат 
		for (int n = 0; n < num_of_experiments; ++n) // пробегаемся по всем экспериментам 
		{
			coord[n] = new double[num_in_neurons + num_out_neurons - 1]; // вычитаем единицу, потому что до этого увеличили количество входных нейронов для одного нейрона смещения 
			for (int i = 0; i < num_in_neurons - 1; ++i)
			{
				coord[n][i] = data.input_values[n][i] * data.max_d;
				Input_Neurons[i] = Input_Neuron(data.input_values[n][i]);  // инициализируем i-й нейрон входного слоя
			}
			Input_Neurons[num_in_neurons - 1] = Input_Neuron(data.input_values[n][num_in_neurons - 1]);  // инициализируем скрытый нейрон входного слоя
			for (int i = 0; i < num_hid_neurons; ++i)
			{
				hid_neu_init(i); // инициализируем i-й нейрон скрытого слоя
			}
			for (int i = 0; i < num_out_neurons; ++i)
			{
				out_neu_init(i); // инициализируем i-й нейрон выходного слоя
				coord[n][num_in_neurons + i - 1] = Output_Neurons[i].output * data.max_d;
			}
			cout << "Look at this results " << coord[n][0] << "               " << coord[n][num_in_neurons - 1] << endl;
		}
		return coord; // возвращаем массив массивов с координатами 
	}
};

