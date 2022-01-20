#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int num_of_arg = 1;
int num_of_val = 1;
int num_of_experiments = 1;

void function(double* arr, double* answers) 
{
	for (int i = 0; i < num_of_val; ++i)
	{
		answers[i] = 0;
		for (int j = 0; j < num_of_arg; ++j)
		{
			answers[i] = arr[j]* arr[j];
		}
	}
}

struct range
{
	double beg = 0;
	double end = 3;
	range(){}
	range(double b, double e):beg(b),end(e){}
	friend istream& operator>> (istream& in, range& r)
	{
		in >> r.beg;
		in >> r.end;
		return in;
	}
};


int main()
{
	string filename;
	cout << "Enter filename where test will be saved:" << endl;
	cin >> filename;
	filename += ".txt";
	ofstream fout;
	fout.open(filename);

	cout << "Enter number of arguments\n";
	cin >> num_of_arg; // считываем количество аргументов
	cout << "Enter number of values\n";
	cin >> num_of_val; // считываем количество значений 
	cout << "Enter number of experiments\n";
	cin >> num_of_experiments; // считываем количество экспериментов 

	fout << num_of_arg << endl; // записываем в файл количество входных значений 
	fout << num_of_val << endl; // записываем в файл количество значений, которые в будущем должны быть угаданы 
	fout << num_of_experiments << endl;// записываем в файл количсетво экспериментов

	double* arr_of_arg = new double[num_of_arg]; // массив для аргументов для одного эксперимента 
	range* arg_range = new range[num_of_arg]; // массив диапазонов значений для аргументов 
	double* arr_of_val = new double[num_of_val]; // массив для значений функции одного эксперимента 

	for (int i = 0; i < num_of_arg; ++i)
	{
		cout << "Enter value range for " << i + 1 << " argument\n"; // считываем диапазон значений для i-ого аргумента 
		cin >> arg_range[i];
	}

	for (int i = 0; i < num_of_experiments; ++i)
	{
		cout << "Experiment number " << i + 1 << endl;
		for (int j = 0; j < num_of_arg; ++j)
		{
			arr_of_arg[j] = arg_range[j].beg + (arg_range[j].end - arg_range[j].beg) * (double)i / num_of_experiments; // записываем значения аргументов для данного эксперимента
			cout << "Argument number " << j + 1 << " = " << arr_of_arg[j] << endl;
			fout << arr_of_arg[j] << endl; // записываем в файл значение j-ого аргумента i-ого эксперимента 
		}
		cout << endl;
		function(arr_of_arg, arr_of_val); // инициализируем массив значений функции
		for (int j = 0; j < num_of_val; ++j)
		{
			cout << j+1 << " value is " << arr_of_val[j] << endl;
			fout << arr_of_val[j] << endl; // записываем в файл j-ое значение i-ого эксперимента 
		}
		cout << endl;
	}
	fout.close();
	return 0;
}