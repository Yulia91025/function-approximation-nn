#include <iostream>
#include "neural_network.cpp"
#include <fstream>
#include <string>
using namespace std;
int main()
{
	string filename; // имя файла
	cout << "Enter filename:" << endl;
	cin >> filename; // считываем имя файла
	filename += ".txt";
	ifstream file(filename);
	// Если мы не можем открыть файл для чтения его содержимого, то выводим сообщение об ошибке
	if (!file)
	{
		cerr << "This file could not be opened for reading!" << endl;
		char c;
		cout << "Please enter any symbol to complete program" << endl;
		cin >> c; // чтобы программа не закрывалась сразу после вывода ошибки
		return -1;
	}
	DATA data(file); // считываем данные 
	Neural_network Net(data); // создаём нейросеть 
	return 0;
};
