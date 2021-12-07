#include <math.h>
inline double activation_function(double in) // функция активации
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
struct Input_Neuron : Neuron  // структура для входных нейронов
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
	Synapse(Neuron& Out_n, Neuron& In_n, double val)
	{
		Out = &Out_n;
		In = &In_n;
		value = val;
	}
};