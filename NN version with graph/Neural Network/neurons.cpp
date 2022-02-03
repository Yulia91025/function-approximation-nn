#include <math.h>
inline double activation_function(double in) // ������� ���������
{
	return (2 / (1 + exp(-2 * in)) - 1); // ��������������� �������
}

//inline double activation_function(double in) // ������� ���������
//{
//	return (1 / (1 + exp(-in))); // ������� 
//}
struct Neuron // ��������� ��� ������� �������� � �������� ������
{
	double input;
	double output;
	double delta = 0; // ���������� ��� ����������� ����� �������� ��������, ��������� � ���� �������� 
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
struct Input_Neuron : Neuron  // ��������� ��� ������� ��������
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
struct Synapse // ��������� ��� ��������
{
	Neuron* Out; // ��������� �� ������, �� �������� ������ ������� 
	Neuron* In; // ��������� �� ������, � ������� ������ ������ 
	double value; // �������� �������
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