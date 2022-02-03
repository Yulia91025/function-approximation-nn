#include <iostream>
#include "neural_network.cpp"
#include <fstream>
#include <string>
#include <GL/freeglut.h>
#include <sstream>
using namespace std;

#include <math.h>

double** coord; // массив координат
string filename; // имя файла

const float g_fieldWidth = 4 * 5 ; // ширина поля
const float g_fieldHeight = 4 * 3; // высота поля
const float g_halfFieldWidth = g_fieldWidth / 2; // половина ширины поля
const float g_halfFieldHeight = g_fieldHeight / 2; // половина высоты поля 
const float g_onePercentOfHalfWidth = (float)g_halfFieldWidth / 100; // один процент от половины ширины поля
const float g_onePercentOfHalfHeight = (float)g_halfFieldHeight / 100; // один процент от половины высоты поля 

void drawText(float x, float y, std::string text, float r, float g, float b)
{
    glColor3f(r, g, b);
    glRasterPos2f(x, y);
    glutBitmapString(GLUT_BITMAP_8_BY_13, (const unsigned char*)text.c_str());
}

void drawLine(
    float x0, float y0, float x1, float y1,
    float r, float g, float b)
{
    glColor3f(r, g, b);
    glBegin(GL_LINES);
    {
        glVertex2f(x0, y0);
        glVertex2f(x1, y1);
    }
    glEnd();
}

void drawGraph(float r, float g, float b)
{
    double xPrev = 0;
    double yPrev = 0;
    for (int i = 0; i < num_in_neurons - 1; ++i) // опа, а тут должен быть выбор по каким переменным график строить !!!!!!!!!
    {
        cout << coord[0][i] << endl;
        xPrev = coord[0][i];
    }
    for (int i = 0; i < num_in_neurons + num_out_neurons - 1; ++i)
    {
        cout << coord[0][i] << endl;
        yPrev = coord[0][i];
    }

    for (int n = 1; n < num_of_experiments; ++n)
    {
        double x0 = xPrev;
        double y0 = yPrev;
        double x1;
        double y1;
        for (int i = 0; i < num_in_neurons - 1; ++i) // опа, а тут должен быть выбор по каким переменным график строить !!!!!!!!!
        {
            cout << coord[n][i] << endl;
            x1 = coord[n][i];
            xPrev = x1;
        }
        for (int i = 0; i < num_in_neurons + num_out_neurons - 1; ++i)
        {
            cout << coord[n][i] << endl;
            y1 = coord[n][i];
            yPrev = y1;
        }

        drawLine(x0, y0, x1, y1, r, g, b); // отрисовываем линию
    }
}

void draw()
{
    glClear(GL_COLOR_BUFFER_BIT);
    ostringstream ostr; // используем строковый поток, чтобы выводить точность на график
    ostr << accuracy; 
    string str = ostr.str();
    drawText(
        g_halfFieldWidth - 50 * g_onePercentOfHalfWidth,
        g_halfFieldHeight - 13 * g_onePercentOfHalfHeight,
        filename, 1, 1, 1); // записываем имя файла
    drawText(
        g_halfFieldWidth - 50 * g_onePercentOfHalfWidth,
        g_halfFieldHeight - 20 * g_onePercentOfHalfHeight,
        str, 1, 1, 1); // записываем точность 

    // X
    drawLine(-g_halfFieldWidth, 0, g_halfFieldWidth, 0, 1, 0, 0);
    drawText(
        g_halfFieldWidth - 10 * g_onePercentOfHalfWidth,
        -10 * g_onePercentOfHalfHeight, "X", 1, 0, 0);

    // Y
    drawLine(0, -g_halfFieldHeight, 0, g_halfFieldHeight, 0, 1, 0);
    drawText(
        -10 * g_onePercentOfHalfWidth,
        g_halfFieldHeight - 15 * g_onePercentOfHalfHeight, "Y", 0, 1, 0);

    //drawCoordinates(1, 2, 1);

    drawGraph(1, 0.5, 0.5); // отрисовываем график - цифры это rgb

    glutSwapBuffers();
}

void setup2DGraphics(double width, double height)
{
    double halfWidth = width / 2;
    double halfHeight = height / 2;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-halfWidth, halfWidth, -halfHeight, halfHeight, 100, -100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

int main(int argc, char** argv)
{
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

    coord = Net.for_graph(data);

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(700, 500); // задаём размеры окна
    glutCreateWindow("Graph"); // название окна 
    setup2DGraphics(g_fieldWidth, g_fieldHeight);
    glutDisplayFunc(draw);
    glutMainLoop();
	return 0;
};