#include <chrono>
using namespace std;
class Timer // класс дл€ измерени€ времени работы программы
{
	using clock_t = chrono::high_resolution_clock;
	using second_t = chrono::duration<double, ratio<1> >;
	chrono::time_point<clock_t> m_beg;
public:
	Timer() : m_beg(clock_t::now()) {}
	double elapsed() const
	{
		return std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
	}
};