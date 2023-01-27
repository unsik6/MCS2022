#include<iostream>
#include<string>
#include <map>
#include <vector>
#include <utility>
#include<time.h>
#include <cstdarg>	// Variable argument

using namespace std;

// BOOLARRAY for iterating combinations
class BOOLARRAY
{
private:
	bool* m_arr;
	int m_size;

public:
	BOOLARRAY(int _size = 1) : m_size(_size)
	{
		m_arr = new bool[m_size]{ false };
	}
	BOOLARRAY(int _size, bool ...) : m_size(_size)
	{
		m_arr = new bool[m_size];
		va_list args;
		va_start(args, _size);

		for (int i = 0; i < m_size; i++)
		{
			m_arr[i] = va_arg(args, bool);
		}

		va_end(args);
	}
	~BOOLARRAY()
	{
		delete[] m_arr;
	}
	
	// ERROR
	bool pauseOUTofRANGE(int _idx)
	{
		if (_idx >= m_size || _idx < 0)
		{
			cerr << "Out of range\n";
			system("PAUSE");
			return true;
		}
		return false;
	}

	// getter/setter
	int getSize() { return m_size; }
	bool getElem(const int _idx)
	{
		if(pauseOUTofRANGE(_idx))
			return false;

		return m_arr[_idx];
	}
	void setElem(const int _idx, const bool _elem)
	{
		if (!pauseOUTofRANGE(_idx))
			m_arr[_idx] = _elem;
	}

	bool isFull()
	{
		for (int i = 0; i < m_size; i++)
			if (!m_arr[i])
				return false;
		return true;
	}
	void switchElem(const int _idx)
	{
		if (!pauseOUTofRANGE(_idx))
			m_arr[_idx] = !m_arr[_idx];
	}

	void next_combination()
	{
		if (isFull())
		{
			cerr << "No more Combination\n";
			return;
		}

		bool carry = true;

		for (int i = 0; i < m_size; i++)
		{
			if (!m_arr[i])
			{
				m_arr[i] = carry;
				return;
			}
			else if (!carry) return;
			else
			{
				m_arr[i] = false;
			}
		}
	}

	friend ostream& operator<< (ostream& _os, const BOOLARRAY& _BA)
	{
		for (int i = 0; i < _BA.m_size; i++)
		{
			_os << _BA.m_arr[i];
		}
		return _os;
	}

};

// Return array of the expected number of common characters about all k.
// In return array A, A[i] = E , where k = i.
double* Naive_Prob_Cal(const double const* _charArr, const int _arrSize, const int _maxK = 512)
{
	// calculate all expect about 1 <= k <= _maxK
	for (int k = 1; k <= _maxK; k++)
	{


		// In Naive solution, check all case from front.
		for (int i = 1; i <= k; i++)
		{

		}
	}
	return 0;
}



int main()
{
	BOOLARRAY a(5, 1, 0, 0, 0, 0);
	cout << a;
	cout << a.getElem(0);
	cout << a.getElem(-1);
	cout << a.getElem(4);
	a.switchElem(0);
	cout << a.getElem(0);
	cout << endl;
	cout << a << '\n';
	for (int i = 0; i < 35; i++)
	{
		a.next_combination();
		cout << a << '\n';
	}
	return 0;
}
