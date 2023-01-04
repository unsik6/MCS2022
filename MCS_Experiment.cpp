#include <iostream>
#include <fstream>
#include <string>
#include <stack>
#include <queue>

#include <algorithm>
#include <vector>
#include <utility>
#include <map>

#include <boost/random.hpp>
#include <time.h>
#include <math.h>

using namespace std;

// The data type to store datas including meta datas.
struct returnPacket
{
	string str1;
	string str2;
	string commsubseq;

	int n;
	int cssLen;
	bool isMaximal;

	double time;
};

returnPacket LCS(string _str1, string _str2);

// The Maximal Common Subsequence algorithms of Yoshifumi Sakai
// Sakai, Y. (2019). Maximal common subsequence algorithms. Theoretical Computer Science, 793, 132-139.
// https://www.sciencedirect.com/science/article/pii/S0304397519304074
struct ghs
{
	int g = -1;
	int h = -1;
	int s = -1;
};
returnPacket MCS_0(string _str1, string _str2);

// The Maximal Common Subsequence algorithms of DongYeop Lee and Joong Chae Na.
// 이동엽, & 나중채. (2022). 더 긴 극대 공통 부분 서열을 찾기 위한 알고리즘. 정보과학회논문지, 49(7), 507-513.
// DongYeop Lee, & Joong Chae Na. (2022), An Improved Algorithm for Finding a Longer Maximal Common Subsequence. Journal of KIISE, 49(7), 507-513.
// https://www.dbpia.co.kr/pdf/pdfView.do?nodeId=NODE11100316&googleIPSandBox=false&mark=0&useDate=&ipRange=false&accessgl=Y&language=ko_KR&hasTopBanner=true
struct leeGhs
{
	int g = 0;
	int h = 0;
	int gs = 0;
	int hs = 0;
};
returnPacket MCS_1(string _str1, string _str2);

// Shin & Sim's
struct myGhs
{
	int g = -1;
	int h = -1;

	bool thisTurn = true; // X = true; Y = false;

	int xS = -1;
	int yS = -1;
};

struct MyGhs
{
	int id = 0;
	int g = 0;
	int h = 0;
	int s = 0;
};

// k Match MCS (using Queue) - Fail
// - Store pending common characters using queues indexed by each common character.
returnPacket MCS_T0(string _str1, string _str2, int _k);

// k match compared MCS (using Deque) - Fail
// - Compare two pending common characters in a deque.
returnPacket MCS_T1(string _str1, string _str2, int _k);

// k match Linear Max Search MCS (using vector)
// The Maximal Common Subsequence algorithms of Hyeonjun Shin and Jeong Seop Sim.
// 신현준, & 심정섭. (2022). 두 문자열의 극대공통부분서열을 찾는 새로운 알고리즘. 2022년 한국소프트웨어종합학술대회 논문집, 1212-1214.
// Hyeonjun Shin, & Jeon Seop Sim. A New Algorithm of Finding a Maximal Common Subsequence of Two Strings. Korea Software Congress 2022, 1212-1214.
struct PQGhs
{
	int id = -1;
	int g = -1;
	int h = -1;
	int s = 0;
};
returnPacket MCS_T2(string _str1, string _str2, int _k);

// k match Two Heap MCS (usin Heap) - Fail
// - using two heaps storing common characters in each string.
struct PQGhs6
{
	int id = -1;
	int g = -1;
	int h = -1;
	int s = 0;

	// push 당시의 iEnd6, jEnd6
	int pushed_iEnd = -1;
	int pushed_jEnd = -1;

	// check 용
	long long pushed_num = -1;
};
returnPacket MCS_T3(string _str1, string _str2, int _k);

bool isMCS(string _str1, string _str2, string _mcs);

// wstring
struct wreturnPacket
{
	wstring str1;
	wstring str2;
	wstring commsubseq;

	int n;
	int cssLen;
	bool isMaximal;

	double time;
};
wreturnPacket LCS(wstring _str1, wstring _str2);
wreturnPacket MCS_0(wstring _str1, wstring _str2);
wreturnPacket MCS_1(wstring _str1, wstring _str2);
wreturnPacket MCS_T2(wstring _str1, wstring _str2, int _k);
bool isMCS(wstring _str1, wstring _str2, wstring _mcs);

// experiment
int funcNOs[4] = { 0, 1, 2, 3 };
int randomAlphabetSizes[11] = { 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048 };
string strX, strY; wstring wstrX, wstrY;
bool ostreamCheck(ofstream& _os, string _name);
bool istreamCheck(ifstream& _is, string _name);
void recordColumn();
void record(ofstream& _os, returnPacket& _rp);
void record(ofstream& _os, wreturnPacket& _rp);
returnPacket experiment_function(int _funcNo, int _k, string _str1, string _str2);
void experiment_RealData(ifstream& _isData, ofstream& _osRecorder, string _dataName, int _sizeIterNum, int _sizeOffset, int _testIterNum);
void experiment_RandomData(int _alphabetSize, ofstream& _osRecorder, int _sizeIterNum, int _sizeOffset, int _testIterNum);

int main()
{
	// DNA
	string filename = "dna.txt";
	ifstream isDNA(filename, ios::in);
	if (!istreamCheck(isDNA, filename)) return -1;
	filename = "dnaResult.txt";
	ofstream osDNA(filename, ios::out);
	if(!ostreamCheck(osDNA, filename)) return -1;

	experiment_RealData(isDNA, osDNA, "DNA", 10, 1000, 100);

	isDNA.clear(); isDNA.close();
	osDNA.clear(); osDNA.close();

	// Protein
	filename = "proteins.txt";
	ifstream isProtein(filename, ios::in);
	if(!istreamCheck(isProtein, filename)) return -1;
	filename = "proteinResult.txt";
	ofstream osProtein(filename, ios::out);
	if (!ostreamCheck(osProtein, filename)) return -1;

	experiment_RealData(isProtein, osProtein, "Protein", 10, 1000, 100);

	isProtein.clear(); isProtein.close();
	osProtein.clear(); osProtein.close();

	// English
	filename = "english.txt";
	ifstream isEng(filename, ios::in);
	if(!istreamCheck(isEng, filename)) return -1;
	filename = "englishResult.txt";
	ofstream osEng(filename, ios::out);
	if(!ostreamCheck(osEng, filename)) return -1;

	experiment_RealData(isEng, osEng, "English", 10, 1000, 100);

	isEng.clear(); isEng.close();
	osEng.clear(); osEng.close();

	// random
	string prefix = "_alphabetResult.txt";
	for (int alphaIdx = 0; alphaIdx < 11; alphaIdx++)
	{
		filename = to_string(randomAlphabetSizes[alphaIdx]) + prefix;
		ofstream osRandom(filename, ios::out);
		if(!ostreamCheck(osRandom, filename)) return -1;

		experiment_RandomData(randomAlphabetSizes[alphaIdx], osRandom, 10, 1000, 100);

		osRandom.clear(); osRandom.close();
	}
	
	return 0;
}

returnPacket LCS(string _str1, string _str2)
{
	// 시간 측정
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	int LCSlength = 0, max;

	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();


	// matrix 생성하는 반복문 최소화를 위한 작업; 더 짧은 것을 str1으로 고정
	if (len1 > len2)
	{
		string tmp = _str1; _str1 = _str2; _str2 = tmp;
		int tmp2 = len1; len1 = len2; len2 = tmp2;
	}

	// LCS 알고리즘을 위해 필요한 행렬
	int** matrix = new int* [len1 + 1];
	for (int i = 0; i < len1 + 1; i++)
	{
		matrix[i] = new int[len2 + 1];
		matrix[i][0] = 0;
	}
	for (int j = 0; j < len2 + 1; j++)
	{
		matrix[0][j] = 0;
	}

	// Find LCS Length
	for (int i = 1; i < len1 + 1; i++)
	{
		max = 0;
		for (int j = 1; j < len2 + 1; j++)
		{
			if (_str1[i - 1] == _str2[j - 1])
			{
				max = matrix[i - 1][j - 1] + 1;
				matrix[i][j] = max;
			}
			else if (matrix[i][j - 1] >= matrix[i - 1][j])
				matrix[i][j] = matrix[i][j - 1];
			else
				matrix[i][j] = matrix[i - 1][j];
		}
		if (LCSlength < max) LCSlength = max;
	}

	// Find LCS
	string lcs = "";
	int curLength = LCSlength, preLength = LCSlength - 1, idxOfLCS = LCSlength - 1;
	int idxJ = len2;

	for (int i = len1; i > 0; i--)
	{
		for (int j = idxJ; j > 0; j--)
		{
			if (matrix[i][j] == curLength
				&& matrix[i - 1][j] == preLength
				&& matrix[i][j - 1] == preLength
				&& matrix[i - 1][j - 1] == preLength)
			{
				curLength--;
				preLength--;
				lcs = _str1[i - 1] + lcs;
				idxJ = j;
				break;
			}
		}
	}

	// 시간 측정
	end_t = clock();
	result_t = (double)(end_t - start_t);

	// 소멸
	for (int i = 0; i < len1 + 1; i++)
	{
		delete[] matrix[i];
	}
	delete[] matrix;

	returnPacket returnPack = { _str1, _str2, lcs, _str1.length() + _str1.length(), lcs.length(), isMCS(_str1, _str2, lcs), result_t };
	return returnPack;
}

// Y Sakai
returnPacket MCS_0(string _str1, string _str2)
{
	// 시간 측정
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();

	string mcs = "";

	stack<ghs> stk;
	ghs init{ -1, -1, 0 };
	stk.push(init);

	// 범위의 오른쪽 끝
	int iEnd = len1; int jEnd = len2;

	while (!stk.empty())
	{
		ghs curGHS = stk.top();
		stk.pop();

		// s가 짝수이면 str1(X)에서 찾는 중이다.
		if (curGHS.s % 2 == 0)
		{
			int searching_g = curGHS.g + curGHS.s / 2 + 1;

			// 아직 str1(X)에서의 탐색 범위가 지정된 범위 내에 존재한다.
			if (searching_g < iEnd)
			{
				ghs next_GHS{ curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);

				// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
				// 새로운 공통문자 offset을 저장
				if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd - 1)
				{
					ghs next_GHS{ searching_g, finding_h, 0 };
					stk.push(next_GHS);
				}
			}
			// str1(X)에서 더이상 탐색할 것이 없고 현재의 공통문자가 시작문자(센티널)이 아닌 경우
			else if (curGHS.g >= 0)
			{
				mcs = _str1[curGHS.g] + mcs;

				iEnd = _str1.find_last_of(_str1[curGHS.g], iEnd - 1);
				jEnd = _str2.find_last_of(_str1[curGHS.g], jEnd - 1);
				if (iEnd == 0 || jEnd == 0 ||
					iEnd == -1 || jEnd == -1 ||
					iEnd == _str1.npos || jEnd == _str2.npos)
					break;
			}
		}
		// s가 홀수이면 str2(Y)에서 찾는 중이다.
		else
		{
			int searching_h = curGHS.h + (curGHS.s + 1) / 2;

			// 아직 str2(Y)에서의 탐색 범위가 지정된 범위 내에 존재한다.
			if (searching_h < jEnd)
			{
				ghs next_GHS{ curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

				// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
				// 새로운 공통문자 offset을 저장
				if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd - 1)
				{
					ghs next_GHS{ finding_g, searching_h, 0 };
					stk.push(next_GHS);
				}
			}
			// str1(X)에서 더이상 탐색할 것이 없고 현재의 공통문자가 시작문자(센티널)이 아닌 경우
			else if (curGHS.g >= 0)
			{
				mcs = _str1[curGHS.g] + mcs;

				iEnd = _str1.find_last_of(_str1[curGHS.g], iEnd - 1);
				jEnd = _str2.find_last_of(_str1[curGHS.g], jEnd - 1);
				if (iEnd == 0 || jEnd == 0 ||
					iEnd == -1 || jEnd == -1 ||
					iEnd == _str1.npos || jEnd == _str2.npos)
					break;
			}
		}
	}

	/*cout << "MCS Sakai Lenght : " << mcs.length() << '\n';
	cout << "MCS Sakai : ";
	cout << (mcs + '\0') << '\n';*/

	// 시간 측정
	end_t = clock();
	result_t = (double)(end_t - start_t);

	returnPacket result = { _str1, _str2, mcs, _str1.length() + _str2.length(), mcs.length(), isMCS(_str1, _str2, mcs), result_t };
	return result;
}

// Lee
returnPacket MCS_1(string _str1, string _str2)
{
	// 시간 측정
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();

	string mcs = "";

	stack<leeGhs> stk;
	leeGhs init{ -1, -1, 1, 1 };
	stk.push(init);

	// 범위의 오른쪽 끝
	int iEnd = len1; int jEnd = len2;

	while (!stk.empty())
	{
		leeGhs curGHS = stk.top();
		stk.pop();

		// X와 Y에서 탐색할 문자의 위치를 구한다.
		int searching_g = curGHS.g + curGHS.gs;
		int searching_h = curGHS.h + curGHS.hs;

		// 두 곳 모두 탐색 범위를 넘지 않는 경우
		if (searching_g < iEnd && searching_h < jEnd)
		{
			int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);
			int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

			// 양쪽 기준의 문자가 다른 문자열에서 나타난 경우
			if ((finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd - 1) &&
				(finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd - 1))
			{
				int minLen_by_g = min((iEnd - searching_g), (jEnd - finding_h));
				int minLen_by_h = min((iEnd - finding_g), (jEnd - searching_h));

				// X에서 찾은 문자를 기준으로 할 때 범위가 더 넓어지는 경우
				if (minLen_by_g >= minLen_by_h)
				{
					leeGhs reGHS{ curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs };
					stk.push(reGHS);

					leeGhs next_GHS{ searching_g, finding_h, 1, 1 };
					stk.push(next_GHS);
				}
				// Y에서 찾은 문자를 기준으로 할 때 범위가 더 넓어지는 경우
				else
				{
					leeGhs reGHS{ curGHS.g, curGHS.h, curGHS.gs, curGHS.hs + 1 };
					stk.push(reGHS);

					leeGhs next_GHS{ finding_g, searching_h, 1, 1 };
					stk.push(next_GHS);
				}
			}
			// X 기준의 문자가 str2에서 나타난 경우만
			else if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd - 1)
			{
				leeGhs reGHS{ curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs + 1 };
				stk.push(reGHS);

				leeGhs next_GHS{ searching_g, finding_h, 1, 1 };
				stk.push(next_GHS);
			}
			// Y 기준의 문자가 str1에서 나타난 경우만
			else if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd - 1)
			{
				leeGhs reGHS{ curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs + 1 };
				stk.push(reGHS);

				leeGhs next_GHS{ finding_g, searching_h, 1, 1 };
				stk.push(next_GHS);
			}
			// 둘다 안 나타난 경우
			else
			{
				leeGhs reGHS{ curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs + 1 };
				stk.push(reGHS);
			}
		}
		else if (curGHS.g >= 0)
		{
			mcs = _str1[curGHS.g] + mcs;

			iEnd = _str1.find_last_of(_str1[curGHS.g], iEnd - 1);
			jEnd = _str2.find_last_of(_str1[curGHS.g], jEnd - 1);
			if (iEnd == 0 || jEnd == 0 ||
				iEnd == -1 || jEnd == -1 ||
				iEnd == _str1.npos || jEnd == _str2.npos)
				break;
		}


	}

	/*cout << "MCS Lee Lenght : " << mcs.length() << '\n';
	cout << "MCS Lee : ";
	cout << (mcs + '\0') << '\n';*/

	// 시간 측정
	end_t = clock();
	result_t = (double)(end_t - start_t);

	returnPacket result = { _str1, _str2, mcs, _str1.length() + _str2.length(), mcs.length(), isMCS(_str1, _str2, mcs), result_t };
	return result;
}

// experiment MCS
// k Match MCS(using Queue) - Fail
returnPacket MCS_T0(string _str1, string _str2, int _k)
{
	// 시간 측정
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();

	string mcs = "";

	stack<MyGhs> stk;
	MyGhs init{ 0, -1, -1, 0 };
	stk.push(init);

	vector<queue<MyGhs>> wthQs;
	queue<MyGhs> newQue;
	wthQs.push_back(newQue);

	// 범위의 오른쪽 끝
	int iEnd = len1; int jEnd = len2;

	while (!stk.empty())
	{

		MyGhs curGHS = stk.top();

		bool isPoped = false;

		// (1) _k만큼 탐색을 했으나 _k 보다 작은 범위 내에서 공통 문자를 찾지 못한 경우
		if (curGHS.s >= _k * 2)
		{
			if (wthQs.size() > curGHS.id)
			{
				while (!wthQs[curGHS.id].empty())
				{
					MyGhs wthGhs = wthQs[curGHS.id].front();
					wthQs[curGHS.id].pop();
					if (wthGhs.g < iEnd && wthGhs.h < jEnd)
					{
						MyGhs next_GHS{ curGHS.id + 1, wthGhs.g, wthGhs.h, 0 };
						stk.push(next_GHS);
						isPoped = true;
						break;
					}
				}
			}
		}
		// (2) X에서 _k만큼 탐색을 진행하기 전에 탐색범위가 종료된 경우
		// (3) Y에서 _k만큼 탐색을 진행하기 전에 탐색범위가 종료된 경우
		// (2)와 (3)의 경우 탐색이 종료될 예정이다. 그러나 탐색이 종료되었을 때, iEnd와 jEnd를 갱신하더라도
		// new_iEnd(new_jEnd)와 old_iEnd(old_jEnd) 사이에 공통문자가 보류된 상태로 남을 수 있다.
		// 즉, 보류 공통문자가 들어가지 않을 수 있으므로 (1)의 경우와 같이 처리해야 한다.
		else if ((curGHS.s % 2 == 0 && curGHS.g + curGHS.s / 2 + 1 >= iEnd) ||
			(curGHS.s % 2 == 1 && curGHS.h + (curGHS.s + 1) / 2 >= jEnd))
		{
			// 보류 공통문자가 없다면 상관 없다.
			if (wthQs.size() > curGHS.id)
			{
				while (!wthQs[curGHS.id].empty())
				{
					MyGhs wthGhs = wthQs[curGHS.id].front();
					wthQs[curGHS.id].pop();
					if (wthGhs.g < iEnd && wthGhs.h < jEnd)
					{
						MyGhs next_GHS{ curGHS.id + 1, wthGhs.g, wthGhs.h, 0 };
						stk.push(next_GHS);
						isPoped = true;
						break;
					}
				}
			}
		}

		if (isPoped) continue;

		stk.pop();

		// s가 짝수이면 str1(X)에서 찾는 중이다.
		if (curGHS.s % 2 == 0)
		{
			int searching_g = curGHS.g + curGHS.s / 2 + 1;

			// 아직 str1(X)에서의 탐색 범위가 지정된 범위 내에 존재한다.
			if (searching_g < iEnd)
			{
				MyGhs next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);

				// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
				// 새로운 공통문자 offset을 저장
				if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd - 1)
				{
					// k 범위 안에서 찾았다.
					if (searching_g <= curGHS.g + _k && finding_h <= curGHS.h + _k)
					{
						MyGhs next_GHS{ curGHS.id + 1, searching_g, finding_h, 0 };
						stk.push(next_GHS);
					}
					// k 범위를 벗어나서 보류한다.
					else
					{
						while (wthQs.size() <= curGHS.id)
						{
							queue<MyGhs> newQue;
							wthQs.push_back(newQue);
						}
						MyGhs wth_GHS{ curGHS.id, searching_g, finding_h, curGHS.s };
						wthQs[curGHS.id].push(wth_GHS);
					}
				}
			}
			// str1(X)에서 더이상 탐색할 것이 없고 현재의 공통문자가 시작문자(센티널)이 아닌 경우
			else if (curGHS.g >= 0)
			{
				mcs = _str1[curGHS.g] + mcs;

				iEnd = _str1.find_last_of(_str1[curGHS.g], iEnd - 1);
				jEnd = _str2.find_last_of(_str1[curGHS.g], jEnd - 1);

				if (iEnd == 0 || jEnd == 0 ||
					iEnd == -1 || jEnd == -1 ||
					iEnd == _str1.npos || jEnd == _str2.npos)
					break;

				// 아직 wthQs[curGHS.id]에 남아 있는 것들은 무조건 g >= iEnd 또는 h >= jEnd일 수밖에 없다.
				// 왜냐하면 아직 남아있다는 것은 축소된 탐색 범위가 _k를 넘지 않아서 빼지 않는 것이고,
				// wthQs[curGHS.id]에 저장되어 있다는 것은 finding_g(h)가 _k범위를 넘기 때문이다.
				// 최악의 경우 _k * 2번 호출하므로 상수 시간이 걸린다.
				if (wthQs.size() > curGHS.id)
					while (!wthQs[curGHS.id].empty())
						wthQs[curGHS.id].pop();
			}
		}
		// s가 홀수이면 str2(Y)에서 찾는 중이다.
		else
		{
			int searching_h = curGHS.h + (curGHS.s + 1) / 2;

			// 아직 str2(Y)에서의 탐색 범위가 지정된 범위 내에 존재한다.
			if (searching_h < jEnd)
			{
				MyGhs next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

				// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
				// 새로운 공통문자 offset을 저장
				if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd - 1)
				{
					// k 범위 안에서 찾았다.
					if (searching_h <= curGHS.h + _k && finding_g <= curGHS.g + _k)
					{
						MyGhs next_GHS{ curGHS.id + 1, finding_g, searching_h, 0 };
						stk.push(next_GHS);
					}
					// k 범위를 벗어나서 보류한다.
					else
					{
						while (wthQs.size() <= curGHS.id)
						{
							queue<MyGhs> newQue;
							wthQs.push_back(newQue);
						}
						MyGhs wth_GHS{ curGHS.id, finding_g, searching_h, curGHS.s };
						wthQs[curGHS.id].push(wth_GHS);
					}

				}
			}
			// str1(X)에서 더이상 탐색할 것이 없고 현재의 공통문자가 시작문자(센티널)이 아닌 경우
			else if (curGHS.g >= 0)
			{
				mcs = _str1[curGHS.g] + mcs;

				if (wthQs.size() > curGHS.id)
					while (!wthQs[curGHS.id].empty())
						wthQs[curGHS.id].pop();

				iEnd = _str1.find_last_of(_str1[curGHS.g], iEnd - 1);
				jEnd = _str2.find_last_of(_str1[curGHS.g], jEnd - 1);
				if (iEnd == 0 || jEnd == 0 ||
					iEnd == -1 || jEnd == -1 ||
					iEnd == _str1.npos || jEnd == _str2.npos)
					break;
			}
		}
	}

	/*cout << "MCS Shin Lenght : " << mcs.length() << '\n';
	cout << "MCS Shin : ";
	cout << (mcs + '\0') << '\n';*/

	// 시간 측정
	end_t = clock();
	result_t = (double)(end_t - start_t);

	returnPacket result = { _str1, _str2, mcs, _str1.length() + _str2.length(), mcs.length(), isMCS(_str1, _str2, mcs), result_t };
	return result;
}

// k match compared MCS (using Deque) - Fail
returnPacket MCS_T1(string _str1, string _str2, int _k)
{
	// 시간 측정
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();

	string mcs = "";

	stack<MyGhs> stk;
	MyGhs init{ 0, -1, -1, 0 };
	stk.push(init);

	vector<deque<MyGhs>> wthQs;
	deque<MyGhs> newQue;
	wthQs.push_back(newQue);

	// 범위의 오른쪽 끝
	int iEnd = len1; int jEnd = len2;

	while (!stk.empty())
	{

		MyGhs curGHS = stk.top();

		bool isPoped = false;

		// (1) _k만큼 탐색을 했으나 _k 보다 작은 범위 내에서 공통 문자를 찾지 못한 경우
		if (curGHS.s >= _k * 2)
		{
			// 보류 공통문자가 없다면 상관 없다.
			if (wthQs.size() > curGHS.id)
			{
				// 가장 앞에 있는 것을 먼저 꺼낸다.
				MyGhs wthGhs = { -1, -1, -1, -1 };
				while (!wthQs[curGHS.id].empty())
				{
					MyGhs tmpGhs = wthQs[curGHS.id].front();
					wthQs[curGHS.id].pop_front();
					if (tmpGhs.g < iEnd && tmpGhs.h < jEnd)
					{
						wthGhs = tmpGhs;
						break;
					}
				}
				// 꺼낸 게 없으면 다음 탐색으로 넘어간다.
				if (wthGhs.id >= 0)
				{
					// 꺼냈음에도 한 개 이상 남아있다면 비교를 위해 꺼낸다.
					MyGhs compared_ghs = { -1, -1, -1, -1 };
					while (!wthQs[curGHS.id].empty())
					{
						MyGhs tmpGhs = wthQs[curGHS.id].front();
						wthQs[curGHS.id].pop_front();
						if (tmpGhs.g < iEnd && tmpGhs.h < jEnd)
						{
							compared_ghs = tmpGhs;
							break;
						}
					}
					// 비교를 위해 꺼낼 것이 없다면 앞서 꺼내둔 것을 사용한다.
					// 만약 비교를 위해 꺼낸 것이 있다면 향후 탐색 범위를 더 크게 만드는 것을 선택하고, 다른 하나는 다시 큐에 집어넣는다.
					if (compared_ghs.id >= 0)
					{
						int minLen_by_wthGHS = min(iEnd - wthGhs.g, jEnd - wthGhs.h);
						int minLen_by_comp = min(iEnd - compared_ghs.g, jEnd - compared_ghs.h);
						if (minLen_by_wthGHS >= minLen_by_comp)
						{
							wthQs[curGHS.id].push_front(compared_ghs);
						}
						else
						{
							wthQs[curGHS.id].push_front(wthGhs);
							wthGhs = compared_ghs;
						}
					}

					MyGhs next_GHS{ curGHS.id + 1, wthGhs.g, wthGhs.h, 0 };
					stk.push(next_GHS);
					isPoped = true;
				}
			}
		}
		// (2) X에서 _k만큼 탐색을 진행하기 전에 탐색범위가 종료된 경우
		// (3) Y에서 _k만큼 탐색을 진행하기 전에 탐색범위가 종료된 경우
		// (2)와 (3)의 경우 탐색이 종료될 예정이다. 그러나 탐색이 종료되었을 때, iEnd와 jEnd를 갱신하더라도
		// new_iEnd(new_jEnd)와 old_iEnd(old_jEnd) 사이에 공통문자가 보류된 상태로 남을 수 있다.
		// 즉, 보류 공통문자가 들어가지 않을 수 있으므로 (1)의 경우와 같이 처리해야 한다.
		else if ((curGHS.s % 2 == 0 && curGHS.g + curGHS.s / 2 + 1 >= iEnd) ||
			(curGHS.s % 2 == 1 && curGHS.h + (curGHS.s + 1) / 2 >= jEnd))
		{
			// 보류 공통문자가 없다면 상관 없다.
			if (wthQs.size() > curGHS.id)
			{
				MyGhs wthGhs = { -1, -1, -1, -1 };
				while (!wthQs[curGHS.id].empty())
				{
					MyGhs tmpGhs = wthQs[curGHS.id].front();
					wthQs[curGHS.id].pop_front();
					if (tmpGhs.g < iEnd && tmpGhs.h < jEnd)
					{
						wthGhs = tmpGhs;
						break;
					}
				}
				if (wthGhs.id >= 0)
				{
					MyGhs compared_ghs = { -1, -1, -1, -1 };
					while (!wthQs[curGHS.id].empty())
					{
						MyGhs tmpGhs = wthQs[curGHS.id].front();
						wthQs[curGHS.id].pop_front();
						if (tmpGhs.g < iEnd && tmpGhs.h < jEnd)
						{
							compared_ghs = tmpGhs;
							break;
						}
					}
					if (compared_ghs.id >= 0)
					{
						int minLen_by_wthGHS = min(iEnd - wthGhs.g, jEnd - wthGhs.h);
						int minLen_by_comp = min(iEnd - compared_ghs.g, jEnd - compared_ghs.h);
						if (minLen_by_wthGHS >= minLen_by_comp)
						{
							wthQs[curGHS.id].push_front(compared_ghs);
						}
						else
						{
							wthQs[curGHS.id].push_front(wthGhs);
							wthGhs = compared_ghs;
						}
					}

					MyGhs next_GHS{ curGHS.id + 1, wthGhs.g, wthGhs.h, 0 };
					stk.push(next_GHS);
					isPoped = true;
				}
			}
		}

		if (isPoped) continue;

		stk.pop();

		// s가 짝수이면 str1(X)에서 찾는 중이다.
		if (curGHS.s % 2 == 0)
		{
			int searching_g = curGHS.g + curGHS.s / 2 + 1;

			// 아직 str1(X)에서의 탐색 범위가 지정된 범위 내에 존재한다.
			if (searching_g < iEnd)
			{
				MyGhs next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);

				// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
				// 새로운 공통문자 offset을 저장
				if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd - 1)
				{
					// k 범위 안에서 찾았다.
					if (searching_g <= curGHS.g + _k && finding_h <= curGHS.h + _k)
					{
						MyGhs next_GHS{ curGHS.id + 1, searching_g, finding_h, 0 };
						stk.push(next_GHS);
					}
					// k 범위를 벗어나서 보류한다.
					else
					{
						while (wthQs.size() <= curGHS.id)
						{
							deque<MyGhs> newQue;
							wthQs.push_back(newQue);
						}
						MyGhs wth_GHS{ curGHS.id, searching_g, finding_h, curGHS.s };
						wthQs[curGHS.id].push_back(wth_GHS);
					}
				}
			}
			// str1(X)에서 더이상 탐색할 것이 없고 현재의 공통문자가 시작문자(센티널)이 아닌 경우
			else if (curGHS.g >= 0)
			{
				mcs = _str1[curGHS.g] + mcs;

				iEnd = _str1.find_last_of(_str1[curGHS.g], iEnd - 1);
				jEnd = _str2.find_last_of(_str1[curGHS.g], jEnd - 1);

				if (iEnd == 0 || jEnd == 0 ||
					iEnd == -1 || jEnd == -1 ||
					iEnd == _str1.npos || jEnd == _str2.npos)
					break;

				// 아직 wthQs[curGHS.id]에 남아 있는 것들은 무조건 g >= iEnd 또는 h >= jEnd일 수밖에 없다.
				// 왜냐하면 아직 남아있다는 것은 축소된 탐색 범위가 _k를 넘지 않아서 빼지 않는 것이고,
				// wthQs[curGHS.id]에 저장되어 있다는 것은 finding_g(h)가 _k범위를 넘기 때문이다.
				// 최악의 경우 _k * 2번 호출하므로 상수 시간이 걸린다.
				if (wthQs.size() > curGHS.id)
					while (!wthQs[curGHS.id].empty())
						wthQs[curGHS.id].pop_front();
			}
		}
		// s가 홀수이면 str2(Y)에서 찾는 중이다.
		else
		{
			int searching_h = curGHS.h + (curGHS.s + 1) / 2;

			// 아직 str2(Y)에서의 탐색 범위가 지정된 범위 내에 존재한다.
			if (searching_h < jEnd)
			{
				MyGhs next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

				// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
				// 새로운 공통문자 offset을 저장
				if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd - 1)
				{
					// k 범위 안에서 찾았다.
					if (searching_h <= curGHS.h + _k && finding_g <= curGHS.g + _k)
					{
						MyGhs next_GHS{ curGHS.id + 1, finding_g, searching_h, 0 };
						stk.push(next_GHS);
					}
					// k 범위를 벗어나서 보류한다.
					else
					{
						while (wthQs.size() <= curGHS.id)
						{
							deque<MyGhs> newQue;
							wthQs.push_back(newQue);
						}
						MyGhs wth_GHS{ curGHS.id, finding_g, searching_h, curGHS.s };
						wthQs[curGHS.id].push_back(wth_GHS);
					}

				}
			}
			// str1(X)에서 더이상 탐색할 것이 없고 현재의 공통문자가 시작문자(센티널)이 아닌 경우
			else if (curGHS.g >= 0)
			{
				mcs = _str1[curGHS.g] + mcs;

				if (wthQs.size() > curGHS.id)
					while (!wthQs[curGHS.id].empty())
						wthQs[curGHS.id].pop_front();

				iEnd = _str1.find_last_of(_str1[curGHS.g], iEnd - 1);
				jEnd = _str2.find_last_of(_str1[curGHS.g], jEnd - 1);
				if (iEnd == 0 || jEnd == 0 ||
					iEnd == -1 || jEnd == -1 ||
					iEnd == _str1.npos || jEnd == _str2.npos)
					break;
			}
		}
	}

	/*cout << "MCS Shin Lenght : " << mcs.length() << '\n';
	cout << "MCS Shin : ";
	cout << (mcs + '\0') << '\n';*/

	// 시간 측정
	end_t = clock();
	result_t = (double)(end_t - start_t);

	returnPacket result = { _str1, _str2, mcs, _str1.length() + _str2.length(), mcs.length(), isMCS(_str1, _str2, mcs), result_t };
	return result;
}

// k match Linear Max Search MCS (using vector)
// The Maximal Common Subsequence algorithms of Hyeonjun Shin and Jeong Seop Sim.
// 신현준, & 심정섭. (2022). 두 문자열의 극대공통부분서열을 찾는 새로운 알고리즘. 2022년 한국소프트웨어종합학술대회 논문집, 1212-1214.
// Hyeonjun Shin, & Jeon Seop Sim. A New Algorithm of Finding a Maximal Common Subsequence of Two Strings. Korea Software Congress 2022, 1212-1214.
int iEnd5, jEnd5;
bool compGHS(const PQGhs _a, const PQGhs _b)
{
	int min_by_a = iEnd5 - _a.g; int max_by_a = jEnd5 - _a.h;
	int min_by_b = iEnd5 - _b.g; int max_by_b = jEnd5 - _b.h;

	if (min_by_a > max_by_a) swap(min_by_a, max_by_a);
	if (min_by_b > max_by_b) swap(min_by_b, max_by_b);

	if (min_by_a > min_by_b)
		return true;
	else if (min_by_a < min_by_b)
		return false;
	else
	{
		return max_by_a >= max_by_b;
	}
}
returnPacket MCS_T2(string _str1, string _str2, int _k)
{
	// 시간 측정
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();

	string mcs = "";

	stack<PQGhs> stk;
	PQGhs init{ 0, -1, -1, 0 };
	stk.push(init);

	vector<vector<PQGhs>> wthPqs;
	vector<PQGhs> init_wthpq;
	wthPqs.push_back(init_wthpq);

	// 범위의 오른쪽 끝
	iEnd5 = len1, jEnd5 = len2;

	// check 용
	int GHS_id = 1;

	while (!stk.empty())
	{
		PQGhs curGHS = stk.top();

		bool isWthPoped = false;

		// (1) _k만큼 탐색을 했으나 _k 보다 작은 범위 내에서 공통 문자를 찾지 못한 경우
		// (2) X에서 _k만큼 탐색을 진행하기 전에 탐색범위가 종료된 경우
		// (3) Y에서 _k만큼 탐색을 진행하기 전에 탐색범위가 종료된 경우
		// (2)와 (3)의 경우 탐색이 종료될 예정이다. 그러나 탐색이 종료되었을 때, iEnd와 jEnd를 갱신하더라도
		// new_iEnd(new_jEnd)와 old_iEnd(old_jEnd) 사이에 공통문자가 보류된 상태로 남을 수 있다.
		// 즉, 보류 공통문자가 들어가지 않을 수 있으므로 (1)의 경우와 같이 처리해야 한다.
		if (curGHS.s >= _k * 2)
		{
			// 보류 공통문자가 없다면 상관 없다.
			if (wthPqs.size() > curGHS.id)
			{
				// 보류 한 것이 없으면 일반 탐색으로 진행
				if (!wthPqs[curGHS.id].empty())
				{
					// 최초 추출 정렬
					// 추출 이후에 엔드 바운드가 달라지지만 그 달라진 엔드 바운드 이후에 기존에 있는 보류 친구들과 엔드 바운드가 다른 보류 친구들이 들어오는 것은 없다. by Lemma 2
					// 가장 탐색 범위가 긴 거 먼저 꺼내기
					// 선형 탐색을 하더라도 엔드 바운드를 넘어버리는 것은 최소 향후 탐색범위가 마이너스이므로 가장 후순위이다.
					PQGhs wthGHS = wthPqs[curGHS.id][0];
					for (int i = 1; i < wthPqs[curGHS.id].size(); i++)
					{
						PQGhs cur = wthPqs[curGHS.id][i];
						if (compGHS(cur, wthGHS))	// cur가 더 길게 만들면 true 반환
						{
							wthGHS = cur;
						}
					}

					// id에 맞는 걸 꺼내므로 End 바운드만 신경쓰면 된다.
					if (wthGHS.g < iEnd5 && wthGHS.h < jEnd5)
					{
						isWthPoped = true;

						// stk에 쌓을 구조체 초기화
						wthGHS.id = GHS_id++;
						wthGHS.s = 0;

						stk.push(wthGHS);
					}
					else
					{
						wthPqs[curGHS.id].clear();
					}
				}
			}
		}
		else if ((curGHS.s % 2 == 0 && curGHS.g + curGHS.s / 2 + 1 >= iEnd5) ||
			(curGHS.s % 2 == 1 && curGHS.h + (curGHS.s + 1) / 2 >= jEnd5))
		{
			// 보류 공통문자가 없다면 상관 없다.
			if (wthPqs.size() > curGHS.id)
			{
				// 보류 한 것이 없으면 일반 탐색으로 진행
				if (!wthPqs[curGHS.id].empty())
				{
					// 최초 추출 정렬
					// 추출 이후에 엔드 바운드가 달라지지만 그 달라진 엔드 바운드 이후에 기존에 있는 보류 친구들과 엔드 바운드가 다른 보류 친구들이 들어오는 것은 없다. by Lemma 2
					// 가장 탐색 범위가 긴 거 먼저 꺼내기
					// 선형 탐색을 하더라도 엔드 바운드를 넘어버리는 것은 최소 향후 탐색범위가 마이너스이므로 가장 후순위이다.
					PQGhs wthGHS = wthPqs[curGHS.id][0];
					for (int i = 1; i < wthPqs[curGHS.id].size(); i++)
					{
						PQGhs cur = wthPqs[curGHS.id][i];
						if (compGHS(cur, wthGHS))	// cur가 더 길게 만들면 true 반환
						{
							wthGHS = cur;
						}
					}

					// id에 맞는 걸 꺼내므로 End 바운드만 신경쓰면 된다.
					if (wthGHS.g < iEnd5 && wthGHS.h < jEnd5)
					{
						isWthPoped = true;

						// stk에 쌓을 구조체 초기화
						wthGHS.id = GHS_id++;
						wthGHS.s = 0;

						stk.push(wthGHS);
					}
					else
					{
						wthPqs[curGHS.id].clear();
					}
				}
			}
		}


		if (isWthPoped)
		{
			curGHS = stk.top();
		}


		stk.pop();

		// s가 짝수이면 str1(X)에서 찾는 중이다.
		if (curGHS.s % 2 == 0)
		{
			int searching_g = curGHS.g + curGHS.s / 2 + 1;

			// 아직 str1(X)에서의 탐색 범위가 지정된 범위 내에 존재한다.
			if (searching_g < iEnd5)
			{
				PQGhs next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);

				// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
				// 새로운 공통문자 offset을 저장
				if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd5 - 1)
				{
					// k 범위 안에서 찾았다.
					if (searching_g <= curGHS.g + _k && finding_h <= curGHS.h + _k)
					{
						PQGhs next_GHS{ GHS_id++, searching_g, finding_h, 0 };
						stk.push(next_GHS);
					}
					// k 범위를 벗어나서 보류한다.
					else
					{
						while (wthPqs.size() <= curGHS.id)
						{
							vector<PQGhs> newQue;
							wthPqs.push_back(newQue);
						}
						PQGhs wth_GHS{ curGHS.id, searching_g, finding_h, curGHS.s };
						wthPqs[curGHS.id].push_back(wth_GHS);
					}
				}
			}
			// str1(X)에서 더이상 탐색할 것이 없고 현재의 공통문자가 시작문자(센티널)이 아닌 경우
			else if (curGHS.g >= 0)
			{
				mcs = _str1[curGHS.g] + mcs;

				iEnd5 = _str1.find_last_of(_str1[curGHS.g], iEnd5 - 1);
				jEnd5 = _str2.find_last_of(_str1[curGHS.g], jEnd5 - 1);
			}
		}
		// s가 홀수이면 str2(Y)에서 찾는 중이다.
		else
		{
			int searching_h = curGHS.h + (curGHS.s + 1) / 2;

			// 아직 str2(Y)에서의 탐색 범위가 지정된 범위 내에 존재한다.
			if (searching_h < jEnd5)
			{
				PQGhs next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

				// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
				// 새로운 공통문자 offset을 저장
				if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd5 - 1)
				{
					// k 범위 안에서 찾았다.
					if (searching_h <= curGHS.h + _k && finding_g <= curGHS.g + _k)
					{
						PQGhs next_GHS{ GHS_id++, finding_g, searching_h, 0 };
						stk.push(next_GHS);
					}
					// k 범위를 벗어나서 보류한다.
					else
					{
						while (wthPqs.size() <= curGHS.id)
						{
							vector<PQGhs> newQue;
							wthPqs.push_back(newQue);
						}
						PQGhs wth_GHS{ curGHS.id, finding_g, searching_h, curGHS.s };
						wthPqs[curGHS.id].push_back(wth_GHS);
					}

				}
			}
			// str1(X)에서 더이상 탐색할 것이 없고 현재의 공통문자가 시작문자(센티널)이 아닌 경우
			else if (curGHS.g >= 0)
			{
				mcs = _str1[curGHS.g] + mcs;

				iEnd5 = _str1.find_last_of(_str1[curGHS.g], iEnd5 - 1);
				jEnd5 = _str2.find_last_of(_str1[curGHS.g], jEnd5 - 1);
			}
		}
	}

	/*cout << "MCS Shin Lenght : " << mcs.length() << '\n';
	cout << "MCS Shin : ";
	cout << (mcs + '\0') << '\n';*/

	// 시간 측정
	end_t = clock();
	result_t = (double)(end_t - start_t);

	returnPacket result = { _str1, _str2, mcs, _str1.length() + _str2.length(), mcs.length(), isMCS(_str1, _str2, mcs), result_t };
	return result;
}

// k match Two Heap MCS (using Heap) - Fail
// - using two heaps storing common characters in each string.
int iEnd6, jEnd6;
struct PQCompX
{
	bool operator()(PQGhs6 _a, PQGhs6 _b)
	{
		int comp_iEnd = (_a.pushed_iEnd > _b.pushed_iEnd ? _a.pushed_iEnd : _b.pushed_iEnd);
		int comp_jEnd = (_a.pushed_jEnd > _b.pushed_jEnd ? _a.pushed_jEnd : _b.pushed_jEnd);
		int minLen_by_a, minLen_by_b;
		int maxLen_by_a, maxLen_by_b;

		if (comp_iEnd - _a.g < comp_jEnd - _a.h)
		{
			minLen_by_a = comp_iEnd - _a.g;
			maxLen_by_a = comp_jEnd - _a.h;
		}
		else
		{
			minLen_by_a = comp_jEnd - _a.h;
			maxLen_by_a = comp_iEnd - _a.g;
		}
		if (comp_iEnd - _b.g < comp_jEnd - _b.h)
		{
			minLen_by_b = comp_iEnd - _b.g;
			maxLen_by_b = comp_jEnd - _b.h;
		}
		else
		{
			minLen_by_b = comp_jEnd - _b.h;
			maxLen_by_b = comp_iEnd - _b.g;
		}

		// 만약 둘다 Y에서 min값을 가질 경우,
		//	추가적인 Lemma에 의하여 X에서의 순서대로 정렬.
		// 그렇지 않을경우, minLen에 의한 정렬로 정렬하면 됨.
		if (minLen_by_a < minLen_by_b)
			return true;
		else if (minLen_by_a > minLen_by_b)
			return false;
		else
		{
			return (_a.g >= _b.g);
		}
	}
};
struct PQCompY
{
	bool operator()(PQGhs6 _a, PQGhs6 _b)
	{
		int comp_iEnd = (_a.pushed_iEnd > _b.pushed_iEnd ? _a.pushed_iEnd : _b.pushed_iEnd);
		int comp_jEnd = (_a.pushed_jEnd > _b.pushed_jEnd ? _a.pushed_jEnd : _b.pushed_jEnd);
		int minLen_by_a, minLen_by_b;
		int maxLen_by_a, maxLen_by_b;

		if (comp_iEnd - _a.g < comp_jEnd - _a.h)
		{
			minLen_by_a = comp_iEnd - _a.g;
			maxLen_by_a = comp_jEnd - _a.h;
		}
		else
		{
			minLen_by_a = comp_jEnd - _a.h;
			maxLen_by_a = comp_iEnd - _a.g;
		}
		if (comp_iEnd - _b.g < comp_jEnd - _b.h)
		{
			minLen_by_b = comp_iEnd - _b.g;
			maxLen_by_b = comp_jEnd - _b.h;
		}
		else
		{
			minLen_by_b = comp_jEnd - _b.h;
			maxLen_by_b = comp_iEnd - _b.g;
		}

		// 만약 둘다 X에서 min값을 가질 경우,
		//	추가적인 Lemma에 의하여 Y에서의 순서대로 정렬.
		// 그렇지 않을경우, minLen에 의한 정렬로 정렬하면 됨.
		if (minLen_by_a < minLen_by_b)
			return true;
		else if (minLen_by_a > minLen_by_b)
			return false;
		else
		{
			return (_a.h >= _b.h);
		}
	}
};
returnPacket MCS_T3(string _str1, string _str2, int _k)
{
	// 시간 측정
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	// 길이 측정
	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();

	string mcs = "";

	stack<PQGhs6> ghStk;
	PQGhs6 init{ 0, -1, -1, 0 };
	ghStk.push(init);
	// End Bound
	iEnd6 = len1; jEnd6 = len2;

	// X에서 발견한 공통문자를 보류하는 경와 Y에서 발견한 공통문자를 보류하는 경우의 힙
	// 각 gh마다 두개의 힙을 저장한다.
	vector<pair< priority_queue<PQGhs6, vector<PQGhs6>, PQCompX >, priority_queue<PQGhs6, vector<PQGhs6>, PQCompY >>> wthTwoHeaps;
	priority_queue<PQGhs6, vector<PQGhs6>, PQCompX > initWthXPQ;
	priority_queue<PQGhs6, vector<PQGhs6>, PQCompY > initWthYPQ;

	wthTwoHeaps.push_back({ initWthXPQ, initWthYPQ });
	// 중복 쿼리를 반복하지 않기 위한 배열
	// 각 배열은 마지막으로 해당 위치에서 쿼리를 진행한 경우를 담고 있다.
	// Worst Case의 경우에만 효과가 있다.
	vector<int> XQueries(len1, -5);
	vector<int> YQueries(len2, -5);

	// check 용
	int GHS_id = 1;
	long long pGHSID = 1;
	PQGhs6 lastGHS;

	while (!ghStk.empty())
	{
		PQGhs6 curGHS = ghStk.top();

		bool isWthPoped = false;

		if (curGHS.s >= _k * 2)
		{
			// 보류 공통문자가 없다면 상관 없다.
			if (wthTwoHeaps.size() > curGHS.id)
			{
				PQGhs6 Xpop, Ypop;
				while (!wthTwoHeaps[curGHS.id].first.empty())
				{
					Xpop = wthTwoHeaps[curGHS.id].first.top();

					if (Xpop.g < iEnd6 && Xpop.h < jEnd6)
					{
						break;
					}
					else
					{
						wthTwoHeaps[curGHS.id].first.pop();
						Xpop.id = -1;
					}
				}
				while (!wthTwoHeaps[curGHS.id].second.empty())
				{
					Ypop = wthTwoHeaps[curGHS.id].second.top();

					if (Ypop.g < iEnd6 && Ypop.h < jEnd6)
					{
						break;
					}
					else
					{
						wthTwoHeaps[curGHS.id].second.pop();
						Ypop.id = -1;
					}
				}
				// 둘다 최대 보류가 나타난 경우
				// 두 보류 공통문자 중 다음 탐색 거리를 더 길게 하는 것을 선택
				if (Xpop.id != -1 && Ypop.id != -1)
				{
					int minLen_by_X = min(iEnd6 - Xpop.g, jEnd6 - Xpop.h);
					int minLen_by_Y = min(iEnd6 - Ypop.g, jEnd6 - Ypop.h);

					// Y가 더 큰 경우 Xpop에 Ypop 넣기 => 최종적으로 Xpop을 스택에 삽입
					if (minLen_by_X < minLen_by_Y)
					{
						Xpop = Ypop;
						wthTwoHeaps[curGHS.id].second.pop();
					}
					else    // X가 더 큰 경우
						wthTwoHeaps[curGHS.id].first.pop();
				}
				// X에서만 최소 보류가 나타난 경우
				else if (Xpop.id != -1)
				{
					wthTwoHeaps[curGHS.id].first.pop();
				}
				// Y에서만 최소 보류가 나타난 경우
				else if (Ypop.id != -1)
				{
					Xpop = Ypop;
					wthTwoHeaps[curGHS.id].second.pop();
				}

				if (Xpop.id != -1)	// 사용 가능한 보류 문자가 존재하는 경우
				{
					isWthPoped = true;

					// ghStk에 쌓을 구조체 초기화
					Xpop.id = GHS_id++;
					Xpop.s = 0;
					Xpop.pushed_num = -2;
					ghStk.push(Xpop);
				}
			}
		}
		else if ((curGHS.s % 2 == 0 && (curGHS.g + curGHS.s / 2 + 1 >= iEnd6 || curGHS.h + (curGHS.s + 2) / 2 >= jEnd6)) ||
			(curGHS.s % 2 == 1 && (curGHS.h + (curGHS.s + 1) / 2 >= jEnd6 || curGHS.g + (curGHS.s + 1) / 2 + 1 >= iEnd6)))
		{
			// 보류 공통문자가 없다면 상관 없다.
			if (wthTwoHeaps.size() > curGHS.id)
			{
				PQGhs6 Xpop, Ypop;
				while (!wthTwoHeaps[curGHS.id].first.empty())
				{
					Xpop = wthTwoHeaps[curGHS.id].first.top();

					if (Xpop.g < iEnd6 && Xpop.h < jEnd6)
					{
						break;
					}
					else
					{
						wthTwoHeaps[curGHS.id].first.pop();
						Xpop.id = -1;
					}
				}
				while (!wthTwoHeaps[curGHS.id].second.empty())
				{
					Ypop = wthTwoHeaps[curGHS.id].second.top();

					if (Ypop.g < iEnd6 && Ypop.h < jEnd6)
					{
						break;
					}
					else
					{
						wthTwoHeaps[curGHS.id].second.pop();
						Ypop.id = -1;
					}
				}
				// 둘다 최대 보류가 나타난 경우
				// 두 보류 공통문자 중 다음 탐색 거리를 더 길게 하는 것을 선택
				if (Xpop.id != -1 && Ypop.id != -1)
				{
					int minLen_by_X = min(iEnd6 - Xpop.g, jEnd6 - Xpop.h);
					int minLen_by_Y = min(iEnd6 - Ypop.g, jEnd6 - Ypop.h);

					// Y가 더 큰 경우 Xpop에 Ypop 넣기 => 최종적으로 Xpop을 스택에 삽입
					if (minLen_by_X < minLen_by_Y)
					{
						Xpop = Ypop;
						wthTwoHeaps[curGHS.id].second.pop();
					}
					else    // X가 더 큰 경우
						wthTwoHeaps[curGHS.id].first.pop();
				}
				// X에서만 최소 보류가 나타난 경우
				else if (Xpop.id != -1)
				{
					wthTwoHeaps[curGHS.id].first.pop();
				}
				// Y에서만 최소 보류가 나타난 경우
				else if (Ypop.id != -1)
				{
					Xpop = Ypop;
					wthTwoHeaps[curGHS.id].second.pop();
				}

				if (Xpop.id != -1)	// 사용 가능한 보류 문자가 존재하는 경우
				{
					isWthPoped = true;

					// ghStk에 쌓을 구조체 초기화
					Xpop.id = GHS_id++;
					Xpop.s = 0;
					Xpop.pushed_num = -2;
					ghStk.push(Xpop);
				}
			}
		}

		if (isWthPoped)
		{
			curGHS = ghStk.top();
		}

		ghStk.pop();

		// s가 짝수이면 str1(X)에서 찾는 중이다.
		if (curGHS.s % 2 == 0)
		{
			int searching_g = curGHS.g + curGHS.s / 2 + 1;

			// 아직 str1(X)에서의 탐색 범위가 지정된 범위 내에 존재한다.
			if (searching_g < iEnd6)
			{
				PQGhs6 next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				ghStk.push(next_GHS);

				int finding_h;
				// 이전에 해두었던 쿼리가 유효한 경우
				if (XQueries[searching_g] < jEnd6 && XQueries[searching_g] > curGHS.h)
				{
					finding_h = XQueries[searching_g];
				}
				else  // 이전에 해두었던 쿼리가 유효하지 않은 경우
				{
					finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);
					XQueries[searching_g] = finding_h;
				}

				// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
				// 새로운 공통문자 offset을 저장
				if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd6 - 1)
				{
					// k 범위 안에서 찾았다.
					if (searching_g <= curGHS.g + _k && finding_h <= curGHS.h + _k)
					{
						PQGhs6 new_GHS{ GHS_id++, searching_g, finding_h, 0 };
						ghStk.push(new_GHS);
					}
					// k 범위를 벗어나서 보류한다.
					else
					{
						while (wthTwoHeaps.size() <= curGHS.id)
						{
							priority_queue<PQGhs6, vector<PQGhs6>, PQCompX > newXPQ;
							priority_queue<PQGhs6, vector<PQGhs6>, PQCompY > newYPQ;

							wthTwoHeaps.push_back({ newXPQ, newYPQ });
						}

						PQGhs6 new_wth_GHS{ curGHS.id, searching_g, finding_h, curGHS.s, iEnd6, jEnd6, pGHSID++ };
						// X에서 찾았으므로
						wthTwoHeaps[curGHS.id].first.push(new_wth_GHS);
					}
				}
			}
			// str1(X)에서 더이상 탐색할 것이 없고 현재의 공통문자가 시작문자(센티널)이 아닌 경우
			else if (curGHS.g >= 0)
			{
				mcs = _str1[curGHS.g] + mcs;
				iEnd6 = _str1.find_last_of(_str1[curGHS.g], iEnd6 - 1);
				jEnd6 = _str2.find_last_of(_str1[curGHS.g], jEnd6 - 1);

				if (iEnd6 == 0 || jEnd6 == 0 ||
					iEnd6 == -1 || jEnd6 == -1 ||
					iEnd6 == _str1.npos || jEnd6 == _str2.npos)
					break;

				lastGHS = curGHS;
				// 아직 wthQs[curGHS.id]에 남아 있는 것들은 무조건 g >= iEnd 또는 h >= jEnd일 수밖에 없다.
				// 왜냐하면 아직 남아있다는 것은 축소된 탐색 범위가 _k를 넘지 않아서 빼지 않는 것이고,
				// wthQs[curGHS.id]에 저장되어 있다는 것은 finding_g(h)가 _k범위를 넘기 때문이다.
				// 최악의 경우 _k * 2번 호출하므로 O(_klog_k) 시간이 걸린다.
				if (wthTwoHeaps.size() > curGHS.id)
				{
					if (!wthTwoHeaps[curGHS.id].first.empty() || !wthTwoHeaps[curGHS.id].second.empty())
					{
						system("PASUE");
					}
					while (!wthTwoHeaps[curGHS.id].first.empty())
						wthTwoHeaps[curGHS.id].first.pop();
					while (!wthTwoHeaps[curGHS.id].second.empty())
						wthTwoHeaps[curGHS.id].second.pop();
				}
			}
		}
		// s가 홀수이면 str2(Y)에서 찾는 중이다.
		else
		{
			int searching_h = curGHS.h + (curGHS.s + 1) / 2;

			// 아직 str2(Y)에서의 탐색 범위가 지정된 범위 내에 존재한다.
			if (searching_h < jEnd6)
			{
				PQGhs6 next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				ghStk.push(next_GHS);

				int finding_g;
				// 이전에 해두었던 쿼리가 유효한 경우
				if (YQueries[searching_h] < iEnd6 && YQueries[searching_h] > curGHS.g)
				{
					finding_g = YQueries[searching_h];
				}
				else  // 이전에 해두었던 쿼리가 유효하지 않은 경우
				{
					finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);
					YQueries[searching_h] = finding_g;
				}
				// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
				// 새로운 공통문자 offset을 저장
				if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd6 - 1)
				{
					// k 범위 안에서 찾았다.
					if (searching_h <= curGHS.h + _k && finding_g <= curGHS.g + _k)
					{
						PQGhs6 new_GHS{ GHS_id++, finding_g, searching_h, 0 };
						ghStk.push(new_GHS);
					}
					// k 범위를 벗어나서 보류한다.
					else
					{
						while (wthTwoHeaps.size() <= curGHS.id)
						{
							priority_queue<PQGhs6, vector<PQGhs6>, PQCompX > newXPQ;
							priority_queue<PQGhs6, vector<PQGhs6>, PQCompY > newYPQ;

							wthTwoHeaps.push_back({ newXPQ, newYPQ });
						}
						PQGhs6 new_wth_GHS{ curGHS.id, finding_g, searching_h, curGHS.s, iEnd6, jEnd6, pGHSID++ };
						// Y에서 찾았으므로
						wthTwoHeaps[curGHS.id].second.push(new_wth_GHS);
					}

				}
			}
			// str1(X)에서 더이상 탐색할 것이 없고 현재의 공통문자가 시작문자(센티널)이 아닌 경우
			else if (curGHS.g >= 0)
			{
				mcs = _str1[curGHS.g] + mcs;
				iEnd6 = _str1.find_last_of(_str1[curGHS.g], iEnd6 - 1);
				jEnd6 = _str2.find_last_of(_str1[curGHS.g], jEnd6 - 1);

				if (curGHS.g != iEnd6 && curGHS.h != jEnd6)
					cout << "IJ wrong" << endl;

				if (iEnd6 == 0 || jEnd6 == 0 ||
					iEnd6 == -1 || jEnd6 == -1 ||
					iEnd6 == _str1.npos || jEnd6 == _str2.npos)
					break;

				lastGHS = curGHS;

				// 아직 wthQs[curGHS.id]에 남아 있는 것들은 무조건 g >= iEnd 또는 h >= jEnd일 수밖에 없다.
				// 왜냐하면 아직 남아있다는 것은 축소된 탐색 범위가 _k를 넘지 않아서 빼지 않는 것이고,
				// wthQs[curGHS.id]에 저장되어 있다는 것은 finding_g(h)가 _k범위를 넘기 때문이다.
				// 최악의 경우 _k * 2번 호출하므로 O(_klog_k) 시간이 걸린다.
				if (wthTwoHeaps.size() > curGHS.id)
				{
					if (!wthTwoHeaps[curGHS.id].first.empty() || !wthTwoHeaps[curGHS.id].second.empty())
					{
						system("PASUE");
					}
					while (!wthTwoHeaps[curGHS.id].first.empty())
						wthTwoHeaps[curGHS.id].first.pop();
					while (!wthTwoHeaps[curGHS.id].second.empty())
						wthTwoHeaps[curGHS.id].second.pop();
				}
			}
		}
	}

	// 시간 측정
	end_t = clock();
	result_t = (double)(end_t - start_t);

	returnPacket result = { _str1, _str2, mcs, _str1.length() + _str2.length(), mcs.length(), isMCS(_str1, _str2, mcs), result_t };
	return result;
}

bool isMCS(string _str1, string _str2, string _mcs)
{
	int g = -1; int h = -1;
	stack<ghs> stk;
	ghs init{ -1, -1 };
	stk.push(init);

	for (int k = 0; k < _mcs.length(); k++)
	{
		g = _str1.find_first_of(_mcs[k], g + 1);
		h = _str2.find_first_of(_mcs[k], h + 1);
		if (g == -1 || g == _str1.npos || h == -1 || h == _str2.npos)
		{
			cout << "isMCS ERROR: There is no character";
			return false;
		}

		ghs pushObj{ g, h, -1 };
		stk.push(pushObj);
	}

	int iEnd = _str1.length(), jEnd = _str2.length();
	while (stk.size() > 1)
	{
		ghs curGH = stk.top();
		stk.pop();

		for (int i = iEnd - 1; i > curGH.g; i--)
		{
			int finding_last_in_Y = _str2.find_last_of(_str1[i], jEnd - 1);
			if (finding_last_in_Y != -1 || finding_last_in_Y != _str2.npos)
			{
				if (finding_last_in_Y > curGH.h)
				{
					cout << "ERROR 3" << endl;
					return false;
				}
			}
		}
		for (int j = jEnd - 1; j > curGH.h; j--)
		{
			int finding_last_in_X = _str1.find_last_of(_str2[j], iEnd - 1);
			if (finding_last_in_X != -1 || finding_last_in_X != _str1.npos)
			{
				if (finding_last_in_X > curGH.g)
				{
					cout << "ERROR 3" << endl;
					return false;
				}
			}
		}
		iEnd = _str1.find_last_of(_str1[curGH.g], iEnd - 1);
		jEnd = _str2.find_last_of(_str2[curGH.h], jEnd - 1);
	}
	return true;
}


// wstring
wreturnPacket LCS(wstring _str1, wstring _str2)
{
	// 시간 측정
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	int LCSlength = 0, max;

	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();


	// matrix 생성하는 반복문 최소화를 위한 작업; 더 짧은 것을 str1으로 고정
	if (len1 > len2)
	{
		wstring tmp = _str1; _str1 = _str2; _str2 = tmp;
		int tmp2 = len1; len1 = len2; len2 = tmp2;
	}

	// LCS 알고리즘을 위해 필요한 행렬
	int** matrix = new int* [len1 + 1];
	for (int i = 0; i < len1 + 1; i++)
	{
		matrix[i] = new int[len2 + 1];
		matrix[i][0] = 0;
	}
	for (int j = 0; j < len2 + 1; j++)
	{
		matrix[0][j] = 0;
	}

	// Find LCS Length
	for (int i = 1; i < len1 + 1; i++)
	{
		max = 0;
		for (int j = 1; j < len2 + 1; j++)
		{
			if (_str1[i - 1] == _str2[j - 1])
			{
				max = matrix[i - 1][j - 1] + 1;
				matrix[i][j] = max;
			}
			else if (matrix[i][j - 1] >= matrix[i - 1][j])
				matrix[i][j] = matrix[i][j - 1];
			else
				matrix[i][j] = matrix[i - 1][j];
		}
		if (LCSlength < max) LCSlength = max;
	}

	// Find LCS
	wstring lcs = L"";
	int curLength = LCSlength, preLength = LCSlength - 1, idxOfLCS = LCSlength - 1;
	int idxJ = len2;

	for (int i = len1; i > 0; i--)
	{
		for (int j = idxJ; j > 0; j--)
		{
			if (matrix[i][j] == curLength
				&& matrix[i - 1][j] == preLength
				&& matrix[i][j - 1] == preLength
				&& matrix[i - 1][j - 1] == preLength)
			{
				curLength--;
				preLength--;
				lcs = _str1[i - 1] + lcs;
				idxJ = j;
				break;
			}
		}
	}

	// 시간 측정
	end_t = clock();
	result_t = (double)(end_t - start_t);

	// 소멸
	for (int i = 0; i < len1 + 1; i++)
	{
		delete[] matrix[i];
	}
	delete[] matrix;

	wreturnPacket returnPack = { _str1, _str2, lcs, _str1.length() + _str1.length(), lcs.length(), isMCS(_str1, _str2, lcs), result_t };
	return returnPack;
}
wreturnPacket MCS_0(wstring _str1, wstring _str2)
{
	// 시간 측정
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();

	wstring mcs = L"";

	stack<ghs> stk;
	ghs init{ -1, -1, 0 };
	stk.push(init);

	// 범위의 오른쪽 끝
	int iEnd = len1; int jEnd = len2;

	while (!stk.empty())
	{
		ghs curGHS = stk.top();
		stk.pop();

		// s가 짝수이면 str1(X)에서 찾는 중이다.
		if (curGHS.s % 2 == 0)
		{
			int searching_g = curGHS.g + curGHS.s / 2 + 1;

			// 아직 str1(X)에서의 탐색 범위가 지정된 범위 내에 존재한다.
			if (searching_g < iEnd)
			{
				ghs next_GHS{ curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);

				// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
				// 새로운 공통문자 offset을 저장
				if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd - 1)
				{
					ghs next_GHS{ searching_g, finding_h, 0 };
					stk.push(next_GHS);
				}
			}
			// str1(X)에서 더이상 탐색할 것이 없고 현재의 공통문자가 시작문자(센티널)이 아닌 경우
			else if (curGHS.g >= 0)
			{
				mcs = _str1[curGHS.g] + mcs;

				iEnd = _str1.find_last_of(_str1[curGHS.g], iEnd - 1);
				jEnd = _str2.find_last_of(_str1[curGHS.g], jEnd - 1);
				if (iEnd == 0 || jEnd == 0 ||
					iEnd == -1 || jEnd == -1 ||
					iEnd == _str1.npos || jEnd == _str2.npos)
					break;
			}
		}
		// s가 홀수이면 str2(Y)에서 찾는 중이다.
		else
		{
			int searching_h = curGHS.h + (curGHS.s + 1) / 2;

			// 아직 str2(Y)에서의 탐색 범위가 지정된 범위 내에 존재한다.
			if (searching_h < jEnd)
			{
				ghs next_GHS{ curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

				// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
				// 새로운 공통문자 offset을 저장
				if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd - 1)
				{
					ghs next_GHS{ finding_g, searching_h, 0 };
					stk.push(next_GHS);
				}
			}
			// str1(X)에서 더이상 탐색할 것이 없고 현재의 공통문자가 시작문자(센티널)이 아닌 경우
			else if (curGHS.g >= 0)
			{
				mcs = _str1[curGHS.g] + mcs;

				iEnd = _str1.find_last_of(_str1[curGHS.g], iEnd - 1);
				jEnd = _str2.find_last_of(_str1[curGHS.g], jEnd - 1);
				if (iEnd == 0 || jEnd == 0 ||
					iEnd == -1 || jEnd == -1 ||
					iEnd == _str1.npos || jEnd == _str2.npos)
					break;
			}
		}
	}

	/*cout << "MCS Sakai Lenght : " << mcs.length() << '\n';
	cout << "MCS Sakai : ";
	cout << (mcs + '\0') << '\n';*/

	// 시간 측정
	end_t = clock();
	result_t = (double)(end_t - start_t);

	wreturnPacket result = { _str1, _str2, mcs, _str1.length() + _str2.length(), mcs.length(), isMCS(_str1, _str2, mcs), result_t };
	return result;
}
wreturnPacket MCS_1(wstring _str1, wstring _str2)
{
	// 시간 측정
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();

	wstring mcs = L"";

	stack<leeGhs> stk;
	leeGhs init{ -1, -1, 1, 1 };
	stk.push(init);

	// 범위의 오른쪽 끝
	int iEnd = len1; int jEnd = len2;

	while (!stk.empty())
	{
		leeGhs curGHS = stk.top();
		stk.pop();

		// X와 Y에서 탐색할 문자의 위치를 구한다.
		int searching_g = curGHS.g + curGHS.gs;
		int searching_h = curGHS.h + curGHS.hs;

		// 두 곳 모두 탐색 범위를 넘지 않는 경우
		if (searching_g < iEnd && searching_h < jEnd)
		{
			int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);
			int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

			// 양쪽 기준의 문자가 다른 문자열에서 나타난 경우
			if ((finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd - 1) &&
				(finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd - 1))
			{
				int minLen_by_g = min((iEnd - searching_g), (jEnd - finding_h));
				int minLen_by_h = min((iEnd - finding_g), (jEnd - searching_h));

				// X에서 찾은 문자를 기준으로 할 때 범위가 더 넓어지는 경우
				if (minLen_by_g >= minLen_by_h)
				{
					leeGhs reGHS{ curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs };
					stk.push(reGHS);

					leeGhs next_GHS{ searching_g, finding_h, 1, 1 };
					stk.push(next_GHS);
				}
				// Y에서 찾은 문자를 기준으로 할 때 범위가 더 넓어지는 경우
				else
				{
					leeGhs reGHS{ curGHS.g, curGHS.h, curGHS.gs, curGHS.hs + 1 };
					stk.push(reGHS);

					leeGhs next_GHS{ finding_g, searching_h, 1, 1 };
					stk.push(next_GHS);
				}
			}
			// X 기준의 문자가 str2에서 나타난 경우만
			else if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd - 1)
			{
				leeGhs reGHS{ curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs + 1 };
				stk.push(reGHS);

				leeGhs next_GHS{ searching_g, finding_h, 1, 1 };
				stk.push(next_GHS);
			}
			// Y 기준의 문자가 str1에서 나타난 경우만
			else if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd - 1)
			{
				leeGhs reGHS{ curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs + 1 };
				stk.push(reGHS);

				leeGhs next_GHS{ finding_g, searching_h, 1, 1 };
				stk.push(next_GHS);
			}
			// 둘다 안 나타난 경우
			else
			{
				leeGhs reGHS{ curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs + 1 };
				stk.push(reGHS);
			}
		}
		else if (curGHS.g >= 0)
		{
			mcs = _str1[curGHS.g] + mcs;

			iEnd = _str1.find_last_of(_str1[curGHS.g], iEnd - 1);
			jEnd = _str2.find_last_of(_str1[curGHS.g], jEnd - 1);
			if (iEnd == 0 || jEnd == 0 ||
				iEnd == -1 || jEnd == -1 ||
				iEnd == _str1.npos || jEnd == _str2.npos)
				break;
		}


	}

	/*cout << "MCS Lee Lenght : " << mcs.length() << '\n';
	cout << "MCS Lee : ";
	cout << (mcs + '\0') << '\n';*/

	// 시간 측정
	end_t = clock();
	result_t = (double)(end_t - start_t);

	wreturnPacket result = { _str1, _str2, mcs, _str1.length() + _str2.length(), mcs.length(), isMCS(_str1, _str2, mcs), result_t };
	return result;
}
wreturnPacket MCS_T2(wstring _str1, wstring _str2, int _k)
{
	// 시간 측정
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();

	wstring mcs = L"";

	stack<PQGhs> stk;
	PQGhs init{ 0, -1, -1, 0 };
	stk.push(init);

	vector<vector<PQGhs>> wthPqs;
	vector<PQGhs> init_wthpq;
	wthPqs.push_back(init_wthpq);

	// 범위의 오른쪽 끝
	iEnd5 = len1, jEnd5 = len2;

	// check 용
	int GHS_id = 1;

	while (!stk.empty())
	{
		PQGhs curGHS = stk.top();

		bool isWthPoped = false;

		// (1) _k만큼 탐색을 했으나 _k 보다 작은 범위 내에서 공통 문자를 찾지 못한 경우
		// (2) X에서 _k만큼 탐색을 진행하기 전에 탐색범위가 종료된 경우
		// (3) Y에서 _k만큼 탐색을 진행하기 전에 탐색범위가 종료된 경우
		// (2)와 (3)의 경우 탐색이 종료될 예정이다. 그러나 탐색이 종료되었을 때, iEnd와 jEnd를 갱신하더라도
		// new_iEnd(new_jEnd)와 old_iEnd(old_jEnd) 사이에 공통문자가 보류된 상태로 남을 수 있다.
		// 즉, 보류 공통문자가 들어가지 않을 수 있으므로 (1)의 경우와 같이 처리해야 한다.
		if (curGHS.s >= _k * 2)
		{
			// 보류 공통문자가 없다면 상관 없다.
			if (wthPqs.size() > curGHS.id)
			{
				// 보류 한 것이 없으면 일반 탐색으로 진행
				if (!wthPqs[curGHS.id].empty())
				{
					// 최초 추출 정렬
					// 추출 이후에 엔드 바운드가 달라지지만 그 달라진 엔드 바운드 이후에 기존에 있는 보류 친구들과 엔드 바운드가 다른 보류 친구들이 들어오는 것은 없다. by Lemma 2
					// 가장 탐색 범위가 긴 거 먼저 꺼내기
					// 선형 탐색을 하더라도 엔드 바운드를 넘어버리는 것은 최소 향후 탐색범위가 마이너스이므로 가장 후순위이다.
					PQGhs wthGHS = wthPqs[curGHS.id][0];
					for (int i = 1; i < wthPqs[curGHS.id].size(); i++)
					{
						PQGhs cur = wthPqs[curGHS.id][i];
						if (compGHS(cur, wthGHS))	// cur가 더 길게 만들면 true 반환
						{
							wthGHS = cur;
						}
					}

					// id에 맞는 걸 꺼내므로 End 바운드만 신경쓰면 된다.
					if (wthGHS.g < iEnd5 && wthGHS.h < jEnd5)
					{
						isWthPoped = true;

						// stk에 쌓을 구조체 초기화
						wthGHS.id = GHS_id++;
						wthGHS.s = 0;

						stk.push(wthGHS);
					}
					else
					{
						wthPqs[curGHS.id].clear();
					}
				}
			}
		}
		else if ((curGHS.s % 2 == 0 && curGHS.g + curGHS.s / 2 + 1 >= iEnd5) ||
			(curGHS.s % 2 == 1 && curGHS.h + (curGHS.s + 1) / 2 >= jEnd5))
		{
			// 보류 공통문자가 없다면 상관 없다.
			if (wthPqs.size() > curGHS.id)
			{
				// 보류 한 것이 없으면 일반 탐색으로 진행
				if (!wthPqs[curGHS.id].empty())
				{
					// 최초 추출 정렬
					// 추출 이후에 엔드 바운드가 달라지지만 그 달라진 엔드 바운드 이후에 기존에 있는 보류 친구들과 엔드 바운드가 다른 보류 친구들이 들어오는 것은 없다. by Lemma 2
					// 가장 탐색 범위가 긴 거 먼저 꺼내기
					// 선형 탐색을 하더라도 엔드 바운드를 넘어버리는 것은 최소 향후 탐색범위가 마이너스이므로 가장 후순위이다.
					PQGhs wthGHS = wthPqs[curGHS.id][0];
					for (int i = 1; i < wthPqs[curGHS.id].size(); i++)
					{
						PQGhs cur = wthPqs[curGHS.id][i];
						if (compGHS(cur, wthGHS))	// cur가 더 길게 만들면 true 반환
						{
							wthGHS = cur;
						}
					}

					// id에 맞는 걸 꺼내므로 End 바운드만 신경쓰면 된다.
					if (wthGHS.g < iEnd5 && wthGHS.h < jEnd5)
					{
						isWthPoped = true;

						// stk에 쌓을 구조체 초기화
						wthGHS.id = GHS_id++;
						wthGHS.s = 0;

						stk.push(wthGHS);
					}
					else
					{
						wthPqs[curGHS.id].clear();
					}
				}
			}
		}


		if (isWthPoped)
		{
			curGHS = stk.top();
		}


		stk.pop();

		// s가 짝수이면 str1(X)에서 찾는 중이다.
		if (curGHS.s % 2 == 0)
		{
			int searching_g = curGHS.g + curGHS.s / 2 + 1;

			// 아직 str1(X)에서의 탐색 범위가 지정된 범위 내에 존재한다.
			if (searching_g < iEnd5)
			{
				PQGhs next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);

				// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
				// 새로운 공통문자 offset을 저장
				if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd5 - 1)
				{
					// k 범위 안에서 찾았다.
					if (searching_g <= curGHS.g + _k && finding_h <= curGHS.h + _k)
					{
						PQGhs next_GHS{ GHS_id++, searching_g, finding_h, 0 };
						stk.push(next_GHS);
					}
					// k 범위를 벗어나서 보류한다.
					else
					{
						while (wthPqs.size() <= curGHS.id)
						{
							vector<PQGhs> newQue;
							wthPqs.push_back(newQue);
						}
						PQGhs wth_GHS{ curGHS.id, searching_g, finding_h, curGHS.s };
						wthPqs[curGHS.id].push_back(wth_GHS);
					}
				}
			}
			// str1(X)에서 더이상 탐색할 것이 없고 현재의 공통문자가 시작문자(센티널)이 아닌 경우
			else if (curGHS.g >= 0)
			{
				mcs = _str1[curGHS.g] + mcs;

				iEnd5 = _str1.find_last_of(_str1[curGHS.g], iEnd5 - 1);
				jEnd5 = _str2.find_last_of(_str1[curGHS.g], jEnd5 - 1);
			}
		}
		// s가 홀수이면 str2(Y)에서 찾는 중이다.
		else
		{
			int searching_h = curGHS.h + (curGHS.s + 1) / 2;

			// 아직 str2(Y)에서의 탐색 범위가 지정된 범위 내에 존재한다.
			if (searching_h < jEnd5)
			{
				PQGhs next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

				// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
				// 새로운 공통문자 offset을 저장
				if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd5 - 1)
				{
					// k 범위 안에서 찾았다.
					if (searching_h <= curGHS.h + _k && finding_g <= curGHS.g + _k)
					{
						PQGhs next_GHS{ GHS_id++, finding_g, searching_h, 0 };
						stk.push(next_GHS);
					}
					// k 범위를 벗어나서 보류한다.
					else
					{
						while (wthPqs.size() <= curGHS.id)
						{
							vector<PQGhs> newQue;
							wthPqs.push_back(newQue);
						}
						PQGhs wth_GHS{ curGHS.id, finding_g, searching_h, curGHS.s };
						wthPqs[curGHS.id].push_back(wth_GHS);
					}

				}
			}
			// str1(X)에서 더이상 탐색할 것이 없고 현재의 공통문자가 시작문자(센티널)이 아닌 경우
			else if (curGHS.g >= 0)
			{
				mcs = _str1[curGHS.g] + mcs;

				iEnd5 = _str1.find_last_of(_str1[curGHS.g], iEnd5 - 1);
				jEnd5 = _str2.find_last_of(_str1[curGHS.g], jEnd5 - 1);
			}
		}
	}

	/*cout << "MCS Shin Lenght : " << mcs.length() << '\n';
	cout << "MCS Shin : ";
	cout << (mcs + '\0') << '\n';*/

	// 시간 측정
	end_t = clock();
	result_t = (double)(end_t - start_t);

	wreturnPacket result = { _str1, _str2, mcs, _str1.length() + _str2.length(), mcs.length(), isMCS(_str1, _str2, mcs), result_t };
	return result;
}
bool isMCS(wstring _str1, wstring _str2, wstring _mcs)
{
	int g = -1; int h = -1;
	stack<ghs> stk;
	ghs init{ -1, -1 };
	stk.push(init);

	for (int k = 0; k < _mcs.length(); k++)
	{
		g = _str1.find_first_of(_mcs[k], g + 1);
		h = _str2.find_first_of(_mcs[k], h + 1);
		if (g == -1 || g == _str1.npos || h == -1 || h == _str2.npos)
		{
			cout << "isMCS ERROR: There is no character";
			return false;
		}

		ghs pushObj{ g, h, -1 };
		stk.push(pushObj);
	}

	int iEnd = _str1.length(), jEnd = _str2.length();
	while (stk.size() > 1)
	{
		ghs curGH = stk.top();
		stk.pop();

		for (int i = iEnd - 1; i > curGH.g; i--)
		{
			int finding_last_in_Y = _str2.find_last_of(_str1[i], jEnd - 1);
			if (finding_last_in_Y != -1 || finding_last_in_Y != _str2.npos)
			{
				if (finding_last_in_Y > curGH.h)
				{
					cout << "ERROR 3" << endl;
					return false;
				}
			}
		}
		for (int j = jEnd - 1; j > curGH.h; j--)
		{
			int finding_last_in_X = _str1.find_last_of(_str2[j], iEnd - 1);
			if (finding_last_in_X != -1 || finding_last_in_X != _str1.npos)
			{
				if (finding_last_in_X > curGH.g)
				{
					cout << "ERROR 3" << endl;
					return false;
				}
			}
		}
		iEnd = _str1.find_last_of(_str1[curGH.g], iEnd - 1);
		jEnd = _str2.find_last_of(_str2[curGH.h], jEnd - 1);
	}
	return true;
}

// experiment
bool ostreamCheck(ofstream& _os, string _name)
{
	if (!_os)
	{
		cerr << _name << ": An error occurred when opening the file" << endl;
		return false;
	}
	else return true;
}
bool istreamCheck(ifstream& _is, string _name)
{
	if (!_is)
	{
		cerr << _name << ": An error occurred when opening the file" << endl;
		return false;
	}
	else return true;
}
void recordColumn(ofstream& _os)
{
	_os << "testNumber\t|X|\t"
		<< "LCS_Len\tcorrectness\tperLen\ttime\t"
		<< "Sakai_Len\tcorrectness\tperLen\ttime\t"
		<< "Lee_Len\tcorrectness\tperLen\ttime\t";

	for (int k = 1; k <= 500; k++)
	{
		_os << k << "_Len\tcorrectness\tperLen\ttime\t";
	}

	_os << '\n';
}
void record(ofstream& _os, returnPacket& _rp)
{
	if (!_rp.isMaximal)
		exit(-1);
	_os << _rp.cssLen << '\t' << _rp.isMaximal << '\t' << (double)((double)_rp.cssLen / (double)_rp.n * 100 * 2) << '\t' << _rp.time << '\t';
}
void record(ofstream& _os, wreturnPacket& _rp)
{
	if (!_rp.isMaximal)
		exit(-1);
	_os << _rp.cssLen << '\t' << _rp.isMaximal << '\t' << (double)((double)_rp.cssLen / (double)_rp.n * 100 * 2) << '\t' << _rp.time << '\t';
}
returnPacket experiment_function(int _funcNo, int _k, string _str1, string _str2)
{
	switch (_funcNo)
	{
	case 0:
		return MCS_T0(_str1, _str2, _k);
	case 1:
		return MCS_T1(_str1, _str2, _k);
	case 2:
		return MCS_T2(_str1, _str2, _k);
	case 3:
		return MCS_T3(_str1, _str2, _k);
	default:
		cerr << "WRONG FUNC\n";
		break;
	}
}
void experiment_RealData(ifstream& _isData, ofstream& _osRecorder, string _dataName, int _sizeIterNum, int _sizeOffset, int _testIterNum)
{
	// 문자열에 대한 문자 집합 분석
	// alphabets[char] = {strX에서의 발생횟수, strY에서의 발생횟수}
	map<char, pair<int, int>> alphabets;

	recordColumn(_osRecorder);

	for (int sizeNum = 0; sizeNum < _sizeIterNum; sizeNum++)
	{
		for (int testNum = 0; testNum < _testIterNum; testNum++)
		{
			// string size
			int size = (sizeNum + 1) * _sizeOffset;

			// random
			boost::uint32_t seed;
			seed = static_cast<boost::uint32_t>(time(0));
			seed = seed + (testNum * sizeNum) % 1049861;
			boost::variate_generator<boost::lagged_fibonacci607, boost::uniform_int<>> rand(boost::lagged_fibonacci607(seed), boost::uniform_int<>(0, size));

			for (int stringNum = 0; stringNum < 2; stringNum++)
			{
				string tmp;
				int startPos = (rand() % 10007) * (rand() % 10007) % 20000003;
				_isData.seekg(startPos, ios::beg);

				char c;
				while (tmp.length() < size)
				{
					if (!_isData)
					{
						_isData.clear();
						_isData.seekg(0, ios::beg);
					}

					c = _isData.get();
					if ((int)c == -1)
						continue;

					tmp = tmp + c;
					if (stringNum == 0)
						alphabets[c].first++;
					else alphabets[c].second++;
				}

				if (stringNum == 0) strX = tmp;
				else strY = tmp;
				tmp.clear();
			}
			// If the lenght of any input string is not equal to test size, then reset the test.
			if (strX.length() != size || strY.length() != size)
			{
				testNum--;
				alphabets.clear();
				continue;
			}

			// 실험 번호, 문자열 사이즈 기록
			_osRecorder << sizeNum * _testIterNum + testNum << '\t' << size << '\t';

			returnPacket caseResult;

			// LCS
			caseResult = LCS(strX, strY);
			record(_osRecorder, caseResult);

			// Y. Sakai
			caseResult = MCS_0(strX, strY);
			record(_osRecorder, caseResult);

			// Lee
			caseResult = MCS_1(strX, strY);
			record(_osRecorder, caseResult);

			// k Linear Searching
			for (int k = 1; k <= 500; k++)
			{
				caseResult = MCS_T2(strX, strY, k);
				record(_osRecorder, caseResult);
			}

			// alphabets
			_osRecorder << "alphabets\t";
			int keyCnt = 0;
			for (map<char, pair<int, int>>::iterator iter = alphabets.begin(); iter != alphabets.end(); iter++)
			{
				_osRecorder << iter->first << '\t' << iter->second.first << '\t' << iter->second.second << '\t';
				keyCnt++;
			}
			_osRecorder << keyCnt << '\t';
			_osRecorder << '\n';

			strX.clear();
			strY.clear();
			alphabets.clear();

			cout << _dataName << sizeNum * _testIterNum + testNum << '\n';
		}
	}
}

void experiment_RandomData(int _alphabetSize, ofstream& _osRecorder, int _sizeIterNum, int _sizeOffset, int _testIterNum)
{
	// 문자열에 대한 문자 집합 분석
	// alphabets[char] = {strX에서의 발생횟수, strY에서의 발생횟수}
	map<wchar_t, pair<int, int>> alphabets;

	recordColumn(_osRecorder);

	for (int sizeNum = 0; sizeNum < _sizeIterNum; sizeNum++)
	{
		for (int testNum = 0; testNum < _testIterNum; testNum++)
		{
			// string size
			int size = (sizeNum + 1) * _sizeOffset;

			// random
			boost::uint32_t seed;
			seed = static_cast<boost::uint32_t>(time(0));
			seed = seed + (testNum * sizeNum) % 1049861;
			boost::variate_generator<boost::lagged_fibonacci607, boost::uniform_int<>> rand(boost::lagged_fibonacci607(seed), boost::uniform_int<>(1, _alphabetSize));

			for (int i = 0; i < size; i++)
			{
				wchar_t wc = (wchar_t)rand();
				wstrX.push_back(wc);

				alphabets[wc].first++;

				wc = (wchar_t)rand();
				wstrY.push_back(wc);

				alphabets[wc].second++;
			}

			if (wstrX.length() != size || wstrY.length() != size)
			{
				testNum--;
				alphabets.clear();
				continue;
			}

			// 실험 번호, 문자열 사이즈 기록
			_osRecorder << sizeNum * _testIterNum + testNum << '\t' << size << '\t';

			wreturnPacket caseResult;

			// LCS
			caseResult = LCS(wstrX, wstrY);
			record(_osRecorder, caseResult);

			// Y. Sakai
			caseResult = MCS_0(wstrX, wstrY);
			record(_osRecorder, caseResult);

			// Lee
			caseResult = MCS_1(wstrX, wstrY);
			record(_osRecorder, caseResult);

			// k Linear Searching
			for (int k = 1; k <= 500; k++)
			{
				caseResult = MCS_T2(wstrX, wstrY, k);
				record(_osRecorder, caseResult);
			}

			// alphabets
			_osRecorder << "alphabets\t";
			int keyCnt = 0;
			for (map<wchar_t, pair<int, int>>::iterator iter = alphabets.begin(); iter != alphabets.end(); iter++)
			{
				_osRecorder << iter->first << '\t' << iter->second.first << '\t' << iter->second.second << '\t';
				keyCnt++;
			}
			_osRecorder << keyCnt << '\t';

			_osRecorder << '\n';

			wstrX.clear();
			wstrY.clear();

			cout << _alphabetSize << ' ' << sizeNum * _testIterNum + testNum << '\n';
		}
	}
}
