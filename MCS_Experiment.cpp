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

// 저널에 추가됨
struct leeGhsA
{
	int g = -1;
	int h = -1;
	int s = -1;
	bool toggle = true;
};
returnPacket MCS_1_A(string _str1, string _str2);	

// k match Linear Max Search MCS (using vector)
// The Maximal Common Subsequence algorithms of Hyeonjun Shin and Jeong Seop Sim.
// 신현준, & 심정섭. (2022). 두 문자열의 극대공통부분서열을 찾는 새로운 알고리즘. 2022년 한국소프트웨어종합학술대회 논문집, 1212-1214.
// Hyeonjun Shin, & Jeon Seop Sim. A New Algorithm of Finding a Maximal Common Subsequence of Two Strings. Korea Software Congress 2022, 1212-1214.
struct kcGHS
{
	int id = -1;
	int g = -1;
	int h = -1;
	int s = 0;
};
returnPacket MCS_T1(string _str1, string _str2, int _k);
returnPacket MCS_T1_1(string _str1, string _str2, int _k);	// stack으로 수정 - 전혀 빠르지 않음: 보통 50만 번에 5만번 꼴로 빠름: 실험에서 제외

// new k Max Search MCS
// k 범위 안에서 가장 긴 것을 선택해서 MCS 진행
// 일단 갈기는 거임.
struct krGHS
{
	int id = -1;
	int g = 0;
	int h = 0;
	int gs = 0;
	int hs = 0;
	bool isKrunned = false;
};
returnPacket MCS_T2(string _str1, string _str2, int _k);


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
wreturnPacket MCS_1_A(wstring _str1, wstring _str2);
wreturnPacket MCS_T1(wstring _str1, wstring _str2, int _k);
wreturnPacket MCS_T1_1(wstring _str1, wstring _str2, int _k);
wreturnPacket MCS_T2(wstring _str1, wstring _str2, int _k);
bool isMCS(wstring _str1, wstring _str2, wstring _mcs);

// experiment
int funcNOs[4] = { 0, 1, 2, 3 };
int randomAlphabetSizes[11] = { 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048 };
string strX, strY; wstring wstrX, wstrY;
bool ostreamCheck(ofstream& _os, string _name);
bool istreamCheck(ifstream& _is, string _name);
void recordColumn(ofstream& _osLen, ofstream& _osT);
void record(ofstream& _osLen, ofstream& _osT, returnPacket& _rp);
void record(ofstream& _osLen, ofstream& _osT, wreturnPacket& _rp);
returnPacket experiment_function(int _funcNo, int _k, string _str1, string _str2);
void experiment_RealData(ifstream& _isData, ofstream& _osRecorder, ofstream& _osTimeRecorder, string _dataName, int _sizeIterNum, int _sizeOffset, int _testIterNum);
void experiment_RandomData(int _alphabetSize, ofstream& _osRecorder, ofstream& _osTimeRecorder, int _sizeIterNum, int _sizeOffset, int _testIterNum);
void experiment_RandomData_Norm(int _alphabetSize, ofstream& _osRecorder, ofstream& _osTimeRecorder, int _sizeIterNum, int _sizeOffset, int _testIterNum);

int main()
{
	// DNA
	string filename = "dna.txt";
	ifstream isDNA(filename, ios::in);
	if (!istreamCheck(isDNA, filename)) return -1;
	filename = "dnaResult.txt";
	ofstream osDNA(filename, ios::out);
	if (!ostreamCheck(osDNA, filename)) return -1;
	filename = "time" + filename;
	ofstream ostDNA(filename, ios::out);
	if (!ostreamCheck(ostDNA, filename)) return -1;

	experiment_RealData(isDNA, osDNA, ostDNA, "DNA", 10, 1000, 100);

	isDNA.clear(); isDNA.close();
	osDNA.clear(); osDNA.close();
	ostDNA.clear(); ostDNA.close();

	// Protein
	filename = "proteins.txt";
	ifstream isProtein(filename, ios::in);
	if(!istreamCheck(isProtein, filename)) return -1;
	filename = "proteinResult.txt";
	ofstream osProtein(filename, ios::out);
	if (!ostreamCheck(osProtein, filename)) return -1;
	filename = "time" + filename;
	ofstream ostProtein(filename, ios::out);
	if (!ostreamCheck(ostProtein, filename)) return -1;

	experiment_RealData(isProtein, osProtein, ostProtein, "Protein", 10, 1000, 100);

	isProtein.clear(); isProtein.close();
	osProtein.clear(); osProtein.close();
	ostProtein.clear(); ostProtein.close();

	// English
	filename = "english.txt";
	ifstream isEng(filename, ios::in);
	if(!istreamCheck(isEng, filename)) return -1;
	filename = "englishResult.txt";
	ofstream osEng(filename, ios::out);
	if(!ostreamCheck(osEng, filename)) return -1;
	filename = "time" + filename;
	ofstream ostEng(filename, ios::out);
	if (!ostreamCheck(ostEng, filename)) return -1;

	experiment_RealData(isEng, osEng, ostEng, "English", 10, 1000, 100);

	isEng.clear(); isEng.close();
	osEng.clear(); osEng.close();
	ostEng.clear(); ostEng.close();

	// random
	string prefix = "_alphabetResult.txt";
	for (int alphaIdx = 0; alphaIdx < 11; alphaIdx++)
	{
		filename = to_string(randomAlphabetSizes[alphaIdx]) + prefix;
		ofstream osRandom(filename, ios::out);
		if(!ostreamCheck(osRandom, filename)) return -1;
		filename = "time" + filename;
		ofstream ostRan(filename, ios::out);
		if (!ostreamCheck(ostRan, filename)) return -1;

		experiment_RandomData(randomAlphabetSizes[alphaIdx], osRandom, ostRan, 10, 1000, 100);

		osRandom.clear(); osRandom.close();
		ostRan.clear(); ostRan.close();
	}
	prefix = "_NormAlphabetResult.txt";
	for (int alphaIdx = 0; alphaIdx < 11; alphaIdx++)
	{
		filename = to_string(randomAlphabetSizes[alphaIdx]) + prefix;
		ofstream osRandom(filename, ios::out);
		if (!ostreamCheck(osRandom, filename)) return -1;
		filename = "time" + filename;
		ofstream ostRan(filename, ios::out);
		if (!ostreamCheck(ostRan, filename)) return -1;

		experiment_RandomData_Norm(randomAlphabetSizes[alphaIdx], osRandom, ostRan, 10, 1000, 100);

		osRandom.clear(); osRandom.close();
		ostRan.clear(); ostRan.close();
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


returnPacket MCS_1_A(string _str1, string _str2)
{
	// 시간 측정
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();

	string mcs = "";

	stack<leeGhsA> stk;
	leeGhsA init{ -1, -1, 0, true };
	stk.push(init);

	// 범위의 오른쪽 끝
	int iEnd = len1; int jEnd = len2;

	while (!stk.empty())
	{
		leeGhsA curGHS = stk.top();
		stk.pop();

		if (curGHS.toggle)
		{
			// s가 짝수이면 str1(X)에서 찾는 중이다.
			if (curGHS.s % 2 == 0)
			{
				int searching_g = curGHS.g + curGHS.s / 2 + 1;

				// 아직 str1(X)에서의 탐색 범위가 지정된 범위 내에 존재한다.
				if (searching_g < iEnd)
				{
					leeGhsA next_GHS{ curGHS.g, curGHS.h, curGHS.s + 1, curGHS.toggle };
					stk.push(next_GHS);

					int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);

					// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
					// 새로운 공통문자 offset을 저장
					if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd - 1)
					{
						// toggle change
						leeGhsA next_GHS{ searching_g, finding_h, 0, !curGHS.toggle };
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
					leeGhsA next_GHS{ curGHS.g, curGHS.h, curGHS.s + 1, curGHS.toggle };
					stk.push(next_GHS);

					int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

					// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
					// 새로운 공통문자 offset을 저장
					if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd - 1)
					{
						leeGhsA next_GHS{ finding_g, searching_h, 0, !curGHS.toggle };
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
		else
		{
			// s가 짝수이면 str1(Y)에서 찾는 중이다.
			if (curGHS.s % 2 == 0)
			{
				int searching_h = curGHS.h + curGHS.s / 2 + 1;

				// 아직 str2(Y)에서의 탐색 범위가 지정된 범위 내에 존재한다.
				if (searching_h < jEnd)
				{
					leeGhsA next_GHS{ curGHS.g, curGHS.h, curGHS.s + 1, curGHS.toggle };
					stk.push(next_GHS);

					int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

					// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
					// 새로운 공통문자 offset을 저장
					if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd - 1)
					{
						leeGhsA next_GHS{ finding_g, searching_h, 0, !curGHS.toggle };
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
			// s가 홀수이면 str2(X)에서 찾는 중이다.
			else
			{
				int searching_g = curGHS.g + (curGHS.s + 1) / 2;

				// 아직 str1(X)에서의 탐색 범위가 지정된 범위 내에 존재한다.
				if (searching_g < iEnd)
				{
					leeGhsA next_GHS{ curGHS.g, curGHS.h, curGHS.s + 1, curGHS.toggle };
					stk.push(next_GHS);

					int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);

					// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
					// 새로운 공통문자 offset을 저장
					if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd - 1)
					{
						leeGhsA next_GHS{ searching_g, finding_h, 0, !curGHS.toggle };
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


// k match Linear Max Search MCS (using vector)
// The Maximal Common Subsequence algorithms of Hyeonjun Shin and Jeong Seop Sim.
// 신현준, & 심정섭. (2022). 두 문자열의 극대공통부분서열을 찾는 새로운 알고리즘. 2022년 한국소프트웨어종합학술대회 논문집, 1212-1214.
// Hyeonjun Shin, & Jeon Seop Sim. A New Algorithm of Finding a Maximal Common Subsequence of Two Strings. Korea Software Congress 2022, 1212-1214.
int iEnd5, jEnd5;
bool compGHS(const kcGHS _a, const kcGHS _b)
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
returnPacket MCS_T1(string _str1, string _str2, int _k)
{
	// 시간 측정
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();

	string mcs = "";

	stack<kcGHS> stk;
	kcGHS init{ 0, -1, -1, 0 };
	stk.push(init);

	vector<vector<kcGHS>> wthPqs;
	vector<kcGHS> init_wthpq;
	wthPqs.push_back(init_wthpq);

	// 범위의 오른쪽 끝
	iEnd5 = len1, jEnd5 = len2;

	// check 용
	int GHS_id = 1;

	while (!stk.empty())
	{
		kcGHS curGHS = stk.top();

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
					kcGHS wthGHS = wthPqs[curGHS.id][0];
					int maxIdx = 0;
					for (int i = 1; i < wthPqs[curGHS.id].size(); i++)
					{
						kcGHS cur = wthPqs[curGHS.id][i];
						if (compGHS(cur, wthGHS))	// cur가 더 길게 만들면 true 반환
						{
							maxIdx = i;
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

						wthPqs[curGHS.id].erase(wthPqs[curGHS.id].begin() + maxIdx);
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
					kcGHS wthGHS = wthPqs[curGHS.id][0];
					int maxIdx = 0;
					for (int i = 1; i < wthPqs[curGHS.id].size(); i++)
					{
						kcGHS cur = wthPqs[curGHS.id][i];
						if (compGHS(cur, wthGHS))	// cur가 더 길게 만들면 true 반환
						{
							maxIdx = i;
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

						wthPqs[curGHS.id].erase(wthPqs[curGHS.id].begin() + maxIdx);
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
				kcGHS next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);

				// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
				// 새로운 공통문자 offset을 저장
				if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd5 - 1)
				{
					// k 범위 안에서 찾았다.
					if (searching_g <= curGHS.g + _k && finding_h <= curGHS.h + _k)
					{
						kcGHS next_GHS{ GHS_id++, searching_g, finding_h, 0 };
						stk.push(next_GHS);
					}
					// k 범위를 벗어나서 보류한다.
					else
					{
						while (wthPqs.size() <= curGHS.id)
						{
							vector<kcGHS> newQue;
							wthPqs.push_back(newQue);
						}
						kcGHS wth_GHS{ curGHS.id, searching_g, finding_h, curGHS.s };
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
				kcGHS next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

				// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
				// 새로운 공통문자 offset을 저장
				if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd5 - 1)
				{
					// k 범위 안에서 찾았다.
					if (searching_h <= curGHS.h + _k && finding_g <= curGHS.g + _k)
					{
						kcGHS next_GHS{ GHS_id++, finding_g, searching_h, 0 };
						stk.push(next_GHS);
					}
					// k 범위를 벗어나서 보류한다.
					else
					{
						while (wthPqs.size() <= curGHS.id)
						{
							vector<kcGHS> newQue;
							wthPqs.push_back(newQue);
						}
						kcGHS wth_GHS{ curGHS.id, finding_g, searching_h, curGHS.s };
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


bool compGHS_1(const ghs _a, const ghs _b)
{
	int min_by_a = iEnd5 - _a.g; int max_by_a = jEnd5 - _a.h;
	int min_by_b = iEnd5 - _b.g; int max_by_b = jEnd5 - _b.h;

	if (min_by_a > max_by_a) swap(min_by_a, max_by_a);
	if (min_by_b > max_by_b) swap(min_by_b, max_by_b);

	if (min_by_a != min_by_b)
		return min_by_a > min_by_b;
	else
		return max_by_a >= max_by_b;
}
returnPacket MCS_T1_1(string _str1, string _str2, int _k)
{
	// 시간 측정
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();

	string mcs = "";

	stack<ghs> stk;
	ghs init{-1, -1, 0 };
	stk.push(init);

	stack<vector<ghs>> wthStk;
	vector<ghs> init_wthVec;
	wthStk.push(init_wthVec);

	// 범위의 오른쪽 끝
	iEnd5 = len1, jEnd5 = len2;

	while (!stk.empty())
	{
		ghs curGHS = stk.top();
		vector<ghs> curWthVec = wthStk.top();

		bool isWthPoped = false;

		// (1) _k만큼 탐색을 했으나 _k 보다 작은 범위 내에서 공통 문자를 찾지 못한 경우
		// (2) X에서 _k만큼 탐색을 진행하기 전에 탐색범위가 종료된 경우
		// (3) Y에서 _k만큼 탐색을 진행하기 전에 탐색범위가 종료된 경우
		// (2)와 (3)의 경우 탐색이 종료될 예정이다. 그러나 탐색이 종료되었을 때, iEnd와 jEnd를 갱신하더라도
		// new_iEnd(new_jEnd)와 old_iEnd(old_jEnd) 사이에 공통문자가 보류된 상태로 남을 수 있다.
		// 즉, 보류 공통문자가 들어가지 않을 수 있으므로 (1)의 경우와 같이 처리해야 한다.
		if (curGHS.s >= _k * 2 || (curGHS.s % 2 == 0 && curGHS.g + curGHS.s / 2 + 1 >= iEnd5) ||
			(curGHS.s % 2 == 1 && curGHS.h + (curGHS.s + 1) / 2 >= jEnd5))
		{
			// 보류한 것이 없으면 일반 탐색으로 진행
			if (!curWthVec.empty())
			{
				ghs maxWth = curWthVec[0];
				int maxIdx = 0;
				for (int i = 1; i < curWthVec.size(); i++)
				{
					ghs cur = curWthVec[i];
					if (compGHS_1(cur, maxWth))
					{
						maxWth = cur;
						maxIdx = i;
					}
				}

				if (maxWth.g < iEnd5 && maxWth.h < jEnd5)
				{
					isWthPoped = true;

					maxWth.s = 0;
					stk.push(maxWth);

					curWthVec.erase(curWthVec.begin() + maxIdx);
					wthStk.pop();
					wthStk.push(curWthVec);
					vector<ghs> new_wthVec;
					wthStk.push(new_wthVec);
				}
				else
				{
					curWthVec.erase(curWthVec.begin() + maxIdx);
					wthStk.pop();
					wthStk.push(curWthVec);
				}
			}
		}

		if (isWthPoped)
		{
			curGHS = stk.top();
			curWthVec = wthStk.top();
		}

		stk.pop();
		wthStk.pop();

		// s가 짝수이면 str1(X)에서 찾는 중이다.
		if (curGHS.s % 2 == 0)
		{
			int searching_g = curGHS.g + curGHS.s / 2 + 1;

			// 아직 str1(X)에서의 탐색 범위가 지정된 범위 내에 존재한다.
			if (searching_g < iEnd5)
			{
				ghs next_GHS{ curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);

				// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
				// 새로운 공통문자 offset을 저장
				if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd5 - 1)
				{
					// k 범위 안에서 찾았다.
					if (searching_g <= curGHS.g + _k && finding_h <= curGHS.h + _k)
					{
						ghs next_GHS{ searching_g, finding_h, 0 };
						stk.push(next_GHS);
						wthStk.push(curWthVec);
						vector<ghs> next_wthVec;
						wthStk.push(next_wthVec);
					}
					// k 범위를 벗어나서 보류한다.
					else
					{
						ghs wth_GHS{ searching_g, finding_h, curGHS.s };
						curWthVec.push_back(wth_GHS);
						wthStk.push(curWthVec);
					}
				}
				else
					wthStk.push(curWthVec);
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
				ghs next_GHS{ curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

				// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
				// 새로운 공통문자 offset을 저장
				if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd5 - 1)
				{
					// k 범위 안에서 찾았다.
					if (searching_h <= curGHS.h + _k && finding_g <= curGHS.g + _k)
					{
						ghs next_GHS{ finding_g, searching_h, 0 };
						stk.push(next_GHS);
						wthStk.push(curWthVec);
						vector<ghs> next_wthVec;
						wthStk.push(next_wthVec);
					}
					// k 범위를 벗어나서 보류한다.
					else
					{
						ghs wth_GHS{ finding_g, searching_h, curGHS.s };
						curWthVec.push_back(wth_GHS);
						wthStk.push(curWthVec);
					}
				}
				else
					wthStk.push(curWthVec);
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


// 실패중
int iEnd6; int jEnd6;
bool compPOS(const pair<int,int> _a, const pair<int,int> _b)
{
	int min_by_a = iEnd6 - _a.first; int max_by_a = jEnd6 - _a.second;
	int min_by_b = iEnd6 - _b.first; int max_by_b = jEnd6 - _b.second;

	if (min_by_a > max_by_a) swap(min_by_a, max_by_a);
	if (min_by_b > max_by_b) swap(min_by_b, max_by_b);

	if (min_by_a != min_by_b)
		return min_by_a > min_by_b;
	else
		return max_by_a >= max_by_b;
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

	stack<krGHS> stk;
	krGHS init{ 0, -1, -1, 1, 1, false };
	stk.push(init);

	vector<vector<pair<int, int>>> wthVec(len1 * 2 + 4);

	// 범위의 오른쪽 끝
	iEnd6 = len1, jEnd6 = len2;
	
	int ghsId = 1;
	while (!stk.empty())
	{
		krGHS curGHS = stk.top();

		// X와 Y에서 탐색할 문자의 위치를 구한다.
		int searching_g = curGHS.g + curGHS.gs;
		int searching_h = curGHS.h + curGHS.hs;

		if (!wthVec[curGHS.id].empty())
		{
			pair<int, int> bestWth = wthVec[curGHS.id][0];
			int bestIdx = 0;
			for (int i = 1; i < wthVec[curGHS.id].size(); i++)
			{
				if (compPOS(wthVec[curGHS.id][i], bestWth))
				{
					bestWth = wthVec[curGHS.id][i];
					bestIdx = i;
				}
			}

			wthVec[curGHS.id].erase(wthVec[curGHS.id].begin() + bestIdx);

			if (bestWth.first < iEnd6 && bestWth.second < jEnd6)
			{
				krGHS new_ghs{ ghsId++, bestWth.first, bestWth.second, 1, 1, false };
				stk.push(new_ghs);
			}
			continue;
		}

		stk.pop();

		// 두 곳 모두 탐색 범위를 넘지 않는 경우
		if (searching_g < iEnd6 && searching_h < jEnd6)
		{
			int finding_h, finding_g;

			if (!curGHS.isKrunned)
			{
				for (int i = 0; i < _k; i++)
				{
					// 탐색한 문자가 밖에 있으면 find_first_of => -1로 걸러짐
					// 탐색에 이용할 문자에 대한 인덱스가 넘을 수도 있음.
					if (searching_g + i >= iEnd6 || searching_h + i >= jEnd6)
						break;

					finding_h = _str2.find_first_of(_str1[searching_g + i], curGHS.h + 1);
					finding_g = _str1.find_first_of(_str2[searching_h + i], curGHS.g + 1);

					if (finding_h != _str2.npos && finding_h != -1 && finding_h < jEnd6)
						wthVec[curGHS.id].push_back({ searching_g + i, finding_h });
					if (finding_g != _str1.npos && finding_g != -1 && finding_g < iEnd6)
						wthVec[curGHS.id].push_back({ finding_g, searching_h + i});
				}
				curGHS.isKrunned = true;
				curGHS.gs = curGHS.hs = _k + 1;
				stk.push(curGHS);
				continue;
			}

			//더 이상 가져다 쓸게 없음이 보장
			finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);
			finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

			// 양쪽 기준의 문자가 다른 문자열에서 나타난 경우
			if ((finding_h != _str2.npos && finding_h != -1 && finding_h < jEnd6) &&
				(finding_g != _str1.npos && finding_g != -1 && finding_g < iEnd6))
			{
				int minLen_by_g = min((iEnd6 - searching_g), (jEnd6 - finding_h));
				int minLen_by_h = min((iEnd6 - finding_g), (jEnd6 - searching_h));

				// X에서 찾은 문자를 기준으로 할 때 범위가 더 넓어지는 경우
				if (minLen_by_g >= minLen_by_h)
				{
					krGHS reGHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs, curGHS.isKrunned };
					stk.push(reGHS);

					krGHS next_GHS{ ghsId++, searching_g, finding_h, 1, 1 };
					stk.push(next_GHS);
				}
				// Y에서 찾은 문자를 기준으로 할 때 범위가 더 넓어지는 경우
				else
				{
					krGHS reGHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.gs, curGHS.hs + 1, curGHS.isKrunned };
					stk.push(reGHS);

					krGHS next_GHS{ghsId++, finding_g, searching_h, 1, 1 };
					stk.push(next_GHS);
				}
			}
			// X 기준의 문자가 str2에서 나타난 경우만
			else if (finding_h != _str2.npos && finding_h != -1 && finding_h < jEnd6)
			{
				krGHS reGHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs + 1, curGHS.isKrunned };
				stk.push(reGHS);

				krGHS next_GHS{ghsId++, searching_g, finding_h, 1, 1 };
				stk.push(next_GHS);
			}
			// Y 기준의 문자가 str1에서 나타난 경우만
			else if (finding_g != _str1.npos && finding_g != -1 && finding_g < iEnd6)
			{
				krGHS reGHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs + 1, curGHS.isKrunned };
				stk.push(reGHS);

				krGHS next_GHS{ghsId++, finding_g, searching_h, 1, 1 };
				stk.push(next_GHS);
			}
			// 둘다 안 나타난 경우
			else
			{
				krGHS reGHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs + 1, curGHS.isKrunned };
				stk.push(reGHS);
			}
		}
		else if (curGHS.g >= 0)
		{
			mcs = _str1[curGHS.g] + mcs;

			iEnd6 = _str1.find_last_of(_str1[curGHS.g], iEnd6 - 1);
			jEnd6 = _str2.find_last_of(_str1[curGHS.g], jEnd6 - 1);
			if (iEnd6 == 0 || jEnd6 == 0 ||
				iEnd6 == -1 || jEnd6 == -1 ||
				iEnd6 == _str1.npos || jEnd6 == _str2.npos)
				break;
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
			//cout << "isMCS ERROR: There is no character\n";
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
					//cout << "isMCS ERROR: Not disjoint" << endl;
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
					//cout << "isMCS ERROR: Not disjoint" << endl;
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

wreturnPacket MCS_1_A(wstring _str1, wstring _str2)
{
	// 시간 측정
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();

	wstring mcs = L"";

	stack<leeGhsA> stk;
	leeGhsA init{ -1, -1, 0, true };
	stk.push(init);

	// 범위의 오른쪽 끝
	int iEnd = len1; int jEnd = len2;

	while (!stk.empty())
	{
		leeGhsA curGHS = stk.top();
		stk.pop();

		if (curGHS.toggle)
		{
			// s가 짝수이면 str1(X)에서 찾는 중이다.
			if (curGHS.s % 2 == 0)
			{
				int searching_g = curGHS.g + curGHS.s / 2 + 1;

				// 아직 str1(X)에서의 탐색 범위가 지정된 범위 내에 존재한다.
				if (searching_g < iEnd)
				{
					leeGhsA next_GHS{ curGHS.g, curGHS.h, curGHS.s + 1, curGHS.toggle };
					stk.push(next_GHS);

					int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);

					// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
					// 새로운 공통문자 offset을 저장
					if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd - 1)
					{
						// toggle change
						leeGhsA next_GHS{ searching_g, finding_h, 0, !curGHS.toggle };
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
					leeGhsA next_GHS{ curGHS.g, curGHS.h, curGHS.s + 1, curGHS.toggle };
					stk.push(next_GHS);

					int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

					// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
					// 새로운 공통문자 offset을 저장
					if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd - 1)
					{
						leeGhsA next_GHS{ finding_g, searching_h, 0, !curGHS.toggle };
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
		else
		{
			// s가 짝수이면 str1(Y)에서 찾는 중이다.
			if (curGHS.s % 2 == 0)
			{
				int searching_h = curGHS.h + curGHS.s / 2 + 1;

				// 아직 str2(Y)에서의 탐색 범위가 지정된 범위 내에 존재한다.
				if (searching_h < jEnd)
				{
					leeGhsA next_GHS{ curGHS.g, curGHS.h, curGHS.s + 1, curGHS.toggle };
					stk.push(next_GHS);

					int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

					// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
					// 새로운 공통문자 offset을 저장
					if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd - 1)
					{
						leeGhsA next_GHS{ finding_g, searching_h, 0, !curGHS.toggle };
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
			// s가 홀수이면 str2(X)에서 찾는 중이다.
			else
			{
				int searching_g = curGHS.g + (curGHS.s + 1) / 2;

				// 아직 str1(X)에서의 탐색 범위가 지정된 범위 내에 존재한다.
				if (searching_g < iEnd)
				{
					leeGhsA next_GHS{ curGHS.g, curGHS.h, curGHS.s + 1, curGHS.toggle };
					stk.push(next_GHS);

					int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);

					// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
					// 새로운 공통문자 offset을 저장
					if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd - 1)
					{
						leeGhsA next_GHS{ searching_g, finding_h, 0, !curGHS.toggle };
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

wreturnPacket MCS_T1(wstring _str1, wstring _str2, int _k)
{
	// 시간 측정
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();

	wstring mcs = L"";

	stack<kcGHS> stk;
	kcGHS init{ 0, -1, -1, 0 };
	stk.push(init);

	vector<vector<kcGHS>> wthPqs;
	vector<kcGHS> init_wthpq;
	wthPqs.push_back(init_wthpq);

	// 범위의 오른쪽 끝
	iEnd5 = len1, jEnd5 = len2;

	// check 용
	int GHS_id = 1;

	while (!stk.empty())
	{
		kcGHS curGHS = stk.top();

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
					kcGHS wthGHS = wthPqs[curGHS.id][0];
					for (int i = 1; i < wthPqs[curGHS.id].size(); i++)
					{
						kcGHS cur = wthPqs[curGHS.id][i];
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
					kcGHS wthGHS = wthPqs[curGHS.id][0];
					for (int i = 1; i < wthPqs[curGHS.id].size(); i++)
					{
						kcGHS cur = wthPqs[curGHS.id][i];
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
				kcGHS next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);

				// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
				// 새로운 공통문자 offset을 저장
				if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd5 - 1)
				{
					// k 범위 안에서 찾았다.
					if (searching_g <= curGHS.g + _k && finding_h <= curGHS.h + _k)
					{
						kcGHS next_GHS{ GHS_id++, searching_g, finding_h, 0 };
						stk.push(next_GHS);
					}
					// k 범위를 벗어나서 보류한다.
					else
					{
						while (wthPqs.size() <= curGHS.id)
						{
							vector<kcGHS> newQue;
							wthPqs.push_back(newQue);
						}
						kcGHS wth_GHS{ curGHS.id, searching_g, finding_h, curGHS.s };
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
				kcGHS next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

				// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
				// 새로운 공통문자 offset을 저장
				if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd5 - 1)
				{
					// k 범위 안에서 찾았다.
					if (searching_h <= curGHS.h + _k && finding_g <= curGHS.g + _k)
					{
						kcGHS next_GHS{ GHS_id++, finding_g, searching_h, 0 };
						stk.push(next_GHS);
					}
					// k 범위를 벗어나서 보류한다.
					else
					{
						while (wthPqs.size() <= curGHS.id)
						{
							vector<kcGHS> newQue;
							wthPqs.push_back(newQue);
						}
						kcGHS wth_GHS{ curGHS.id, finding_g, searching_h, curGHS.s };
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

wreturnPacket MCS_T1_1(wstring _str1, wstring _str2, int _k)
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

	stack<vector<ghs>> wthStk;
	vector<ghs> init_wthVec;
	wthStk.push(init_wthVec);

	// 범위의 오른쪽 끝
	iEnd5 = len1, jEnd5 = len2;

	// check 용
	int GHS_id = 1;

	while (!stk.empty())
	{
		ghs curGHS = stk.top();
		vector<ghs> curWthVec = wthStk.top();

		bool isWthPoped = false;

		// (1) _k만큼 탐색을 했으나 _k 보다 작은 범위 내에서 공통 문자를 찾지 못한 경우
		// (2) X에서 _k만큼 탐색을 진행하기 전에 탐색범위가 종료된 경우
		// (3) Y에서 _k만큼 탐색을 진행하기 전에 탐색범위가 종료된 경우
		// (2)와 (3)의 경우 탐색이 종료될 예정이다. 그러나 탐색이 종료되었을 때, iEnd와 jEnd를 갱신하더라도
		// new_iEnd(new_jEnd)와 old_iEnd(old_jEnd) 사이에 공통문자가 보류된 상태로 남을 수 있다.
		// 즉, 보류 공통문자가 들어가지 않을 수 있으므로 (1)의 경우와 같이 처리해야 한다.
		if (curGHS.s >= _k * 2 || (curGHS.s % 2 == 0 && curGHS.g + curGHS.s / 2 + 1 >= iEnd5) ||
			(curGHS.s % 2 == 1 && curGHS.h + (curGHS.s + 1) / 2 >= jEnd5))
		{
			// 보류한 것이 없으면 일반 탐색으로 진행
			if (!curWthVec.empty())
			{
				ghs maxWth = curWthVec[0];
				int maxIdx = 0;
				for (int i = 1; i < curWthVec.size(); i++)
				{
					ghs cur = curWthVec[i];
					if (compGHS_1(cur, maxWth))
					{
						maxWth = cur;
						maxIdx = i;
					}
				}

				if (maxWth.g < iEnd5 && maxWth.h < jEnd5)
				{
					isWthPoped = true;

					maxWth.s = 0;
					stk.push(maxWth);

					curWthVec.erase(curWthVec.begin() + maxIdx);
					wthStk.pop();
					wthStk.push(curWthVec);
					vector<ghs> new_wthVec;
					wthStk.push(new_wthVec);
				}
				else
				{
					curWthVec.erase(curWthVec.begin() + maxIdx);
					wthStk.pop();
					wthStk.push(curWthVec);
				}

			}
		}

		if (isWthPoped)
		{
			curGHS = stk.top();
			curWthVec = wthStk.top();
		}

		stk.pop();
		wthStk.pop();

		// s가 짝수이면 str1(X)에서 찾는 중이다.
		if (curGHS.s % 2 == 0)
		{
			int searching_g = curGHS.g + curGHS.s / 2 + 1;

			// 아직 str1(X)에서의 탐색 범위가 지정된 범위 내에 존재한다.
			if (searching_g < iEnd5)
			{
				ghs next_GHS{ curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);

				// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
				// 새로운 공통문자 offset을 저장
				if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd5 - 1)
				{
					// k 범위 안에서 찾았다.
					if (searching_g <= curGHS.g + _k && finding_h <= curGHS.h + _k)
					{
						ghs next_GHS{ searching_g, finding_h, 0 };
						stk.push(next_GHS);
						wthStk.push(curWthVec);
						vector<ghs> next_wthVec;
						wthStk.push(next_wthVec);
					}
					// k 범위를 벗어나서 보류한다.
					else
					{
						ghs wth_GHS{ searching_g, finding_h, curGHS.s };
						curWthVec.push_back(wth_GHS);
						wthStk.push(curWthVec);
					}
				}
				else
					wthStk.push(curWthVec);
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
				ghs next_GHS{ curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

				// 지정된 str2(Y)의 범위 내에서 공통 문자를 찾음
				// 새로운 공통문자 offset을 저장
				if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd5 - 1)
				{
					// k 범위 안에서 찾았다.
					if (searching_h <= curGHS.h + _k && finding_g <= curGHS.g + _k)
					{
						ghs next_GHS{ finding_g, searching_h, 0 };
						stk.push(next_GHS);
						wthStk.push(curWthVec);
						vector<ghs> next_wthVec;
						wthStk.push(next_wthVec);
					}
					// k 범위를 벗어나서 보류한다.
					else
					{
						ghs wth_GHS{ finding_g, searching_h, curGHS.s };
						curWthVec.push_back(wth_GHS);
						wthStk.push(curWthVec);
					}
				}
				else
					wthStk.push(curWthVec);
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

wreturnPacket MCS_T2(wstring _str1, wstring _str2, int _k)
{
	// 시간 측정
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();

	wstring mcs = L"";

	stack<krGHS> stk;
	krGHS init{ 0, -1, -1, 1, 1, false };
	stk.push(init);

	vector<vector<pair<int, int>>> wthVec(len1 * 2 + 4);

	// 범위의 오른쪽 끝
	iEnd6 = len1, jEnd6 = len2;

	int ghsId = 1;
	while (!stk.empty())
	{
		krGHS curGHS = stk.top();

		// X와 Y에서 탐색할 문자의 위치를 구한다.
		int searching_g = curGHS.g + curGHS.gs;
		int searching_h = curGHS.h + curGHS.hs;

		if (!wthVec[curGHS.id].empty())
		{
			pair<int, int> bestWth = wthVec[curGHS.id][0];
			int bestIdx = 0;
			for (int i = 1; i < wthVec[curGHS.id].size(); i++)
			{
				if (compPOS(wthVec[curGHS.id][i], bestWth))
				{
					bestWth = wthVec[curGHS.id][i];
					bestIdx = i;
				}
			}

			wthVec[curGHS.id].erase(wthVec[curGHS.id].begin() + bestIdx);

			if (bestWth.first < iEnd6 && bestWth.second < jEnd6)
			{
				krGHS new_ghs{ ghsId++, bestWth.first, bestWth.second, 1, 1, false };
				stk.push(new_ghs);
			}
			continue;
		}

		stk.pop();

		// 두 곳 모두 탐색 범위를 넘지 않는 경우
		if (searching_g < iEnd6 && searching_h < jEnd6)
		{
			int finding_h, finding_g;

			if (!curGHS.isKrunned)
			{
				for (int i = 0; i < _k; i++)
				{
					// 탐색한 문자가 밖에 있으면 find_first_of => -1로 걸러짐
					// 탐색에 이용할 문자에 대한 인덱스가 넘을 수도 있음.
					if (searching_g + i >= iEnd6 || searching_h + i >= jEnd6)
						break;

					finding_h = _str2.find_first_of(_str1[searching_g + i], curGHS.h + 1);
					finding_g = _str1.find_first_of(_str2[searching_h + i], curGHS.g + 1);

					if (finding_h != _str2.npos && finding_h != -1 && finding_h < jEnd6)
						wthVec[curGHS.id].push_back({ searching_g + i, finding_h });
					if (finding_g != _str1.npos && finding_g != -1 && finding_g < iEnd6)
						wthVec[curGHS.id].push_back({ finding_g, searching_h + i });
				}
				curGHS.isKrunned = true;
				curGHS.gs = curGHS.hs = _k + 1;
				stk.push(curGHS);
				continue;
			}

			//더 이상 가져다 쓸게 없음이 보장
			finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);
			finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

			// 양쪽 기준의 문자가 다른 문자열에서 나타난 경우
			if ((finding_h != _str2.npos && finding_h != -1 && finding_h < jEnd6) &&
				(finding_g != _str1.npos && finding_g != -1 && finding_g < iEnd6))
			{
				int minLen_by_g = min((iEnd6 - searching_g), (jEnd6 - finding_h));
				int minLen_by_h = min((iEnd6 - finding_g), (jEnd6 - searching_h));

				// X에서 찾은 문자를 기준으로 할 때 범위가 더 넓어지는 경우
				if (minLen_by_g >= minLen_by_h)
				{
					krGHS reGHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs, curGHS.isKrunned };
					stk.push(reGHS);

					krGHS next_GHS{ ghsId++, searching_g, finding_h, 1, 1 };
					stk.push(next_GHS);
				}
				// Y에서 찾은 문자를 기준으로 할 때 범위가 더 넓어지는 경우
				else
				{
					krGHS reGHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.gs, curGHS.hs + 1, curGHS.isKrunned };
					stk.push(reGHS);

					krGHS next_GHS{ ghsId++, finding_g, searching_h, 1, 1 };
					stk.push(next_GHS);
				}
			}
			// X 기준의 문자가 str2에서 나타난 경우만
			else if (finding_h != _str2.npos && finding_h != -1 && finding_h < jEnd6)
			{
				krGHS reGHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs + 1, curGHS.isKrunned };
				stk.push(reGHS);

				krGHS next_GHS{ ghsId++, searching_g, finding_h, 1, 1 };
				stk.push(next_GHS);
			}
			// Y 기준의 문자가 str1에서 나타난 경우만
			else if (finding_g != _str1.npos && finding_g != -1 && finding_g < iEnd6)
			{
				krGHS reGHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs + 1, curGHS.isKrunned };
				stk.push(reGHS);

				krGHS next_GHS{ ghsId++, finding_g, searching_h, 1, 1 };
				stk.push(next_GHS);
			}
			// 둘다 안 나타난 경우
			else
			{
				krGHS reGHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs + 1, curGHS.isKrunned };
				stk.push(reGHS);
			}
		}
		else if (curGHS.g >= 0)
		{
			mcs = _str1[curGHS.g] + mcs;

			iEnd6 = _str1.find_last_of(_str1[curGHS.g], iEnd6 - 1);
			jEnd6 = _str2.find_last_of(_str1[curGHS.g], jEnd6 - 1);
			if (iEnd6 == 0 || jEnd6 == 0 ||
				iEnd6 == -1 || jEnd6 == -1 ||
				iEnd6 == _str1.npos || jEnd6 == _str2.npos)
				break;
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
			//cout << "isMCS ERROR: There is no character";
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
				//	cout << "ERROR 3" << endl;
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
					//cout << "ERROR 3" << endl;
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
void recordColumn(ofstream& _osLen, ofstream& _osT)
{
	_osLen << "testNumber\t|X|\t"
		<< "LCS_Len\t" << "Sakai_Len\t"	<< "LeeP_Len\t" << "LeeA_Len\t";

	for (int k = 1; k <= 128; k++)
	{
		_osLen << "kc_" << k << "_Len\t";
	}
	for (int k = 1; k <= 128; k++)
	{
		_osLen << "kr_" << k << "_Len\t";
	}

	_osLen << "alphabet_Size\tVx\tVy\n";

	_osT << "testNumber\t|X|\t"
		<< "LCS_Len\t" << "Sakai_Len\t" << "LeeP_Len\t" << "LeeA_Len\t";

	for (int k = 1; k <= 128; k++)
	{
		_osT << k << "_Len\t";
	}
	for (int k = 1; k <= 128; k++)
	{
		_osT << "kr_" << k << "_Len\t";
	}

	_osT << "\n";
}
// older Recorder: record all data
//void record(ofstream& _os, returnPacket& _rp)
//{
//	if (!_rp.isMaximal)
//		exit(-1);
//	_os << _rp.cssLen << '\t' << _rp.isMaximal << '\t' << (double)((double)_rp.cssLen / (double)_rp.n * 100 * 2) << '\t' << _rp.time << '\t';
//}
//void record(ofstream& _os, wreturnPacket& _rp)
//{
//	if (!_rp.isMaximal)
//		exit(-1);
//	_os << _rp.cssLen << '\t' << _rp.isMaximal << '\t' << (double)((double)_rp.cssLen / (double)_rp.n * 100 * 2) << '\t' << _rp.time << '\t';
//}

void record(ofstream& _osLen, ofstream& _osT, returnPacket& _rp)
{
	_osLen << (double)((double)_rp.cssLen / (double)_rp.n * 100 * 2) << '\t';
	_osT << _rp.time << '\t';
}
void record(ofstream& _osLen, ofstream& _osT, wreturnPacket& _rp)
{
	_osLen << (double)((double)_rp.cssLen / (double)_rp.n * 100 * 2) << '\t';
	_osT << _rp.time << '\t';
}

void experiment_RealData(ifstream& _isData, ofstream& _osRecorder, ofstream& _osTimeRecorder, string _dataName, int _sizeIterNum, int _sizeOffset, int _testIterNum)
{
	// 문자열에 대한 문자 집합 분석
	// alphabets[char] = {strX에서의 발생횟수, strY에서의 발생횟수}
	map<char, pair<int, int>> alphabets;

	recordColumn(_osRecorder, _osTimeRecorder);

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
				_isData.seekg(startPos, ios::cur);

				char c;
				while (tmp.length() < size)
				{
					if (!_isData)
					{
						_isData.clear();
						_isData.seekg(0, ios::beg);
					}

					c = _isData.get();
					while(c == '\n' || c == '\t')
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
			_osTimeRecorder << sizeNum * _testIterNum + testNum << '\t' << size << '\t';

			returnPacket caseResult;

			// LCS
			caseResult = LCS(strX, strY);
			record(_osRecorder, _osTimeRecorder, caseResult);
			if (!caseResult.isMaximal)
			{
				cout << "LCS" << caseResult.n << "is not maximal\n";
				exit(1);
			}

			// Y. Sakai
			caseResult = MCS_0(strX, strY);
			record(_osRecorder, _osTimeRecorder, caseResult);
			if (!caseResult.isMaximal)
			{
				cout << "Sakai" << caseResult.n << "is not maximal\n";
				exit(1);
			}

			// LeeP
			caseResult = MCS_1(strX, strY);
			record(_osRecorder, _osTimeRecorder, caseResult);
			if (!caseResult.isMaximal)
			{
				cout << "LeeP" << caseResult.n << "is not maximal\n";
				exit(1);
			}

			// LeeA
			caseResult = MCS_1_A(strX, strY);
			record(_osRecorder, _osTimeRecorder, caseResult);
			if (!caseResult.isMaximal)
			{
				cout << "LeeA" << caseResult.n << "is not maximal\n";
				exit(1);
			}

			// k Linear Searching
			for (int k = 1; k <= 128; k++)
			{
				caseResult = MCS_T1(strX, strY, k);
				if (!caseResult.isMaximal)
				{
					cout << "kc" << caseResult.n << "is not maximal\n";
					exit(1);
				}
				record(_osRecorder, _osTimeRecorder, caseResult);
			}

			// Except kr because it is not good algorithm
			//// k Range
			//for (int k = 1; k <= 128; k++)
			//{
			//	caseResult = MCS_T2(strX, strY, k);
			//	if (!caseResult.isMaximal)
			//	{
			//		cout << "kr" << caseResult.n << "is not maximal\n";
			//		exit(1);
			//	}
			//	record(_osRecorder, _osTimeRecorder, caseResult);
			//}
			

			// alphabets
			_osRecorder << alphabets.size() << '\t';
			vector<long double> possibilities[2];
			long double ePs[2]{ 0 };
			for (map<char, pair<int, int>>::iterator iter = alphabets.begin(); iter != alphabets.end(); iter++)
			{
				long double XP = (long double)(iter->second.first / (long double)size);
				long double YP = (long double)(iter->second.first / (long double)size);
				possibilities[0].push_back(XP);
				possibilities[1].push_back(YP);
				ePs[0] += XP;
				ePs[1] += YP;
			}
			ePs[0] /= (long double)alphabets.size();
			ePs[1] /= (long double)alphabets.size();

			long double vPs[2]{ 0 };
			for (int i = 0; i < possibilities[0].size(); i++)
			{
				vPs[0] += pow(ePs[0] - possibilities[0][i], 2);
				vPs[1] += pow(ePs[1] - possibilities[1][i], 2);
			}
			vPs[0] /= (long double)alphabets.size();
			vPs[1] /= (long double)alphabets.size();
			_osRecorder << vPs[0] << '\t' << vPs[1];
			_osRecorder << '\n';
			_osTimeRecorder << '\n';

			strX.clear();
			strY.clear();
			alphabets.clear();

			cout << _dataName << sizeNum * _testIterNum + testNum << '\n';
		}
	}
}

void experiment_RandomData(int _alphabetSize, ofstream& _osRecorder, ofstream& _osTimeRecorder, int _sizeIterNum, int _sizeOffset, int _testIterNum)
{
	// 문자열에 대한 문자 집합 분석
	// alphabets[char] = {strX에서의 발생횟수, strY에서의 발생횟수}
	map<wchar_t, pair<int, int>> alphabets;

	recordColumn(_osRecorder, _osTimeRecorder);

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
			boost::variate_generator<boost::mt19937, boost::uniform_int<>> rand(boost::mt19937(seed), boost::uniform_int<>(1, _alphabetSize));

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
			_osTimeRecorder << sizeNum * _testIterNum + testNum << '\t' << size << '\t';

			wreturnPacket caseResult;

			// LCS
			caseResult = LCS(wstrX, wstrY);
			if (!caseResult.isMaximal)
			{
				cout << "LCS" << caseResult.n << "is not maximal\n";
				exit(1);
			}
			record(_osRecorder, _osTimeRecorder, caseResult);


			// Y. Sakai
			caseResult = MCS_0(wstrX, wstrY);
			if (!caseResult.isMaximal)
			{
				cout << "Sakai" << caseResult.n << "is not maximal\n";
				exit(1);
			}
			record(_osRecorder, _osTimeRecorder, caseResult);

			// LeeP
			caseResult = MCS_1(wstrX, wstrY);
			record(_osRecorder, _osTimeRecorder, caseResult);
			if (!caseResult.isMaximal)
			{
				cout << "LeeP" << caseResult.n << "is not maximal\n";
				exit(1);
			}

			// LeeA
			caseResult = MCS_1_A(wstrX, wstrY);
			record(_osRecorder, _osTimeRecorder, caseResult);
			if (!caseResult.isMaximal)
			{
				cout << "LeeA" << caseResult.n << "is not maximal\n";
				exit(1);
			}

			// k Linear Searching
			for (int k = 1; k <= 128; k++)
			{
				caseResult = MCS_T1(wstrX, wstrY, k);
				if (!caseResult.isMaximal)
				{
					cout << "kc" << caseResult.n << "is not maximal\n";
					exit(1);
				}
				record(_osRecorder, _osTimeRecorder, caseResult);
			}

			// Except kr because it is not good algorithm
			//// k Range
			//for (int k = 1; k <= 128; k++)
			//{
			//	caseResult = MCS_T2(wstrX, wstrY, k);
			//	if (!caseResult.isMaximal)
			//	{
			//		cout << "kr" << caseResult.n << "is not maximal\n";
			//		exit(1);
			//	}
			//	record(_osRecorder, _osTimeRecorder, caseResult);
			//}

			// alphabets
			_osRecorder << alphabets.size() << '\t';
			vector<long double> possibilities[2];
			long double ePs[2]{ 0 };
			for (map<wchar_t, pair<int, int>>::iterator iter = alphabets.begin(); iter != alphabets.end(); iter++)
			{
				long double XP = (long double)(iter->second.first / (long double)size);
				long double YP = (long double)(iter->second.first / (long double)size);
				possibilities[0].push_back(XP);
				possibilities[1].push_back(YP);
				ePs[0] += XP;
				ePs[1] += YP;
			}
			ePs[0] /= (long double)alphabets.size();
			ePs[1] /= (long double)alphabets.size();

			long double vPs[2]{ 0 };
			for (int i = 0; i < possibilities[0].size(); i++)
			{
				vPs[0] += pow(ePs[0] - possibilities[0][i], 2);
				vPs[1] += pow(ePs[1] - possibilities[1][i], 2);
			}
			vPs[0] /= (long double)alphabets.size();
			vPs[1] /= (long double)alphabets.size();
			_osRecorder << vPs[0] << '\t' << vPs[1];

			_osRecorder << '\n';
			_osTimeRecorder << '\n';

			wstrX.clear();
			wstrY.clear();

			cout << _alphabetSize << ' ' << sizeNum * _testIterNum + testNum << '\n';
		}
	}
}

void experiment_RandomData_Norm(int _alphabetSize, ofstream& _osRecorder, ofstream& _osTimeRecorder, int _sizeIterNum, int _sizeOffset, int _testIterNum)
{
	// 문자열에 대한 문자 집합 분석
	// alphabets[char] = {strX에서의 발생횟수, strY에서의 발생횟수}
	map<wchar_t, pair<int, int>> alphabets;

	recordColumn(_osRecorder, _osTimeRecorder);

	for (int sizeNum = 0; sizeNum < _sizeIterNum; sizeNum++)
	{
		for (int testNum = 0; testNum < _testIterNum; testNum++)
		{
			// string size
			int size = (sizeNum + 1) * _sizeOffset;

			// random - Norm(0.0, 1,0)
			boost::uint32_t seed;
			seed = static_cast<boost::uint32_t>(time(0));
			seed = seed + (testNum * sizeNum) % 1049861;
			boost::variate_generator<boost::mt19937, boost::normal_distribution<>> rand(boost::mt19937(seed), boost::normal_distribution<>(0.0, 1.0));

			double* keys = new double[size];
			int* iKeys = new int[size];
			double rmin = 0, rmax = 0;

			for (int i = 0; i < size; i++)
			{
				keys[i] = rand();
				if (keys[i] < rmin) rmin = keys[i];
				if (keys[i] > rmax) rmax = keys[i];
			}

			// convert Norm(0.0, 1.0) to [1:_alphabetSize] integers
			double off = (abs(rmin) > abs(rmax) ? abs(rmin) : abs(rmax));
			double rOff = _alphabetSize / 2;
			off = rOff / off;

			for (int i = 0; i < size; i++)
			{
				keys[i] *= off;
				keys[i] += _alphabetSize / 2 + 1;

				// boundary out => make new char
				while (keys[i] < 1 || keys[i] > _alphabetSize)
				{
					double key = rand();
					keys[i] = key;
					keys[i] *= off;
					keys[i] += _alphabetSize / 2 + 1;
				}

				iKeys[i] = keys[i];
			}

			for (int i = 0; i < size; i++)
			{
				wchar_t wc = (wchar_t)iKeys[i];
				wstrX.push_back(wc);
				alphabets[wc].first++;
			}


			// string Y
			rmin = 0, rmax = 0;

			for (int i = 0; i < size; i++)
			{
				keys[i] = rand();
				if (keys[i] < rmin) rmin = keys[i];
				if (keys[i] > rmax) rmax = keys[i];
			}

			// convert Norm(0.0, 1.0) to [1:_alphabetSize] integers
			off = (abs(rmin) > abs(rmax) ? abs(rmin) : abs(rmax));
			rOff = _alphabetSize / 2;
			off = rOff / off;

			for (int i = 0; i < size; i++)
			{
				keys[i] *= off;
				keys[i] += _alphabetSize / 2 + 1;


				// boundary out => make new char
				while (keys[i] < 1 || keys[i] > _alphabetSize)
				{
					double key = rand();
					keys[i] = key;
					keys[i] *= off;
					keys[i] += _alphabetSize / 2 + 1;
				}

				iKeys[i] = keys[i];
			}

			for (int i = 0; i < size; i++)
			{
				wchar_t wc = (wchar_t)iKeys[i];
				wstrY.push_back(wc);
				alphabets[wc].second++;
			}
			delete[] keys;
			delete[] iKeys;

			if (wstrX.length() != size || wstrY.length() != size)
			{
				testNum--;
				alphabets.clear();
				continue;
			}

			// 실험 번호, 문자열 사이즈 기록
			_osRecorder << sizeNum * _testIterNum + testNum << '\t' << size << '\t';
			_osTimeRecorder << sizeNum * _testIterNum + testNum << '\t' << size << '\t';

			wreturnPacket caseResult;

			// LCS
			caseResult = LCS(wstrX, wstrY);
			if (!caseResult.isMaximal)
			{
				cout << "LCS" << caseResult.n << "is not maximal\n";
				exit(1);
			}
			record(_osRecorder, _osTimeRecorder, caseResult);


			// Y. Sakai
			caseResult = MCS_0(wstrX, wstrY);
			if (!caseResult.isMaximal)
			{
				cout << "Sakai" << caseResult.n << "is not maximal\n";
				exit(1);
			}
			record(_osRecorder, _osTimeRecorder, caseResult);

			// LeeP
			caseResult = MCS_1(wstrX, wstrY);
			record(_osRecorder, _osTimeRecorder, caseResult);
			if (!caseResult.isMaximal)
			{
				cout << "LeeP" << caseResult.n << "is not maximal\n";
				exit(1);
			}

			// LeeA
			caseResult = MCS_1_A(wstrX, wstrY);
			record(_osRecorder, _osTimeRecorder, caseResult);
			if (!caseResult.isMaximal)
			{
				cout << "LeeA" << caseResult.n << "is not maximal\n";
				exit(1);
			}

			// k Linear Searching
			for (int k = 1; k <= 128; k++)
			{
				caseResult = MCS_T1(wstrX, wstrY, k);
				if (!caseResult.isMaximal)
				{
					cout << "kc" << caseResult.n << "is not maximal\n";
					exit(1);
				}
				record(_osRecorder, _osTimeRecorder, caseResult);
			}

			// Except kr because it is not good algorithm
			//// k Range
			//for (int k = 1; k <= 128; k++)
			//{
			//	caseResult = MCS_T2(wstrX, wstrY, k);
			//	if (!caseResult.isMaximal)
			//	{
			//		cout << "kr" << caseResult.n << "is not maximal\n";
			//		exit(1);
			//	}
			//	record(_osRecorder, _osTimeRecorder, caseResult);
			//}

			// alphabets
			_osRecorder << alphabets.size() << '\t';
			vector<long double> possibilities[2];
			long double ePs[2]{ 0 };
			for (map<wchar_t, pair<int, int>>::iterator iter = alphabets.begin(); iter != alphabets.end(); iter++)
			{
				long double XP = (long double)(iter->second.first / (long double)size);
				long double YP = (long double)(iter->second.first / (long double)size);
				possibilities[0].push_back(XP);
				possibilities[1].push_back(YP);
				ePs[0] += XP;
				ePs[1] += YP;
			}
			ePs[0] /= (long double)alphabets.size();
			ePs[1] /= (long double)alphabets.size();

			long double vPs[2]{ 0 };
			for (int i = 0; i < possibilities[0].size(); i++)
			{
				vPs[0] += pow(ePs[0] - possibilities[0][i], 2);
				vPs[1] += pow(ePs[1] - possibilities[1][i], 2);
			}
			vPs[0] /= (long double)alphabets.size();
			vPs[1] /= (long double)alphabets.size();
			_osRecorder << vPs[0] << '\t' << vPs[1];

			_osRecorder << '\n';
			_osTimeRecorder << '\n';

			wstrX.clear();
			wstrY.clear();

			cout << _alphabetSize << ' ' << sizeNum * _testIterNum + testNum << '\n';
		}
	}
}