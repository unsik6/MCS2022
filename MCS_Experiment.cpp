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
// �̵���, & ����ä. (2022). �� �� �ش� ���� �κ� ������ ã�� ���� �˰���. ��������ȸ����, 49(7), 507-513.
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
// ������, & ������. (2022). �� ���ڿ��� �ش����κм����� ã�� ���ο� �˰���. 2022�� �ѱ�����Ʈ���������м���ȸ ����, 1212-1214.
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

	// push ����� iEnd6, jEnd6
	int pushed_iEnd = -1;
	int pushed_jEnd = -1;

	// check ��
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
	for (int alphaIdx = 0; alphaIdx < 11; alphaIdx)
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
	// �ð� ����
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	int LCSlength = 0, max;

	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();


	// matrix �����ϴ� �ݺ��� �ּ�ȭ�� ���� �۾�; �� ª�� ���� str1���� ����
	if (len1 > len2)
	{
		string tmp = _str1; _str1 = _str2; _str2 = tmp;
		int tmp2 = len1; len1 = len2; len2 = tmp2;
	}

	// LCS �˰����� ���� �ʿ��� ���
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

	// �ð� ����
	end_t = clock();
	result_t = (double)(end_t - start_t);

	// �Ҹ�
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
	// �ð� ����
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();

	string mcs = "";

	stack<ghs> stk;
	ghs init{ -1, -1, 0 };
	stk.push(init);

	// ������ ������ ��
	int iEnd = len1; int jEnd = len2;

	while (!stk.empty())
	{
		ghs curGHS = stk.top();
		stk.pop();

		// s�� ¦���̸� str1(X)���� ã�� ���̴�.
		if (curGHS.s % 2 == 0)
		{
			int searching_g = curGHS.g + curGHS.s / 2 + 1;

			// ���� str1(X)������ Ž�� ������ ������ ���� ���� �����Ѵ�.
			if (searching_g < iEnd)
			{
				ghs next_GHS{ curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);

				// ������ str2(Y)�� ���� ������ ���� ���ڸ� ã��
				// ���ο� ���빮�� offset�� ����
				if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd - 1)
				{
					ghs next_GHS{ searching_g, finding_h, 0 };
					stk.push(next_GHS);
				}
			}
			// str1(X)���� ���̻� Ž���� ���� ���� ������ ���빮�ڰ� ���۹���(��Ƽ��)�� �ƴ� ���
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
		// s�� Ȧ���̸� str2(Y)���� ã�� ���̴�.
		else
		{
			int searching_h = curGHS.h + (curGHS.s + 1) / 2;

			// ���� str2(Y)������ Ž�� ������ ������ ���� ���� �����Ѵ�.
			if (searching_h < jEnd)
			{
				ghs next_GHS{ curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

				// ������ str2(Y)�� ���� ������ ���� ���ڸ� ã��
				// ���ο� ���빮�� offset�� ����
				if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd - 1)
				{
					ghs next_GHS{ finding_g, searching_h, 0 };
					stk.push(next_GHS);
				}
			}
			// str1(X)���� ���̻� Ž���� ���� ���� ������ ���빮�ڰ� ���۹���(��Ƽ��)�� �ƴ� ���
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

	// �ð� ����
	end_t = clock();
	result_t = (double)(end_t - start_t);

	returnPacket result = { _str1, _str2, mcs, _str1.length() + _str2.length(), mcs.length(), isMCS(_str1, _str2, mcs), result_t };
	return result;
}

// Lee
returnPacket MCS_1(string _str1, string _str2)
{
	// �ð� ����
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();

	string mcs = "";

	stack<leeGhs> stk;
	leeGhs init{ -1, -1, 1, 1 };
	stk.push(init);

	// ������ ������ ��
	int iEnd = len1; int jEnd = len2;

	while (!stk.empty())
	{
		leeGhs curGHS = stk.top();
		stk.pop();

		// X�� Y���� Ž���� ������ ��ġ�� ���Ѵ�.
		int searching_g = curGHS.g + curGHS.gs;
		int searching_h = curGHS.h + curGHS.hs;

		// �� �� ��� Ž�� ������ ���� �ʴ� ���
		if (searching_g < iEnd && searching_h < jEnd)
		{
			int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);
			int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

			// ���� ������ ���ڰ� �ٸ� ���ڿ����� ��Ÿ�� ���
			if ((finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd - 1) &&
				(finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd - 1))
			{
				int minLen_by_g = min((iEnd - searching_g), (jEnd - finding_h));
				int minLen_by_h = min((iEnd - finding_g), (jEnd - searching_h));

				// X���� ã�� ���ڸ� �������� �� �� ������ �� �о����� ���
				if (minLen_by_g >= minLen_by_h)
				{
					leeGhs reGHS{ curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs };
					stk.push(reGHS);

					leeGhs next_GHS{ searching_g, finding_h, 1, 1 };
					stk.push(next_GHS);
				}
				// Y���� ã�� ���ڸ� �������� �� �� ������ �� �о����� ���
				else
				{
					leeGhs reGHS{ curGHS.g, curGHS.h, curGHS.gs, curGHS.hs + 1 };
					stk.push(reGHS);

					leeGhs next_GHS{ finding_g, searching_h, 1, 1 };
					stk.push(next_GHS);
				}
			}
			// X ������ ���ڰ� str2���� ��Ÿ�� ��츸
			else if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd - 1)
			{
				leeGhs reGHS{ curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs + 1 };
				stk.push(reGHS);

				leeGhs next_GHS{ searching_g, finding_h, 1, 1 };
				stk.push(next_GHS);
			}
			// Y ������ ���ڰ� str1���� ��Ÿ�� ��츸
			else if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd - 1)
			{
				leeGhs reGHS{ curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs + 1 };
				stk.push(reGHS);

				leeGhs next_GHS{ finding_g, searching_h, 1, 1 };
				stk.push(next_GHS);
			}
			// �Ѵ� �� ��Ÿ�� ���
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

	// �ð� ����
	end_t = clock();
	result_t = (double)(end_t - start_t);

	returnPacket result = { _str1, _str2, mcs, _str1.length() + _str2.length(), mcs.length(), isMCS(_str1, _str2, mcs), result_t };
	return result;
}

// experiment MCS
// k Match MCS(using Queue) - Fail
returnPacket MCS_T0(string _str1, string _str2, int _k)
{
	// �ð� ����
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

	// ������ ������ ��
	int iEnd = len1; int jEnd = len2;

	while (!stk.empty())
	{

		MyGhs curGHS = stk.top();

		bool isPoped = false;

		// (1) _k��ŭ Ž���� ������ _k ���� ���� ���� ������ ���� ���ڸ� ã�� ���� ���
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
		// (2) X���� _k��ŭ Ž���� �����ϱ� ���� Ž�������� ����� ���
		// (3) Y���� _k��ŭ Ž���� �����ϱ� ���� Ž�������� ����� ���
		// (2)�� (3)�� ��� Ž���� ����� �����̴�. �׷��� Ž���� ����Ǿ��� ��, iEnd�� jEnd�� �����ϴ���
		// new_iEnd(new_jEnd)�� old_iEnd(old_jEnd) ���̿� ���빮�ڰ� ������ ���·� ���� �� �ִ�.
		// ��, ���� ���빮�ڰ� ���� ���� �� �����Ƿ� (1)�� ���� ���� ó���ؾ� �Ѵ�.
		else if ((curGHS.s % 2 == 0 && curGHS.g + curGHS.s / 2 + 1 >= iEnd) ||
			(curGHS.s % 2 == 1 && curGHS.h + (curGHS.s + 1) / 2 >= jEnd))
		{
			// ���� ���빮�ڰ� ���ٸ� ��� ����.
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

		// s�� ¦���̸� str1(X)���� ã�� ���̴�.
		if (curGHS.s % 2 == 0)
		{
			int searching_g = curGHS.g + curGHS.s / 2 + 1;

			// ���� str1(X)������ Ž�� ������ ������ ���� ���� �����Ѵ�.
			if (searching_g < iEnd)
			{
				MyGhs next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);

				// ������ str2(Y)�� ���� ������ ���� ���ڸ� ã��
				// ���ο� ���빮�� offset�� ����
				if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd - 1)
				{
					// k ���� �ȿ��� ã�Ҵ�.
					if (searching_g <= curGHS.g + _k && finding_h <= curGHS.h + _k)
					{
						MyGhs next_GHS{ curGHS.id + 1, searching_g, finding_h, 0 };
						stk.push(next_GHS);
					}
					// k ������ ����� �����Ѵ�.
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
			// str1(X)���� ���̻� Ž���� ���� ���� ������ ���빮�ڰ� ���۹���(��Ƽ��)�� �ƴ� ���
			else if (curGHS.g >= 0)
			{
				mcs = _str1[curGHS.g] + mcs;

				iEnd = _str1.find_last_of(_str1[curGHS.g], iEnd - 1);
				jEnd = _str2.find_last_of(_str1[curGHS.g], jEnd - 1);

				if (iEnd == 0 || jEnd == 0 ||
					iEnd == -1 || jEnd == -1 ||
					iEnd == _str1.npos || jEnd == _str2.npos)
					break;

				// ���� wthQs[curGHS.id]�� ���� �ִ� �͵��� ������ g >= iEnd �Ǵ� h >= jEnd�� ���ۿ� ����.
				// �ֳ��ϸ� ���� �����ִٴ� ���� ��ҵ� Ž�� ������ _k�� ���� �ʾƼ� ���� �ʴ� ���̰�,
				// wthQs[curGHS.id]�� ����Ǿ� �ִٴ� ���� finding_g(h)�� _k������ �ѱ� �����̴�.
				// �־��� ��� _k * 2�� ȣ���ϹǷ� ��� �ð��� �ɸ���.
				if (wthQs.size() > curGHS.id)
					while (!wthQs[curGHS.id].empty())
						wthQs[curGHS.id].pop();
			}
		}
		// s�� Ȧ���̸� str2(Y)���� ã�� ���̴�.
		else
		{
			int searching_h = curGHS.h + (curGHS.s + 1) / 2;

			// ���� str2(Y)������ Ž�� ������ ������ ���� ���� �����Ѵ�.
			if (searching_h < jEnd)
			{
				MyGhs next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

				// ������ str2(Y)�� ���� ������ ���� ���ڸ� ã��
				// ���ο� ���빮�� offset�� ����
				if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd - 1)
				{
					// k ���� �ȿ��� ã�Ҵ�.
					if (searching_h <= curGHS.h + _k && finding_g <= curGHS.g + _k)
					{
						MyGhs next_GHS{ curGHS.id + 1, finding_g, searching_h, 0 };
						stk.push(next_GHS);
					}
					// k ������ ����� �����Ѵ�.
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
			// str1(X)���� ���̻� Ž���� ���� ���� ������ ���빮�ڰ� ���۹���(��Ƽ��)�� �ƴ� ���
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

	// �ð� ����
	end_t = clock();
	result_t = (double)(end_t - start_t);

	returnPacket result = { _str1, _str2, mcs, _str1.length() + _str2.length(), mcs.length(), isMCS(_str1, _str2, mcs), result_t };
	return result;
}

// k match compared MCS (using Deque) - Fail
returnPacket MCS_T1(string _str1, string _str2, int _k)
{
	// �ð� ����
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

	// ������ ������ ��
	int iEnd = len1; int jEnd = len2;

	while (!stk.empty())
	{

		MyGhs curGHS = stk.top();

		bool isPoped = false;

		// (1) _k��ŭ Ž���� ������ _k ���� ���� ���� ������ ���� ���ڸ� ã�� ���� ���
		if (curGHS.s >= _k * 2)
		{
			// ���� ���빮�ڰ� ���ٸ� ��� ����.
			if (wthQs.size() > curGHS.id)
			{
				// ���� �տ� �ִ� ���� ���� ������.
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
				// ���� �� ������ ���� Ž������ �Ѿ��.
				if (wthGhs.id >= 0)
				{
					// ���������� �� �� �̻� �����ִٸ� �񱳸� ���� ������.
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
					// �񱳸� ���� ���� ���� ���ٸ� �ռ� ������ ���� ����Ѵ�.
					// ���� �񱳸� ���� ���� ���� �ִٸ� ���� Ž�� ������ �� ũ�� ����� ���� �����ϰ�, �ٸ� �ϳ��� �ٽ� ť�� ����ִ´�.
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
		// (2) X���� _k��ŭ Ž���� �����ϱ� ���� Ž�������� ����� ���
		// (3) Y���� _k��ŭ Ž���� �����ϱ� ���� Ž�������� ����� ���
		// (2)�� (3)�� ��� Ž���� ����� �����̴�. �׷��� Ž���� ����Ǿ��� ��, iEnd�� jEnd�� �����ϴ���
		// new_iEnd(new_jEnd)�� old_iEnd(old_jEnd) ���̿� ���빮�ڰ� ������ ���·� ���� �� �ִ�.
		// ��, ���� ���빮�ڰ� ���� ���� �� �����Ƿ� (1)�� ���� ���� ó���ؾ� �Ѵ�.
		else if ((curGHS.s % 2 == 0 && curGHS.g + curGHS.s / 2 + 1 >= iEnd) ||
			(curGHS.s % 2 == 1 && curGHS.h + (curGHS.s + 1) / 2 >= jEnd))
		{
			// ���� ���빮�ڰ� ���ٸ� ��� ����.
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

		// s�� ¦���̸� str1(X)���� ã�� ���̴�.
		if (curGHS.s % 2 == 0)
		{
			int searching_g = curGHS.g + curGHS.s / 2 + 1;

			// ���� str1(X)������ Ž�� ������ ������ ���� ���� �����Ѵ�.
			if (searching_g < iEnd)
			{
				MyGhs next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);

				// ������ str2(Y)�� ���� ������ ���� ���ڸ� ã��
				// ���ο� ���빮�� offset�� ����
				if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd - 1)
				{
					// k ���� �ȿ��� ã�Ҵ�.
					if (searching_g <= curGHS.g + _k && finding_h <= curGHS.h + _k)
					{
						MyGhs next_GHS{ curGHS.id + 1, searching_g, finding_h, 0 };
						stk.push(next_GHS);
					}
					// k ������ ����� �����Ѵ�.
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
			// str1(X)���� ���̻� Ž���� ���� ���� ������ ���빮�ڰ� ���۹���(��Ƽ��)�� �ƴ� ���
			else if (curGHS.g >= 0)
			{
				mcs = _str1[curGHS.g] + mcs;

				iEnd = _str1.find_last_of(_str1[curGHS.g], iEnd - 1);
				jEnd = _str2.find_last_of(_str1[curGHS.g], jEnd - 1);

				if (iEnd == 0 || jEnd == 0 ||
					iEnd == -1 || jEnd == -1 ||
					iEnd == _str1.npos || jEnd == _str2.npos)
					break;

				// ���� wthQs[curGHS.id]�� ���� �ִ� �͵��� ������ g >= iEnd �Ǵ� h >= jEnd�� ���ۿ� ����.
				// �ֳ��ϸ� ���� �����ִٴ� ���� ��ҵ� Ž�� ������ _k�� ���� �ʾƼ� ���� �ʴ� ���̰�,
				// wthQs[curGHS.id]�� ����Ǿ� �ִٴ� ���� finding_g(h)�� _k������ �ѱ� �����̴�.
				// �־��� ��� _k * 2�� ȣ���ϹǷ� ��� �ð��� �ɸ���.
				if (wthQs.size() > curGHS.id)
					while (!wthQs[curGHS.id].empty())
						wthQs[curGHS.id].pop_front();
			}
		}
		// s�� Ȧ���̸� str2(Y)���� ã�� ���̴�.
		else
		{
			int searching_h = curGHS.h + (curGHS.s + 1) / 2;

			// ���� str2(Y)������ Ž�� ������ ������ ���� ���� �����Ѵ�.
			if (searching_h < jEnd)
			{
				MyGhs next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

				// ������ str2(Y)�� ���� ������ ���� ���ڸ� ã��
				// ���ο� ���빮�� offset�� ����
				if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd - 1)
				{
					// k ���� �ȿ��� ã�Ҵ�.
					if (searching_h <= curGHS.h + _k && finding_g <= curGHS.g + _k)
					{
						MyGhs next_GHS{ curGHS.id + 1, finding_g, searching_h, 0 };
						stk.push(next_GHS);
					}
					// k ������ ����� �����Ѵ�.
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
			// str1(X)���� ���̻� Ž���� ���� ���� ������ ���빮�ڰ� ���۹���(��Ƽ��)�� �ƴ� ���
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

	// �ð� ����
	end_t = clock();
	result_t = (double)(end_t - start_t);

	returnPacket result = { _str1, _str2, mcs, _str1.length() + _str2.length(), mcs.length(), isMCS(_str1, _str2, mcs), result_t };
	return result;
}

// k match Linear Max Search MCS (using vector)
// The Maximal Common Subsequence algorithms of Hyeonjun Shin and Jeong Seop Sim.
// ������, & ������. (2022). �� ���ڿ��� �ش����κм����� ã�� ���ο� �˰���. 2022�� �ѱ�����Ʈ���������м���ȸ ����, 1212-1214.
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
	// �ð� ����
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

	// ������ ������ ��
	iEnd5 = len1, jEnd5 = len2;

	// check ��
	int GHS_id = 1;

	while (!stk.empty())
	{
		PQGhs curGHS = stk.top();

		bool isWthPoped = false;

		// (1) _k��ŭ Ž���� ������ _k ���� ���� ���� ������ ���� ���ڸ� ã�� ���� ���
		// (2) X���� _k��ŭ Ž���� �����ϱ� ���� Ž�������� ����� ���
		// (3) Y���� _k��ŭ Ž���� �����ϱ� ���� Ž�������� ����� ���
		// (2)�� (3)�� ��� Ž���� ����� �����̴�. �׷��� Ž���� ����Ǿ��� ��, iEnd�� jEnd�� �����ϴ���
		// new_iEnd(new_jEnd)�� old_iEnd(old_jEnd) ���̿� ���빮�ڰ� ������ ���·� ���� �� �ִ�.
		// ��, ���� ���빮�ڰ� ���� ���� �� �����Ƿ� (1)�� ���� ���� ó���ؾ� �Ѵ�.
		if (curGHS.s >= _k * 2)
		{
			// ���� ���빮�ڰ� ���ٸ� ��� ����.
			if (wthPqs.size() > curGHS.id)
			{
				// ���� �� ���� ������ �Ϲ� Ž������ ����
				if (!wthPqs[curGHS.id].empty())
				{
					// ���� ���� ����
					// ���� ���Ŀ� ���� �ٿ�尡 �޶������� �� �޶��� ���� �ٿ�� ���Ŀ� ������ �ִ� ���� ģ����� ���� �ٿ�尡 �ٸ� ���� ģ������ ������ ���� ����. by Lemma 2
					// ���� Ž�� ������ �� �� ���� ������
					// ���� Ž���� �ϴ��� ���� �ٿ�带 �Ѿ������ ���� �ּ� ���� Ž�������� ���̳ʽ��̹Ƿ� ���� �ļ����̴�.
					PQGhs wthGHS = wthPqs[curGHS.id][0];
					for (int i = 1; i < wthPqs[curGHS.id].size(); i++)
					{
						PQGhs cur = wthPqs[curGHS.id][i];
						if (compGHS(cur, wthGHS))	// cur�� �� ��� ����� true ��ȯ
						{
							wthGHS = cur;
						}
					}

					// id�� �´� �� �����Ƿ� End �ٿ�常 �Ű澲�� �ȴ�.
					if (wthGHS.g < iEnd5 && wthGHS.h < jEnd5)
					{
						isWthPoped = true;

						// stk�� ���� ����ü �ʱ�ȭ
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
			// ���� ���빮�ڰ� ���ٸ� ��� ����.
			if (wthPqs.size() > curGHS.id)
			{
				// ���� �� ���� ������ �Ϲ� Ž������ ����
				if (!wthPqs[curGHS.id].empty())
				{
					// ���� ���� ����
					// ���� ���Ŀ� ���� �ٿ�尡 �޶������� �� �޶��� ���� �ٿ�� ���Ŀ� ������ �ִ� ���� ģ����� ���� �ٿ�尡 �ٸ� ���� ģ������ ������ ���� ����. by Lemma 2
					// ���� Ž�� ������ �� �� ���� ������
					// ���� Ž���� �ϴ��� ���� �ٿ�带 �Ѿ������ ���� �ּ� ���� Ž�������� ���̳ʽ��̹Ƿ� ���� �ļ����̴�.
					PQGhs wthGHS = wthPqs[curGHS.id][0];
					for (int i = 1; i < wthPqs[curGHS.id].size(); i++)
					{
						PQGhs cur = wthPqs[curGHS.id][i];
						if (compGHS(cur, wthGHS))	// cur�� �� ��� ����� true ��ȯ
						{
							wthGHS = cur;
						}
					}

					// id�� �´� �� �����Ƿ� End �ٿ�常 �Ű澲�� �ȴ�.
					if (wthGHS.g < iEnd5 && wthGHS.h < jEnd5)
					{
						isWthPoped = true;

						// stk�� ���� ����ü �ʱ�ȭ
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

		// s�� ¦���̸� str1(X)���� ã�� ���̴�.
		if (curGHS.s % 2 == 0)
		{
			int searching_g = curGHS.g + curGHS.s / 2 + 1;

			// ���� str1(X)������ Ž�� ������ ������ ���� ���� �����Ѵ�.
			if (searching_g < iEnd5)
			{
				PQGhs next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);

				// ������ str2(Y)�� ���� ������ ���� ���ڸ� ã��
				// ���ο� ���빮�� offset�� ����
				if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd5 - 1)
				{
					// k ���� �ȿ��� ã�Ҵ�.
					if (searching_g <= curGHS.g + _k && finding_h <= curGHS.h + _k)
					{
						PQGhs next_GHS{ GHS_id++, searching_g, finding_h, 0 };
						stk.push(next_GHS);
					}
					// k ������ ����� �����Ѵ�.
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
			// str1(X)���� ���̻� Ž���� ���� ���� ������ ���빮�ڰ� ���۹���(��Ƽ��)�� �ƴ� ���
			else if (curGHS.g >= 0)
			{
				mcs = _str1[curGHS.g] + mcs;

				iEnd5 = _str1.find_last_of(_str1[curGHS.g], iEnd5 - 1);
				jEnd5 = _str2.find_last_of(_str1[curGHS.g], jEnd5 - 1);
			}
		}
		// s�� Ȧ���̸� str2(Y)���� ã�� ���̴�.
		else
		{
			int searching_h = curGHS.h + (curGHS.s + 1) / 2;

			// ���� str2(Y)������ Ž�� ������ ������ ���� ���� �����Ѵ�.
			if (searching_h < jEnd5)
			{
				PQGhs next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

				// ������ str2(Y)�� ���� ������ ���� ���ڸ� ã��
				// ���ο� ���빮�� offset�� ����
				if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd5 - 1)
				{
					// k ���� �ȿ��� ã�Ҵ�.
					if (searching_h <= curGHS.h + _k && finding_g <= curGHS.g + _k)
					{
						PQGhs next_GHS{ GHS_id++, finding_g, searching_h, 0 };
						stk.push(next_GHS);
					}
					// k ������ ����� �����Ѵ�.
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
			// str1(X)���� ���̻� Ž���� ���� ���� ������ ���빮�ڰ� ���۹���(��Ƽ��)�� �ƴ� ���
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

	// �ð� ����
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

		// ���� �Ѵ� Y���� min���� ���� ���,
		//	�߰����� Lemma�� ���Ͽ� X������ ������� ����.
		// �׷��� �������, minLen�� ���� ���ķ� �����ϸ� ��.
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

		// ���� �Ѵ� X���� min���� ���� ���,
		//	�߰����� Lemma�� ���Ͽ� Y������ ������� ����.
		// �׷��� �������, minLen�� ���� ���ķ� �����ϸ� ��.
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
	// �ð� ����
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	// ���� ����
	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();

	string mcs = "";

	stack<PQGhs6> ghStk;
	PQGhs6 init{ 0, -1, -1, 0 };
	ghStk.push(init);
	// End Bound
	iEnd6 = len1; jEnd6 = len2;

	// X���� �߰��� ���빮�ڸ� �����ϴ� ��� Y���� �߰��� ���빮�ڸ� �����ϴ� ����� ��
	// �� gh���� �ΰ��� ���� �����Ѵ�.
	vector<pair< priority_queue<PQGhs6, vector<PQGhs6>, PQCompX >, priority_queue<PQGhs6, vector<PQGhs6>, PQCompY >>> wthTwoHeaps;
	priority_queue<PQGhs6, vector<PQGhs6>, PQCompX > initWthXPQ;
	priority_queue<PQGhs6, vector<PQGhs6>, PQCompY > initWthYPQ;

	wthTwoHeaps.push_back({ initWthXPQ, initWthYPQ });
	// �ߺ� ������ �ݺ����� �ʱ� ���� �迭
	// �� �迭�� ���������� �ش� ��ġ���� ������ ������ ��츦 ��� �ִ�.
	// Worst Case�� ��쿡�� ȿ���� �ִ�.
	vector<int> XQueries(len1, -5);
	vector<int> YQueries(len2, -5);

	// check ��
	int GHS_id = 1;
	long long pGHSID = 1;
	PQGhs6 lastGHS;

	while (!ghStk.empty())
	{
		PQGhs6 curGHS = ghStk.top();

		bool isWthPoped = false;

		if (curGHS.s >= _k * 2)
		{
			// ���� ���빮�ڰ� ���ٸ� ��� ����.
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
				// �Ѵ� �ִ� ������ ��Ÿ�� ���
				// �� ���� ���빮�� �� ���� Ž�� �Ÿ��� �� ��� �ϴ� ���� ����
				if (Xpop.id != -1 && Ypop.id != -1)
				{
					int minLen_by_X = min(iEnd6 - Xpop.g, jEnd6 - Xpop.h);
					int minLen_by_Y = min(iEnd6 - Ypop.g, jEnd6 - Ypop.h);

					// Y�� �� ū ��� Xpop�� Ypop �ֱ� => ���������� Xpop�� ���ÿ� ����
					if (minLen_by_X < minLen_by_Y)
					{
						Xpop = Ypop;
						wthTwoHeaps[curGHS.id].second.pop();
					}
					else    // X�� �� ū ���
						wthTwoHeaps[curGHS.id].first.pop();
				}
				// X������ �ּ� ������ ��Ÿ�� ���
				else if (Xpop.id != -1)
				{
					wthTwoHeaps[curGHS.id].first.pop();
				}
				// Y������ �ּ� ������ ��Ÿ�� ���
				else if (Ypop.id != -1)
				{
					Xpop = Ypop;
					wthTwoHeaps[curGHS.id].second.pop();
				}

				if (Xpop.id != -1)	// ��� ������ ���� ���ڰ� �����ϴ� ���
				{
					isWthPoped = true;

					// ghStk�� ���� ����ü �ʱ�ȭ
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
			// ���� ���빮�ڰ� ���ٸ� ��� ����.
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
				// �Ѵ� �ִ� ������ ��Ÿ�� ���
				// �� ���� ���빮�� �� ���� Ž�� �Ÿ��� �� ��� �ϴ� ���� ����
				if (Xpop.id != -1 && Ypop.id != -1)
				{
					int minLen_by_X = min(iEnd6 - Xpop.g, jEnd6 - Xpop.h);
					int minLen_by_Y = min(iEnd6 - Ypop.g, jEnd6 - Ypop.h);

					// Y�� �� ū ��� Xpop�� Ypop �ֱ� => ���������� Xpop�� ���ÿ� ����
					if (minLen_by_X < minLen_by_Y)
					{
						Xpop = Ypop;
						wthTwoHeaps[curGHS.id].second.pop();
					}
					else    // X�� �� ū ���
						wthTwoHeaps[curGHS.id].first.pop();
				}
				// X������ �ּ� ������ ��Ÿ�� ���
				else if (Xpop.id != -1)
				{
					wthTwoHeaps[curGHS.id].first.pop();
				}
				// Y������ �ּ� ������ ��Ÿ�� ���
				else if (Ypop.id != -1)
				{
					Xpop = Ypop;
					wthTwoHeaps[curGHS.id].second.pop();
				}

				if (Xpop.id != -1)	// ��� ������ ���� ���ڰ� �����ϴ� ���
				{
					isWthPoped = true;

					// ghStk�� ���� ����ü �ʱ�ȭ
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

		// s�� ¦���̸� str1(X)���� ã�� ���̴�.
		if (curGHS.s % 2 == 0)
		{
			int searching_g = curGHS.g + curGHS.s / 2 + 1;

			// ���� str1(X)������ Ž�� ������ ������ ���� ���� �����Ѵ�.
			if (searching_g < iEnd6)
			{
				PQGhs6 next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				ghStk.push(next_GHS);

				int finding_h;
				// ������ �صξ��� ������ ��ȿ�� ���
				if (XQueries[searching_g] < jEnd6 && XQueries[searching_g] > curGHS.h)
				{
					finding_h = XQueries[searching_g];
				}
				else  // ������ �صξ��� ������ ��ȿ���� ���� ���
				{
					finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);
					XQueries[searching_g] = finding_h;
				}

				// ������ str2(Y)�� ���� ������ ���� ���ڸ� ã��
				// ���ο� ���빮�� offset�� ����
				if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd6 - 1)
				{
					// k ���� �ȿ��� ã�Ҵ�.
					if (searching_g <= curGHS.g + _k && finding_h <= curGHS.h + _k)
					{
						PQGhs6 new_GHS{ GHS_id++, searching_g, finding_h, 0 };
						ghStk.push(new_GHS);
					}
					// k ������ ����� �����Ѵ�.
					else
					{
						while (wthTwoHeaps.size() <= curGHS.id)
						{
							priority_queue<PQGhs6, vector<PQGhs6>, PQCompX > newXPQ;
							priority_queue<PQGhs6, vector<PQGhs6>, PQCompY > newYPQ;

							wthTwoHeaps.push_back({ newXPQ, newYPQ });
						}

						PQGhs6 new_wth_GHS{ curGHS.id, searching_g, finding_h, curGHS.s, iEnd6, jEnd6, pGHSID++ };
						// X���� ã�����Ƿ�
						wthTwoHeaps[curGHS.id].first.push(new_wth_GHS);
					}
				}
			}
			// str1(X)���� ���̻� Ž���� ���� ���� ������ ���빮�ڰ� ���۹���(��Ƽ��)�� �ƴ� ���
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
				// ���� wthQs[curGHS.id]�� ���� �ִ� �͵��� ������ g >= iEnd �Ǵ� h >= jEnd�� ���ۿ� ����.
				// �ֳ��ϸ� ���� �����ִٴ� ���� ��ҵ� Ž�� ������ _k�� ���� �ʾƼ� ���� �ʴ� ���̰�,
				// wthQs[curGHS.id]�� ����Ǿ� �ִٴ� ���� finding_g(h)�� _k������ �ѱ� �����̴�.
				// �־��� ��� _k * 2�� ȣ���ϹǷ� O(_klog_k) �ð��� �ɸ���.
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
		// s�� Ȧ���̸� str2(Y)���� ã�� ���̴�.
		else
		{
			int searching_h = curGHS.h + (curGHS.s + 1) / 2;

			// ���� str2(Y)������ Ž�� ������ ������ ���� ���� �����Ѵ�.
			if (searching_h < jEnd6)
			{
				PQGhs6 next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				ghStk.push(next_GHS);

				int finding_g;
				// ������ �صξ��� ������ ��ȿ�� ���
				if (YQueries[searching_h] < iEnd6 && YQueries[searching_h] > curGHS.g)
				{
					finding_g = YQueries[searching_h];
				}
				else  // ������ �صξ��� ������ ��ȿ���� ���� ���
				{
					finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);
					YQueries[searching_h] = finding_g;
				}
				// ������ str2(Y)�� ���� ������ ���� ���ڸ� ã��
				// ���ο� ���빮�� offset�� ����
				if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd6 - 1)
				{
					// k ���� �ȿ��� ã�Ҵ�.
					if (searching_h <= curGHS.h + _k && finding_g <= curGHS.g + _k)
					{
						PQGhs6 new_GHS{ GHS_id++, finding_g, searching_h, 0 };
						ghStk.push(new_GHS);
					}
					// k ������ ����� �����Ѵ�.
					else
					{
						while (wthTwoHeaps.size() <= curGHS.id)
						{
							priority_queue<PQGhs6, vector<PQGhs6>, PQCompX > newXPQ;
							priority_queue<PQGhs6, vector<PQGhs6>, PQCompY > newYPQ;

							wthTwoHeaps.push_back({ newXPQ, newYPQ });
						}
						PQGhs6 new_wth_GHS{ curGHS.id, finding_g, searching_h, curGHS.s, iEnd6, jEnd6, pGHSID++ };
						// Y���� ã�����Ƿ�
						wthTwoHeaps[curGHS.id].second.push(new_wth_GHS);
					}

				}
			}
			// str1(X)���� ���̻� Ž���� ���� ���� ������ ���빮�ڰ� ���۹���(��Ƽ��)�� �ƴ� ���
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

				// ���� wthQs[curGHS.id]�� ���� �ִ� �͵��� ������ g >= iEnd �Ǵ� h >= jEnd�� ���ۿ� ����.
				// �ֳ��ϸ� ���� �����ִٴ� ���� ��ҵ� Ž�� ������ _k�� ���� �ʾƼ� ���� �ʴ� ���̰�,
				// wthQs[curGHS.id]�� ����Ǿ� �ִٴ� ���� finding_g(h)�� _k������ �ѱ� �����̴�.
				// �־��� ��� _k * 2�� ȣ���ϹǷ� O(_klog_k) �ð��� �ɸ���.
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

	// �ð� ����
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
	// �ð� ����
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	int LCSlength = 0, max;

	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();


	// matrix �����ϴ� �ݺ��� �ּ�ȭ�� ���� �۾�; �� ª�� ���� str1���� ����
	if (len1 > len2)
	{
		wstring tmp = _str1; _str1 = _str2; _str2 = tmp;
		int tmp2 = len1; len1 = len2; len2 = tmp2;
	}

	// LCS �˰����� ���� �ʿ��� ���
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

	// �ð� ����
	end_t = clock();
	result_t = (double)(end_t - start_t);

	// �Ҹ�
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
	// �ð� ����
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();

	wstring mcs = L"";

	stack<ghs> stk;
	ghs init{ -1, -1, 0 };
	stk.push(init);

	// ������ ������ ��
	int iEnd = len1; int jEnd = len2;

	while (!stk.empty())
	{
		ghs curGHS = stk.top();
		stk.pop();

		// s�� ¦���̸� str1(X)���� ã�� ���̴�.
		if (curGHS.s % 2 == 0)
		{
			int searching_g = curGHS.g + curGHS.s / 2 + 1;

			// ���� str1(X)������ Ž�� ������ ������ ���� ���� �����Ѵ�.
			if (searching_g < iEnd)
			{
				ghs next_GHS{ curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);

				// ������ str2(Y)�� ���� ������ ���� ���ڸ� ã��
				// ���ο� ���빮�� offset�� ����
				if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd - 1)
				{
					ghs next_GHS{ searching_g, finding_h, 0 };
					stk.push(next_GHS);
				}
			}
			// str1(X)���� ���̻� Ž���� ���� ���� ������ ���빮�ڰ� ���۹���(��Ƽ��)�� �ƴ� ���
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
		// s�� Ȧ���̸� str2(Y)���� ã�� ���̴�.
		else
		{
			int searching_h = curGHS.h + (curGHS.s + 1) / 2;

			// ���� str2(Y)������ Ž�� ������ ������ ���� ���� �����Ѵ�.
			if (searching_h < jEnd)
			{
				ghs next_GHS{ curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

				// ������ str2(Y)�� ���� ������ ���� ���ڸ� ã��
				// ���ο� ���빮�� offset�� ����
				if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd - 1)
				{
					ghs next_GHS{ finding_g, searching_h, 0 };
					stk.push(next_GHS);
				}
			}
			// str1(X)���� ���̻� Ž���� ���� ���� ������ ���빮�ڰ� ���۹���(��Ƽ��)�� �ƴ� ���
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

	// �ð� ����
	end_t = clock();
	result_t = (double)(end_t - start_t);

	wreturnPacket result = { _str1, _str2, mcs, _str1.length() + _str2.length(), mcs.length(), isMCS(_str1, _str2, mcs), result_t };
	return result;
}
wreturnPacket MCS_1(wstring _str1, wstring _str2)
{
	// �ð� ����
	clock_t start_t, end_t;
	double result_t;
	start_t = clock();

	int len1 = (int)_str1.length();
	int len2 = (int)_str2.length();

	wstring mcs = L"";

	stack<leeGhs> stk;
	leeGhs init{ -1, -1, 1, 1 };
	stk.push(init);

	// ������ ������ ��
	int iEnd = len1; int jEnd = len2;

	while (!stk.empty())
	{
		leeGhs curGHS = stk.top();
		stk.pop();

		// X�� Y���� Ž���� ������ ��ġ�� ���Ѵ�.
		int searching_g = curGHS.g + curGHS.gs;
		int searching_h = curGHS.h + curGHS.hs;

		// �� �� ��� Ž�� ������ ���� �ʴ� ���
		if (searching_g < iEnd && searching_h < jEnd)
		{
			int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);
			int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

			// ���� ������ ���ڰ� �ٸ� ���ڿ����� ��Ÿ�� ���
			if ((finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd - 1) &&
				(finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd - 1))
			{
				int minLen_by_g = min((iEnd - searching_g), (jEnd - finding_h));
				int minLen_by_h = min((iEnd - finding_g), (jEnd - searching_h));

				// X���� ã�� ���ڸ� �������� �� �� ������ �� �о����� ���
				if (minLen_by_g >= minLen_by_h)
				{
					leeGhs reGHS{ curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs };
					stk.push(reGHS);

					leeGhs next_GHS{ searching_g, finding_h, 1, 1 };
					stk.push(next_GHS);
				}
				// Y���� ã�� ���ڸ� �������� �� �� ������ �� �о����� ���
				else
				{
					leeGhs reGHS{ curGHS.g, curGHS.h, curGHS.gs, curGHS.hs + 1 };
					stk.push(reGHS);

					leeGhs next_GHS{ finding_g, searching_h, 1, 1 };
					stk.push(next_GHS);
				}
			}
			// X ������ ���ڰ� str2���� ��Ÿ�� ��츸
			else if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd - 1)
			{
				leeGhs reGHS{ curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs + 1 };
				stk.push(reGHS);

				leeGhs next_GHS{ searching_g, finding_h, 1, 1 };
				stk.push(next_GHS);
			}
			// Y ������ ���ڰ� str1���� ��Ÿ�� ��츸
			else if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd - 1)
			{
				leeGhs reGHS{ curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs + 1 };
				stk.push(reGHS);

				leeGhs next_GHS{ finding_g, searching_h, 1, 1 };
				stk.push(next_GHS);
			}
			// �Ѵ� �� ��Ÿ�� ���
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

	// �ð� ����
	end_t = clock();
	result_t = (double)(end_t - start_t);

	wreturnPacket result = { _str1, _str2, mcs, _str1.length() + _str2.length(), mcs.length(), isMCS(_str1, _str2, mcs), result_t };
	return result;
}
wreturnPacket MCS_T2(wstring _str1, wstring _str2, int _k)
{
	// �ð� ����
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

	// ������ ������ ��
	iEnd5 = len1, jEnd5 = len2;

	// check ��
	int GHS_id = 1;

	while (!stk.empty())
	{
		PQGhs curGHS = stk.top();

		bool isWthPoped = false;

		// (1) _k��ŭ Ž���� ������ _k ���� ���� ���� ������ ���� ���ڸ� ã�� ���� ���
		// (2) X���� _k��ŭ Ž���� �����ϱ� ���� Ž�������� ����� ���
		// (3) Y���� _k��ŭ Ž���� �����ϱ� ���� Ž�������� ����� ���
		// (2)�� (3)�� ��� Ž���� ����� �����̴�. �׷��� Ž���� ����Ǿ��� ��, iEnd�� jEnd�� �����ϴ���
		// new_iEnd(new_jEnd)�� old_iEnd(old_jEnd) ���̿� ���빮�ڰ� ������ ���·� ���� �� �ִ�.
		// ��, ���� ���빮�ڰ� ���� ���� �� �����Ƿ� (1)�� ���� ���� ó���ؾ� �Ѵ�.
		if (curGHS.s >= _k * 2)
		{
			// ���� ���빮�ڰ� ���ٸ� ��� ����.
			if (wthPqs.size() > curGHS.id)
			{
				// ���� �� ���� ������ �Ϲ� Ž������ ����
				if (!wthPqs[curGHS.id].empty())
				{
					// ���� ���� ����
					// ���� ���Ŀ� ���� �ٿ�尡 �޶������� �� �޶��� ���� �ٿ�� ���Ŀ� ������ �ִ� ���� ģ����� ���� �ٿ�尡 �ٸ� ���� ģ������ ������ ���� ����. by Lemma 2
					// ���� Ž�� ������ �� �� ���� ������
					// ���� Ž���� �ϴ��� ���� �ٿ�带 �Ѿ������ ���� �ּ� ���� Ž�������� ���̳ʽ��̹Ƿ� ���� �ļ����̴�.
					PQGhs wthGHS = wthPqs[curGHS.id][0];
					for (int i = 1; i < wthPqs[curGHS.id].size(); i++)
					{
						PQGhs cur = wthPqs[curGHS.id][i];
						if (compGHS(cur, wthGHS))	// cur�� �� ��� ����� true ��ȯ
						{
							wthGHS = cur;
						}
					}

					// id�� �´� �� �����Ƿ� End �ٿ�常 �Ű澲�� �ȴ�.
					if (wthGHS.g < iEnd5 && wthGHS.h < jEnd5)
					{
						isWthPoped = true;

						// stk�� ���� ����ü �ʱ�ȭ
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
			// ���� ���빮�ڰ� ���ٸ� ��� ����.
			if (wthPqs.size() > curGHS.id)
			{
				// ���� �� ���� ������ �Ϲ� Ž������ ����
				if (!wthPqs[curGHS.id].empty())
				{
					// ���� ���� ����
					// ���� ���Ŀ� ���� �ٿ�尡 �޶������� �� �޶��� ���� �ٿ�� ���Ŀ� ������ �ִ� ���� ģ����� ���� �ٿ�尡 �ٸ� ���� ģ������ ������ ���� ����. by Lemma 2
					// ���� Ž�� ������ �� �� ���� ������
					// ���� Ž���� �ϴ��� ���� �ٿ�带 �Ѿ������ ���� �ּ� ���� Ž�������� ���̳ʽ��̹Ƿ� ���� �ļ����̴�.
					PQGhs wthGHS = wthPqs[curGHS.id][0];
					for (int i = 1; i < wthPqs[curGHS.id].size(); i++)
					{
						PQGhs cur = wthPqs[curGHS.id][i];
						if (compGHS(cur, wthGHS))	// cur�� �� ��� ����� true ��ȯ
						{
							wthGHS = cur;
						}
					}

					// id�� �´� �� �����Ƿ� End �ٿ�常 �Ű澲�� �ȴ�.
					if (wthGHS.g < iEnd5 && wthGHS.h < jEnd5)
					{
						isWthPoped = true;

						// stk�� ���� ����ü �ʱ�ȭ
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

		// s�� ¦���̸� str1(X)���� ã�� ���̴�.
		if (curGHS.s % 2 == 0)
		{
			int searching_g = curGHS.g + curGHS.s / 2 + 1;

			// ���� str1(X)������ Ž�� ������ ������ ���� ���� �����Ѵ�.
			if (searching_g < iEnd5)
			{
				PQGhs next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);

				// ������ str2(Y)�� ���� ������ ���� ���ڸ� ã��
				// ���ο� ���빮�� offset�� ����
				if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd5 - 1)
				{
					// k ���� �ȿ��� ã�Ҵ�.
					if (searching_g <= curGHS.g + _k && finding_h <= curGHS.h + _k)
					{
						PQGhs next_GHS{ GHS_id++, searching_g, finding_h, 0 };
						stk.push(next_GHS);
					}
					// k ������ ����� �����Ѵ�.
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
			// str1(X)���� ���̻� Ž���� ���� ���� ������ ���빮�ڰ� ���۹���(��Ƽ��)�� �ƴ� ���
			else if (curGHS.g >= 0)
			{
				mcs = _str1[curGHS.g] + mcs;

				iEnd5 = _str1.find_last_of(_str1[curGHS.g], iEnd5 - 1);
				jEnd5 = _str2.find_last_of(_str1[curGHS.g], jEnd5 - 1);
			}
		}
		// s�� Ȧ���̸� str2(Y)���� ã�� ���̴�.
		else
		{
			int searching_h = curGHS.h + (curGHS.s + 1) / 2;

			// ���� str2(Y)������ Ž�� ������ ������ ���� ���� �����Ѵ�.
			if (searching_h < jEnd5)
			{
				PQGhs next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
				stk.push(next_GHS);

				int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

				// ������ str2(Y)�� ���� ������ ���� ���ڸ� ã��
				// ���ο� ���빮�� offset�� ����
				if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd5 - 1)
				{
					// k ���� �ȿ��� ã�Ҵ�.
					if (searching_h <= curGHS.h + _k && finding_g <= curGHS.g + _k)
					{
						PQGhs next_GHS{ GHS_id++, finding_g, searching_h, 0 };
						stk.push(next_GHS);
					}
					// k ������ ����� �����Ѵ�.
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
			// str1(X)���� ���̻� Ž���� ���� ���� ������ ���빮�ڰ� ���۹���(��Ƽ��)�� �ƴ� ���
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

	// �ð� ����
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
void record(ofstream& _os, returnPacket& _rp)
{
	_os << _rp.cssLen << '\t' << _rp.isMaximal << '\t' << (double)((double)_rp.cssLen / (double)_rp.n * 100 * 2) << '\t' << _rp.time << '\t';
}
void record(ofstream& _os, wreturnPacket& _rp)
{
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
	// ���ڿ��� ���� ���� ���� �м�
	// alphabets[char] = {strX������ �߻�Ƚ��, strY������ �߻�Ƚ��}
	map<char, pair<int, int>> alphabets;

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

			// ���� ��ȣ, ���ڿ� ������ ���
			_osRecorder << sizeNum * _testIterNum + testNum << '\t' << size << '\t';
			_osRecorder << alphabets.size();

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
			for (map<char, pair<int, int>>::iterator iter = alphabets.begin(); iter != alphabets.end(); iter++)
				_osRecorder << iter->first << '\t' << iter->second.first << '\t' << iter->second.second << '\t';

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
	// ���ڿ��� ���� ���� ���� �м�
	// alphabets[char] = {strX������ �߻�Ƚ��, strY������ �߻�Ƚ��}
	map<wchar_t, pair<int, int>> alphabets;

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

			// ���� ��ȣ, ���ڿ� ������ ���
			_osRecorder << sizeNum * _testIterNum + testNum << '\t' << size << '\t';
			_osRecorder << alphabets.size();

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
			for (map<wchar_t, pair<int, int>>::iterator iter = alphabets.begin(); iter != alphabets.end(); iter++)
				_osRecorder << iter->first << '\t' << iter->second.first << '\t' << iter->second.second << '\t';

			_osRecorder << '\n';

			wstrX.clear();
			wstrY.clear();

			cout << _alphabetSize << ' ' << sizeNum * _testIterNum + testNum << '\n';
		}
	}
}