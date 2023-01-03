#pragma once
#include <iostream>
#include <string.h>

#include <stack>
#include <vector>

namespace MCS
{
	// the Longest Common Subsequence
	std::string LCS(std::string _str1, std::string _str2)
	{
		int LCSlength = 0, max;

		int len1 = (int)_str1.length();
		int len2 = (int)_str2.length();


		// Optimization: To make less iteration to build the Dynamic Programming matrix.
		if (len1 > len2)
		{
			std::string tmp = _str1; _str1 = _str2; _str2 = tmp;
			int tmp2 = len1; len1 = len2; len2 = tmp2;
		}

		// the Dynamic Prograamming Matrix
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
		std::string lcs = "";
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

		return lcs;
	}

	// The Maximal Common Subsequence algorithms of Yoshifumi Sakai
	// Sakai, Y. (2019). Maximal common subsequence algorithms. Theoretical Computer Science, 793, 132-139.
	// https://www.sciencedirect.com/science/article/pii/S0304397519304074
	
	// The data type to store the number of finding a common character and each position of common characters in each string.
	struct S_GHS
	{
		int g = -1;
		int h = -1;
		int s = -1;
	};

	std::string Y_Sakai_MCS(std::string _str1, std::string _str2)
	{
		int len1 = (int)_str1.length();
		int len2 = (int)_str2.length();

		std::string mcs = "";

		std::stack<S_GHS> stk;
		S_GHS init{ -1, -1, 0 };
		stk.push(init);

		// Indices of each right border of two input strings.
		int iEnd = len1; int jEnd = len2;

		while (!stk.empty())
		{
			S_GHS curGHS = stk.top();
			stk.pop();

			// If 's' is even, find a common character in '_str1'.
			if (curGHS.s % 2 == 0)
			{
				int searching_g = curGHS.g + curGHS.s / 2 + 1;

				// The character, target to find a common character, is still in the range to find a common character.
				if (searching_g < iEnd)
				{
					S_GHS next_GHS{ curGHS.g, curGHS.h, curGHS.s + 1 };
					stk.push(next_GHS);

					int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);

					// Find a common character in '_str2'.
					// Change each left border of two input strings to the new common character.
					if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd - 1)
					{
						S_GHS next_GHS{ searching_g, finding_h, 0 };
						stk.push(next_GHS);
					}
				}
				// There is no more character in the range of '_str1', and the top of stack does not store sentinel position.
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
			// If 's' is odd, find a common character in '_str2'.
			else
			{
				int searching_h = curGHS.h + (curGHS.s + 1) / 2;

				// The character, target to find a common character, is still in the range to find a common character.
				if (searching_h < jEnd)
				{
					S_GHS next_GHS{ curGHS.g, curGHS.h, curGHS.s + 1 };
					stk.push(next_GHS);

					int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

					// Find a common character in '_str1'.
					// Change each left border of two input strings to the new common character.
					if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd - 1)
					{
						S_GHS next_GHS{ finding_g, searching_h, 0 };
						stk.push(next_GHS);
					}
				}
				// There is no more character in the range of '_str2', and the top of stack does not store sentinel position.
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
		return mcs;
	}

	
	// The Maximal Common Subsequence algorithms of DongYeop Lee and Joong Chae Na.
	// 이동엽, & 나중채. (2022). 더 긴 극대 공통 부분 서열을 찾기 위한 알고리즘. 정보과학회논문지, 49(7), 507-513.
	// DongYeop Lee, & Joong Chae Na. (2022), An Improved Algorithm for Finding a Longer Maximal Common Subsequence. Journal of KIISE, 49(7), 507-513.
	// https://www.dbpia.co.kr/pdf/pdfView.do?nodeId=NODE11100316&googleIPSandBox=false&mark=0&useDate=&ipRange=false&accessgl=Y&language=ko_KR&hasTopBanner=true

	// The data type to store the number of finding a common character and each position of common characters in each string.
	struct LN_GHS
	{
		int g = 0;
		int h = 0;
		int gs = 0;
		int hs = 0;
	};
	std::string Lee_Na_MCS(std::string _str1, std::string _str2)
	{
		int len1 = (int)_str1.length();
		int len2 = (int)_str2.length();

		std::string mcs = "";

		std::stack<LN_GHS> stk;
		LN_GHS init{ -1, -1, 1, 1 };
		stk.push(init);

		// indices of each right border of two input strings.
		int iEnd = len1; int jEnd = len2;

		while (!stk.empty())
		{
			LN_GHS curGHS = stk.top();
			stk.pop();

			// Get the positions to find a common character in both strings.
			int searching_g = curGHS.g + curGHS.gs;
			int searching_h = curGHS.h + curGHS.hs;

			// The characters, target to find a common character, are still in the range to find a common character in both strings.
			if (searching_g < iEnd && searching_h < jEnd)
			{
				int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);
				int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

				// Find common characters in both strings.
				if ((finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd - 1) &&
					(finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd - 1))
				{
					int minLen_by_g = std::min((iEnd - searching_g), (jEnd - finding_h));
					int minLen_by_h = std::min((iEnd - finding_g), (jEnd - searching_h));

					// The max number of trying to find next common characters by the common character found using a character in '_str1' is more greater than another.
					if (minLen_by_g >= minLen_by_h)
					{
						LN_GHS reGHS{ curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs };
						stk.push(reGHS);

						LN_GHS next_GHS{ searching_g, finding_h, 1, 1 };
						stk.push(next_GHS);
					}
					// The max number of trying to find next common characters by the common character found using a character in '_str2' is more greater than another.
					else
					{
						LN_GHS reGHS{ curGHS.g, curGHS.h, curGHS.gs, curGHS.hs + 1 };
						stk.push(reGHS);

						LN_GHS next_GHS{ finding_g, searching_h, 1, 1 };
						stk.push(next_GHS);
					}
				}
				// A common character is found in not '_str2' but '_str1'.
				else if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd - 1)
				{
					LN_GHS reGHS{ curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs + 1 };
					stk.push(reGHS);

					LN_GHS next_GHS{ searching_g, finding_h, 1, 1 };
					stk.push(next_GHS);
				}
				// A common character is found in not '_str1' but '_str2'.
				else if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd - 1)
				{
					LN_GHS reGHS{ curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs + 1 };
					stk.push(reGHS);

					LN_GHS next_GHS{ finding_g, searching_h, 1, 1 };
					stk.push(next_GHS);
				}
				// There is no common character found in this stage.
				else
				{
					LN_GHS reGHS{ curGHS.g, curGHS.h, curGHS.gs + 1, curGHS.hs + 1 };
					stk.push(reGHS);
				}
			}
			// There is no more character in the range of both strings, and the top of stack does not store sentinel position.
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
		return mcs;
	}


	// The Maximal Common Subsequence algorithms of Hyeonjun Shin and Jeong Seop Sim.
	// 신현준, & 심정섭. (2022). 두 문자열의 극대공통부분서열을 찾는 새로운 알고리즘. 2022년 한국소프트웨어종합학술대회 논문집, 1212-1214.
	// Hyeonjun Shin, & Jeon Seop Sim. A New Algorithm of Finding a Maximal Common Subsequence of Two Strings. Korea Software Congress 2022, 1212-1214.

	// Indices of each right border of two input strings.
	int iEnd_SS, jEnd_SS;
	// The data type to store the number of finding a common character and each position of common characters in each string.
	struct SS_GHS
	{
		int id = -1;
		int g = -1;
		int h = -1;
		int s = 0;
	};
	// The function to compare two common characters returns TRUE if the first argument(SS_GHS) makes more number of finding next common characters than the second one.
	bool compGHS(const SS_GHS _a, const SS_GHS _b)
	{
		int min_by_a = iEnd_SS - _a.g; int max_by_a = jEnd_SS - _a.h;
		int min_by_b = iEnd_SS - _b.g; int max_by_b = jEnd_SS - _b.h;

		if (min_by_a > max_by_a) std::swap(min_by_a, max_by_a);
		if (min_by_b > max_by_b) std::swap(min_by_b, max_by_b);

		if (min_by_a > min_by_b)
			return true;
		else if (min_by_a < min_by_b)
			return false;
		else
		{
			return max_by_a >= max_by_b;
		}
	}
	std::string Shin_Sim_MCS(std::string _str1, std::string _str2, int _k)
	{
		int len1 = (int)_str1.length();
		int len2 = (int)_str2.length();

		std::string mcs = "";

		std::stack<SS_GHS> stk;
		SS_GHS init{ 0, -1, -1, 0 };
		stk.push(init);

		std::vector<std::vector<SS_GHS>> wthArrays;
		std::vector<SS_GHS> init_wthArr;
		wthArrays.push_back(init_wthArr);

		// Initiate indices of each right border of two input strings.
		iEnd_SS = len1, jEnd_SS = len2;

		// The index of arrays to store pending common characters.
		int GHS_id = 1;

		while (!stk.empty())
		{
			SS_GHS curGHS = stk.top();

			bool isWthPoped = false;

			// The stage to use pending common characters.
			// Case 1: There is no common characters whose distance is less than '_k'.
			if (curGHS.s >= _k * 2)
			{
				// There is no common pending character.
				if (wthArrays.size() > curGHS.id)
				{
					// There is no pending common character.
					if (!wthArrays[curGHS.id].empty())
					{
						// Find the pending common character that makes makes more number of finding next common characters.
						SS_GHS wthGHS = wthArrays[curGHS.id][0];
						for (int i = 1; i < wthArrays[curGHS.id].size(); i++)
						{
							SS_GHS cur = wthArrays[curGHS.id][i];
							if (compGHS(cur, wthGHS))
								wthGHS = cur;
						}

						// Check whether the common character is in the range of both strings.
						if (wthGHS.g < iEnd_SS && wthGHS.h < jEnd_SS)
						{
							isWthPoped = true;

							// push the information of the common character choosed to the stack.
							wthGHS.id = GHS_id++;
							wthGHS.s = 0;
							stk.push(wthGHS);
						}
						// If the common character is not in the range of both strings, then all common characters in the array is not in the range.
						// So, erase all.
						else
						{
							wthArrays[curGHS.id].clear();
						}
					}
				}
			}
			// Case 2: There is no more characters that can be target of finding a common character before '_k' distance.
			// Do same thing in Case 1.
			else if ((curGHS.s % 2 == 0 && curGHS.g + curGHS.s / 2 + 1 >= iEnd_SS) || (curGHS.s % 2 == 1 && curGHS.h + (curGHS.s + 1) / 2 >= jEnd_SS))
			{
				if (wthArrays.size() > curGHS.id)
				{
					if (!wthArrays[curGHS.id].empty())
					{
						SS_GHS wthGHS = wthArrays[curGHS.id][0];
						for (int i = 1; i < wthArrays[curGHS.id].size(); i++)
						{
							SS_GHS cur = wthArrays[curGHS.id][i];
							if (compGHS(cur, wthGHS))
								wthGHS = cur;
						}

						if (wthGHS.g < iEnd_SS && wthGHS.h < jEnd_SS)
						{
							isWthPoped = true;

							wthGHS.id = GHS_id++;
							wthGHS.s = 0;

							stk.push(wthGHS);
						}
						else
						{
							wthArrays[curGHS.id].clear();
						}
					}
				}
			}

			// If a pending common character is used, then change 'curGHS'.
			if (isWthPoped)
				curGHS = stk.top();

			stk.pop();

			// If 's' is even, find a common character in '_str1'.
			if (curGHS.s % 2 == 0)
			{
				int searching_g = curGHS.g + curGHS.s / 2 + 1;

				// The character, target to find a common character, is still in the range to find a common character.
				if (searching_g < iEnd_SS)
				{
					SS_GHS next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
					stk.push(next_GHS);

					int finding_h = _str2.find_first_of(_str1[searching_g], curGHS.h + 1);

					// Find a common character in '_str2'.
					// Change each left border of two input strings to the new common character.
					if (finding_h != _str2.npos && finding_h != -1 && finding_h <= jEnd_SS - 1)
					{
						// Find a common character in '_k' distance.
						// Choose the common charcter found as a charcter of MCS.
						if (searching_g <= curGHS.g + _k && finding_h <= curGHS.h + _k)
						{
							SS_GHS next_GHS{ GHS_id++, searching_g, finding_h, 0 };
							stk.push(next_GHS);
						}
						// Find a common character in more distance than '_k'.
						else
						{
							// If there is no array indexed by current left border character, then make new array.
							while (wthArrays.size() <= curGHS.id)
							{
								std::vector<SS_GHS> newQue;
								wthArrays.push_back(newQue);
							}

							// Push the common character found in this stage to the array indexed by current left border character.
							SS_GHS wth_GHS{ curGHS.id, searching_g, finding_h, curGHS.s };
							wthArrays[curGHS.id].push_back(wth_GHS);
						}
					}
				}
				// There is no more character in the range of '_str1', and the top of stack does not store sentinel position.
				else if (curGHS.g >= 0)
				{
					mcs = _str1[curGHS.g] + mcs;

					iEnd_SS = _str1.find_last_of(_str1[curGHS.g], iEnd_SS - 1);
					jEnd_SS = _str2.find_last_of(_str1[curGHS.g], jEnd_SS - 1);
				}
			}
			// If 's' is odd, find a common character in '_str2'.
			else
			{
				int searching_h = curGHS.h + (curGHS.s + 1) / 2;

				// The character, target to find a common character, is still in the range to find a common character.
				if (searching_h < jEnd_SS)
				{
					SS_GHS next_GHS{ curGHS.id, curGHS.g, curGHS.h, curGHS.s + 1 };
					stk.push(next_GHS);

					int finding_g = _str1.find_first_of(_str2[searching_h], curGHS.g + 1);

					// Find a common character in '_str1'.
					// Change each left border of two input strings to the new common character.
					if (finding_g != _str1.npos && finding_g != -1 && finding_g <= iEnd_SS - 1)
					{
						// Find a common character in '_k' distance.
						// Choose the common charcter found as a charcter of MCS.
						if (searching_h <= curGHS.h + _k && finding_g <= curGHS.g + _k)
						{
							SS_GHS next_GHS{ GHS_id++, finding_g, searching_h, 0 };
							stk.push(next_GHS);
						}
						// Find a common character in more distance than '_k'.
						else
						{
							// If there is no array indexed by current left border character, then make new array.
							while (wthArrays.size() <= curGHS.id)
							{
								std::vector<SS_GHS> newQue;
								wthArrays.push_back(newQue);
							}

							// Push the common character found in this stage to the array indexed by current left border character.
							SS_GHS wth_GHS{ curGHS.id, finding_g, searching_h, curGHS.s };
							wthArrays[curGHS.id].push_back(wth_GHS);
						}

					}
				}
				// There is no more character in the range of '_str2', and the top of stack does not store sentinel position.
				else if (curGHS.g >= 0)
				{
					mcs = _str1[curGHS.g] + mcs;

					iEnd_SS = _str1.find_last_of(_str1[curGHS.g], iEnd_SS - 1);
					jEnd_SS = _str2.find_last_of(_str1[curGHS.g], jEnd_SS - 1);
				}
			}
		}
		return mcs;
	}

	// Naive algorithm to determine whether the third argument(_mcs) is a mcs of the first one(_str1) and the second one(_str2).
	bool isMCS(std::string _str1, std::string _str2, std::string _mcs)
	{
		int g = -1; int h = -1;
		std::stack<S_GHS> stk;
		S_GHS init{ -1, -1 };
		stk.push(init);

		for (int k = 0; k < _mcs.length(); k++)
		{
			g = _str1.find_first_of(_mcs[k], g + 1);
			h = _str2.find_first_of(_mcs[k], h + 1);
			if (g == -1 || g == _str1.npos || h == -1 || h == _str2.npos)
			{
				return false;
			}

			S_GHS pushObj{ g, h, -1 };
			stk.push(pushObj);
		}

		int iEnd = _str1.length(), jEnd = _str2.length();
		while (stk.size() > 1)
		{
			S_GHS curGH = stk.top();
			stk.pop();

			for (int i = iEnd - 1; i > curGH.g; i--)
			{
				int finding_last_in_Y = _str2.find_last_of(_str1[i], jEnd - 1);
				if (finding_last_in_Y != -1 || finding_last_in_Y != _str2.npos)
				{
					if (finding_last_in_Y > curGH.h)
					{
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
						return false;
					}
				}
			}
			iEnd = _str1.find_last_of(_str1[curGH.g], iEnd - 1);
			jEnd = _str2.find_last_of(_str2[curGH.h], jEnd - 1);
		}
		return true;
	}
}