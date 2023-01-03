# A New Algorithm of Finding a Maximal Common Subsequence of Two Strings

-   Hyeonjun Shin, & Jeon Seop Sim. A New Algorithm of Finding a Maximal Common Subsequence of Two Strings. Korea Software Congress 2022, 1212-1214.
- This paper is submitted at Korea Software Congress 2022 of Korean Institute of Information Scientists and Engineers. And, this paper won the Top Theorical Computer Science Paper Award.


## Maximal Common Subsequences
- MCS(Maximal Common Subsequences) is common subsequences that is no longer a common subsequence when any character is inserted.
- This was proposed by Campbell B. Fraser and Robert W. Irving in 1975[[1]](https://www.sciencedirect.com/science/article/pii/S0890540196900115).
- Yoshifumi Sakai proposed a (sub)linearithmic-time, linear-space algorithm for finding a maximal common subsequence of two strings and a linear-time algorithm for determining if a common subsequence of two strings is maximal[[2]](https://www.sciencedirect.com/science/article/pii/S0304397519304074).
- DongYeop Lee and Joon Chae Na proposed a algorithm finding a longer MCS than Y. Sakai's[[3]](https://www.dbpia.co.kr/pdf/pdfView.do?nodeId=NODE11035975&googleIPSandBox=false&mark=0&useDate=&ipRange=false&accessgl=Y&language=ko_KR&hasTopBanner=true)[[4]](https://www.dbpia.co.kr/pdf/pdfView.do?nodeId=NODE11100316&googleIPSandBox=false&mark=0&useDate=&ipRange=false&accessgl=Y&language=ko_KR&hasTopBanner=true).

## A New Algorithm
- I and my professor proposed a new algorithm finding a longer MCS than Y.Sakai's and Lee & Na's.
- This algorithm check _k_ common characters and choose one that make next searching range more longer.
- Its time complexity and space complexity are equal to Y. Sakai's.





---
[[1]](https://www.sciencedirect.com/science/article/pii/S0890540196900115) Fraser, Campbell B., Robert W. Irving, and Martin Middendorf. "Maximal common subsequences and minimal common supersequences." _information and computation_ 124.2 (1996): 145-153.

[[2]](https://www.sciencedirect.com/science/article/pii/S0304397519304074)Sakai, Yoshifumi. "Maximal common subsequence algorithms." _Theoretical Computer Science_ 793 (2019): 132-139.
[[3]](https://www.dbpia.co.kr/pdf/pdfView.do?nodeId=NODE11035975&googleIPSandBox=false&mark=0&useDate=&ipRange=false&accessgl=Y&language=ko_KR&hasTopBanner=true)이동엽, and 나중채. "극대 공통 부분 서열 알고리즘을 개선하기 위한 간단한 전략." _한국정보과학회 학술발표논문집_ (2021): 1149-1151.
[[4]](https://www.dbpia.co.kr/pdf/pdfView.do?nodeId=NODE11100316&googleIPSandBox=false&mark=0&useDate=&ipRange=false&accessgl=Y&language=ko_KR&hasTopBanner=true)이동엽, and 나중채. "더 긴 극대 공통 부분 서열을 찾기 위한 알고리즘." _정보과학회논문지_ 49.7 (2022): 507-513.
