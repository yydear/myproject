#include<iostream>
using namespace std;

#include<iostream>
#include<ctime>
using namespace std;


#include<cmath>
#include<stdlib.h>
#include<stdio.h>


//-------------------------------------------------------------------MixColumn------------------------------------------------------------------------------------------------------//
unsigned int Xto2X(unsigned int* x)
{
	(*x) = (*x) << 1;//x左移1位
	if ((*x) <= 127)
	{
		return (*x);//如果b7=0，即x小于127，直接返回常数
	}
	else
	{
		return (*x ^ 0x100) ^ 0x1b;
	}
}

void Matrix_device(const unsigned int in[4][4], unsigned int out[4][4])
{
	for (int i = 0; i < 4; i += 1)
	{
		for (int j = 0; j < 4; j += 1)
		{
			out[i][j] = in[j][i];
		}
	}
}

//在[min,max]范围内获取一个随机整数
int getRand(int min, int max) {
	return (rand() % (max - min + 1)) + min;
}

void MixColumns(const unsigned int in[4][4], unsigned int out[4][4])
{
	unsigned int temp1[4][4];
	unsigned int temp2[4][4];
	//先将输入矩阵进行转置（取连续地址）
	Matrix_device(in, temp1);


	for (int i = 0; i < 4; i += 1)
	{
		for (int j = 0; j < 4; j += 1)
		{
			unsigned int x(temp1[i][j] ^ temp1[i][(j + 1) % 4]);
			unsigned int x2x = Xto2X(&x);
			temp2[i][j] = x2x ^ temp1[i][(j + 1) % 4] ^ temp1[i][(j + 2) % 4] ^ temp1[i][(j + 3) % 4];
		}
	}

	//先将输入矩阵进行转置（取连续地址）
	Matrix_device(temp2, out);
}
//-------------------------------------------------------------------MixColumn------------------------------------------------------------------------------------------------------//

//-------------------------------------------------------------------GetKey------------------------------------------------------------------------------------------------------//
void SR(unsigned int mat[4], unsigned int sr[4])//ShiftRow 
{
	for (int i = 1; i < 4; i += 1) {
		sr[i] = mat[i - 1];
	}
	sr[0] = mat[3];
}
void BRO(unsigned int bro[8], unsigned int n) //get Binary reverse order
{
	unsigned int m = n;
	for (int i = 0; i < 8; i += 1) {
		bro[i] = m % 2;
		m /= 2;
	}
}
unsigned int RTO(unsigned int bro[8])//Return to the Original Number
{
	int sum = 0;
	for (int i = 0; i < 8; i += 1) {
		sum += pow(2, i) * bro[i];
	}
	return sum;
}
unsigned int GetInverse(unsigned int n)//通过查表法求模逆
{
	int table[256];//原根指数正表
	table[0] = 1;
	for (int i = 1; i < 256; i += 1) {
		table[i] = (table[i - 1] << 1) ^ table[i - 1];
		if (table[i] & 0x100) {
			table[i] ^= 0x11B;
		}
	}
	int arc_table[256];//离散对数表
	for (int i = 0; i < 256; ++i) {
		arc_table[table[i]] = i;
	}
	int inverse[256];
	inverse[0] = 0;
	for (int i = 1; i < 256; ++i) {
		int j = 255 - arc_table[i];//费马小定理
		inverse[i] = table[j];
	}
	return inverse[n];
}
unsigned int GetS(unsigned int n)//STable lookup
{
	unsigned int t = GetInverse(n);
	unsigned int bro[8];
	BRO(bro, t);
	unsigned int p[8];
	unsigned int a[8][8] = {
		{1,0,0,0,1,1,1,1},
		{1,1,0,0,0,1,1,1},
		{1,1,1,0,0,0,1,1},
		{1,1,1,1,0,0,0,1},
		{1,1,1,1,1,0,0,0},
		{0,1,1,1,1,1,0,0},
		{0,0,1,1,1,1,1,0},
		{0,0,0,1,1,1,1,1},
	};
	unsigned int b[8] = { 1,1,0,0,0,1,1,0 };
	for (int i = 0; i < 8; i++) {
		int sum = 0;
		for (int j = 0; j < 8; j++) {
			sum ^= a[i][j] * bro[j];
		}
		p[i] = sum ^ b[i];
	}
	unsigned int m = RTO(p);
	return m;
}
void TSP(unsigned int p[4][4], unsigned int q[4][4])//Transportation 后是前的转置
{
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			q[j][i] = p[i][j];
		}
	}
}
void T(unsigned int LastKey[4][4], unsigned int thk[4][4], int n)//生成轮密钥第一列的过程函数 
{
	unsigned int Rcon[10][4] = {
		{0x01,0x00,0x00,0x00},
		{0x02,0x00,0x00,0x00},
		{0x04,0x00,0x00,0x00},
		{0x08,0x00,0x00,0x00},
		{0x10,0x00,0x00,0x00},
		{0x20,0x00,0x00,0x00},
		{0x40,0x00,0x00,0x00},
		{0x80,0x00,0x00,0x00},
		{0x1b,0x00,0x00,0x00},
		{0x36,0x00,0x00,0x00},
	};
	unsigned int lsk[4][4];
	TSP(LastKey, lsk);//转置，把列处理转换为行处理
	unsigned int l3[4];
	SR(lsk[3], l3);//对密钥的最后一列进行移位
	for (int i = 0; i < 4; i++) {
		thk[0][i] = GetS(l3[i]) ^ Rcon[n - 1][i] ^ lsk[0][i];//移位后的列先进行S盒置换，后与密钥第一列异或，后放入预处理密钥的第一列
	}
}
void GetKey(unsigned int LastKey[4][4], unsigned int KeyNext[4][4], int n) {
	unsigned int thk[4][4];//转置前的密钥的承接
	T(LastKey, thk, n);
	for (int i = 1; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			thk[i][j] = thk[i - 1][j] ^ LastKey[j][i];
		}
	}
	TSP(thk, KeyNext);
}
//-------------------------------------------------------------------GetKey------------------------------------------------------------------------------------------------------//




void ShiftRows(const unsigned int in[4][4], unsigned int out[4][4])
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			//每行位移量和i有关，向左移动i位，故减i
			int temp = (j - i + 4) % 4;
			out[i][temp] = in[i][j];
		}
	}
}
/*两种不同实现的SubBytes可用来比较时间*/
void SubBytes1(unsigned int in[4][4], unsigned int out[4][4])
{
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			out[i][j] = GetS(in[i][j]);
		}
	}
}
void SubBytes2(unsigned int ST[16][16], unsigned int in[4][4], unsigned int out[4][4])
{
	for (int i = 0; i < 16; i++) {
		for (int j = 0; j < 16; j++) {
			int a = in[i][j] >> 4;
			int b = in[i][j] & 0x0f;
			out[i][j] = ST[a][b];
		}
	}
}
/*两种不同实现的SubBytes可用来比较时间*/
void GetST(unsigned int ST[16][16])//一键生成S表
{
	for (int i = 0; i < 16; i++) {
		int h = i << 4;
		for (int j = 0; j < 16; j++) {
			ST[i][j] = GetS(h ^ j);
		}
	}
}
void AddRoundKey(unsigned int p[4][4], unsigned int k[4][4]) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			p[i][j] ^= k[i][j];
		}
	}
}
void Equ(unsigned int LastKey[4][4], unsigned int lsk[4][4])//equal 后和前相等
{
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			lsk[i][j] = LastKey[i][j];
		}
	}
}


//===============================================//AES的综合实现===================================================
void AES(const unsigned int p[4][4], unsigned int k[4][4], unsigned int Out[4][4])
{
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			Out[i][j] = p[i][j];
		}
	}
	unsigned int st[16][16];
	GetST(st);
	AddRoundKey(Out, k);//白化
	unsigned int re[4][4];
	unsigned int lk[4][4];//更新存放上轮密钥
	unsigned int nk[4][4];//更新存放下轮密钥
	Equ(k, lk);
	for (int i = 1; i < 10; i++)//前九轮加密
	{
		SubBytes1(Out, re);
		//SubBytes2(st, Out, re);
		ShiftRows(re, Out);
		MixColumns(Out, re);
		GetKey(lk, nk, i);
		AddRoundKey(re, nk);
		Equ(re, Out);
		Equ(nk, lk);
	}
	SubBytes1(Out, re);
	//SubBytes2(st, Out, re);
	ShiftRows(re, Out);
	GetKey(lk, nk, 10);
	AddRoundKey(re, nk);
}


//字符串编码函数：输入一个最大长度为16的字符串，通过ASCLL码转换，
//输出一个4*4矩阵，字符串不足16长则转换后自动补0
void String_To_NumberMat(const char in[], unsigned int out[4][4]) {
	int len = strlen(in);
	for (int i = 0; i < 4; i += 1) {
		for (int j = 0; j < 4; j += 1) {
			int x = i * 4 + j;
			if (x < len) {
				out[i][j] = (unsigned int)in[x];
			}
			else {
				out[i][j] = 0;
			}

		}
	}
}


int main(void)
{

	{
		char st[] = "ILIKEPLAYCOMPUTERGAME";//在此处放入明文
		unsigned int Pt_in[4][4];
		unsigned int Key[4][4];
		unsigned int Ct_out[4][4];

		String_To_NumberMat(st, Pt_in);//将st编译成4*4整数矩阵
		String_To_NumberMat(st, Key);

		printf("明文=\n");
		for (int i = 0; i < 4; i += 1) {
			for (int j = 0; j < 4; j += 1) {
				printf("%x", Pt_in[i][j]);
				printf("\t");
			}
			printf("\n");
		}

		printf("密钥=\n");
		for (int i = 0; i < 4; i += 1) {
			for (int j = 0; j < 4; j += 1) {
				printf("%x", Key[i][j]);
				printf("\t");
			}
			printf("\n");
		}

		AES(Pt_in, Key, Ct_out);

		printf("密文=\n");
		for (int i = 0; i < 4; i += 1) {
			for (int j = 0; j < 4; j += 1) {
				printf("%x", Ct_out[i][j]);
				printf("\t");
			}
			printf("\n");
		}

	}
}