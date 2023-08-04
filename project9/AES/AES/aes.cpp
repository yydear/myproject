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
	(*x) = (*x) << 1;//x����1λ
	if ((*x) <= 127)
	{
		return (*x);//���b7=0����xС��127��ֱ�ӷ��س���
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

//��[min,max]��Χ�ڻ�ȡһ���������
int getRand(int min, int max) {
	return (rand() % (max - min + 1)) + min;
}

void MixColumns(const unsigned int in[4][4], unsigned int out[4][4])
{
	unsigned int temp1[4][4];
	unsigned int temp2[4][4];
	//�Ƚ�����������ת�ã�ȡ������ַ��
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

	//�Ƚ�����������ת�ã�ȡ������ַ��
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
unsigned int GetInverse(unsigned int n)//ͨ�������ģ��
{
	int table[256];//ԭ��ָ������
	table[0] = 1;
	for (int i = 1; i < 256; i += 1) {
		table[i] = (table[i - 1] << 1) ^ table[i - 1];
		if (table[i] & 0x100) {
			table[i] ^= 0x11B;
		}
	}
	int arc_table[256];//��ɢ������
	for (int i = 0; i < 256; ++i) {
		arc_table[table[i]] = i;
	}
	int inverse[256];
	inverse[0] = 0;
	for (int i = 1; i < 256; ++i) {
		int j = 255 - arc_table[i];//����С����
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
void TSP(unsigned int p[4][4], unsigned int q[4][4])//Transportation ����ǰ��ת��
{
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			q[j][i] = p[i][j];
		}
	}
}
void T(unsigned int LastKey[4][4], unsigned int thk[4][4], int n)//��������Կ��һ�еĹ��̺��� 
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
	TSP(LastKey, lsk);//ת�ã����д���ת��Ϊ�д���
	unsigned int l3[4];
	SR(lsk[3], l3);//����Կ�����һ�н�����λ
	for (int i = 0; i < 4; i++) {
		thk[0][i] = GetS(l3[i]) ^ Rcon[n - 1][i] ^ lsk[0][i];//��λ������Ƚ���S���û���������Կ��һ����򣬺����Ԥ������Կ�ĵ�һ��
	}
}
void GetKey(unsigned int LastKey[4][4], unsigned int KeyNext[4][4], int n) {
	unsigned int thk[4][4];//ת��ǰ����Կ�ĳн�
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
			//ÿ��λ������i�йأ������ƶ�iλ���ʼ�i
			int temp = (j - i + 4) % 4;
			out[i][temp] = in[i][j];
		}
	}
}
/*���ֲ�ͬʵ�ֵ�SubBytes�������Ƚ�ʱ��*/
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
/*���ֲ�ͬʵ�ֵ�SubBytes�������Ƚ�ʱ��*/
void GetST(unsigned int ST[16][16])//һ������S��
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
void Equ(unsigned int LastKey[4][4], unsigned int lsk[4][4])//equal ���ǰ���
{
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			lsk[i][j] = LastKey[i][j];
		}
	}
}


//===============================================//AES���ۺ�ʵ��===================================================
void AES(const unsigned int p[4][4], unsigned int k[4][4], unsigned int Out[4][4])
{
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			Out[i][j] = p[i][j];
		}
	}
	unsigned int st[16][16];
	GetST(st);
	AddRoundKey(Out, k);//�׻�
	unsigned int re[4][4];
	unsigned int lk[4][4];//���´��������Կ
	unsigned int nk[4][4];//���´��������Կ
	Equ(k, lk);
	for (int i = 1; i < 10; i++)//ǰ���ּ���
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


//�ַ������뺯��������һ����󳤶�Ϊ16���ַ�����ͨ��ASCLL��ת����
//���һ��4*4�����ַ�������16����ת�����Զ���0
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
		char st[] = "ILIKEPLAYCOMPUTERGAME";//�ڴ˴���������
		unsigned int Pt_in[4][4];
		unsigned int Key[4][4];
		unsigned int Ct_out[4][4];

		String_To_NumberMat(st, Pt_in);//��st�����4*4��������
		String_To_NumberMat(st, Key);

		printf("����=\n");
		for (int i = 0; i < 4; i += 1) {
			for (int j = 0; j < 4; j += 1) {
				printf("%x", Pt_in[i][j]);
				printf("\t");
			}
			printf("\n");
		}

		printf("��Կ=\n");
		for (int i = 0; i < 4; i += 1) {
			for (int j = 0; j < 4; j += 1) {
				printf("%x", Key[i][j]);
				printf("\t");
			}
			printf("\n");
		}

		AES(Pt_in, Key, Ct_out);

		printf("����=\n");
		for (int i = 0; i < 4; i += 1) {
			for (int j = 0; j < 4; j += 1) {
				printf("%x", Ct_out[i][j]);
				printf("\t");
			}
			printf("\n");
		}

	}
}