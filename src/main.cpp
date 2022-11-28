#include "rc4.h"
#include "hc256.h"
#include <iostream>
#include <algorithm>
#include <bitset>
#include <set>
#include <vector>
#include <chrono>
#include <random>
#include <thread>
#include <functional>
#include <numeric>
#include <string>

using namespace std;
using namespace chrono;

const int MAX_L = 1000000;
const int MAX_SIZE = 2048;

const double Z_critical_normal_distribution = 3.090;
const double X2_critical_value = 6.6349;

int bitCount[MAX_SIZE + 1][MAX_SIZE + 1];
int rejections[MAX_SIZE];
bitset<MAX_SIZE> D[MAX_L];//input
bitset<MAX_SIZE> output[MAX_L];//output

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

bitset<MAX_SIZE> stream_cipher_RC4(bitset<MAX_SIZE>input,int n,int m){
	unsigned char key[256], K[256], S[256];	
	fill(key, key + 256, 0);
	fill(K, K + 256, 0);
	fill(S, S + 256, 0);	
	unsigned int key_len = n/8;
	int round, keystream_len = m/8;
	for (int i = 0; i < key_len; i++) {
		int pos = (key_len - i - 1)*8;
		bitset<8>temp = bitset<8>(0);
		for (int j = 0; j < 8; j++)
			temp[j] = input[pos + j];
		unsigned char uc = temp.to_ulong();
		key[i] = uc;
	}
	rc4_expandkey(key_len, key, K);
    rc4_ksa(S, K);
	int i = 0, j = 0;
	bitset<MAX_SIZE>result  = bitset<MAX_SIZE>(0);
    for(round = 0; round < keystream_len; round ++) {
        unsigned char v = rc4_prga(i, j, S);
		bitset<8>_byte = bitset<8>(v);
		int pos = (keystream_len - round - 1)*8;
		for (int k = 0; k < 8; k++)
			result[pos + k] = _byte[k];
    }
	return result;
}

bitset<MAX_SIZE> stream_cipher_HC256(bitset<MAX_SIZE>input,int n,int m){//256 bits input key
	const int keystream_len = m/32;int round;
	uint32 key[8],iv[8],DATA[16],P[1024],Q[1024];
	uint32 X[16],Y[16];
	uint32 counter2048; // counter2048 = i mod 2048;

	fill(key, key + 8, 0);
	fill(iv, iv + 8, 0);
	fill(DATA, DATA + keystream_len, 0);

	unsigned int key_len = 8;	// 256 bit key
	for (int it = 0; it < key_len; it++) {
		int pos = (key_len - it - 1)*32;
		bitset<32> temp = bitset<32>(0);
		for (int j = 0; j < 32; j++)
			temp[j] = input[pos + j];		
		key[it] = temp.to_ulong();
	}
    initialization(key,iv,P,Q,counter2048,X,Y);
	encrypt(DATA,counter2048,P,X,Q,Y);

	bitset<MAX_SIZE> result  = bitset<MAX_SIZE>(0);
    for(round = 0; round < keystream_len; round ++) {
		bitset<32>_byte = bitset<32> (DATA[round]);
		int pos = (keystream_len - round - 1)*32;
		for (int k = 0; k < 32; k++)
			result[pos + k] = _byte[k];
    }
	return result;
}

static void parallelFor(unsigned n,function<void(int begin,int end)>func){

    unsigned k = std::thread::hardware_concurrency();
    unsigned n_threads = k == 0 ? 8 : k;
    unsigned chunk_size = n / n_threads;
    unsigned chunk_remainder = n % n_threads;

    std::vector< std::thread > my_threads(n_threads);

	for(unsigned i = 0; i < n_threads; ++i)	{
		int start = i * chunk_size;
		my_threads[i] = std::thread(func, start, start+chunk_size);
	}

    int start = n_threads * chunk_size;
    func( start, start + chunk_remainder);

    std::for_each(my_threads.begin(), my_threads.end(), std::mem_fn(&std::thread::join));
}

int ComputeBIC(int i,bitset<MAX_SIZE>(*cipher)(bitset<MAX_SIZE>,int n,int m),int n,int m,int l){
	int t = 0; double expected = l/2,_d = l/4;
	vector<bitset<MAX_SIZE>> v;
	for (int r = 0; r < l; r++) {
		bitset<MAX_SIZE> X_r_i = D[r];
		X_r_i.flip(i);
		bitset<MAX_SIZE> Y_r_i = cipher(X_r_i,n,m);
		bitset<MAX_SIZE> V_r_i = Y_r_i ^ output[r];
		v.push_back(V_r_i);
	}
	for (int j = 1; j < m; j++) {
		for (int k = 0; k < j; k++) {
			double observed = 0.0;
			for (int ll = 0; ll < l; ll++) observed+=v[ll][j]^v[ll][k];
			double chiSqTest = ((observed - expected) * (observed - expected))/_d;
			if (chiSqTest > X2_critical_value) t++;
		}
	}
	return t;
}

bool BIC_stream_cipher___Parallel(bitset<MAX_SIZE>(*cipher)(bitset<MAX_SIZE>,int n,int m),double a1,int n,int m,int l,bitset<MAX_SIZE>D[]){
	double T = 0;
	parallelFor(n, [&](int start, int end){
	    for(int i = start; i < end; ++i)
	        rejections[i] = ComputeBIC(i,cipher,n,m,l);
	});
	T = accumulate(rejections, rejections + n, T);
	cout << "Numero de rechazos T = " << T << endl;
	double C_m_2 = (m * (m - 1))/2, _n = n, _m = m;
	cout << "valor alpha_1 = " << a1 << endl;
	double E_T_H0 =  a1 * _n * C_m_2;
	cout << "valor esperado E(T|H0) = " << E_T_H0 << endl;
	double variance = E_T_H0 * (1 - a1);
	cout << "varianza = " << variance << endl;
	double Z_T = (T - E_T_H0) / sqrt(variance);
	cout << "valor Z_T = " << Z_T << endl;
	return Z_T <= Z_critical_normal_distribution;
}

bool BIC_stream_cipher_Sequential(bitset<MAX_SIZE>(*cipher)(bitset<MAX_SIZE>,int n,int m),double a1,int n,int m,int l,bitset<MAX_SIZE>D[]){
	double T = 0,expected=l/2,_d=l/4;
	for (int i = 0; i < n; i++) {
		vector<bitset<MAX_SIZE>> v;
		for (int r = 0; r < l; r++) {
			bitset<MAX_SIZE> X_r_i = D[r];
			X_r_i.flip(i);
			bitset<MAX_SIZE> Y_r_i = cipher(X_r_i,n,m);
			bitset<MAX_SIZE> V_r_i = Y_r_i ^ output[r];
			v.push_back(V_r_i);
		}
		for (int j = 1; j < m; j++) {
			for (int k = 0; k < j; k++) {
				double observed = 0.0;
				for (int ll = 0; ll < l; ll++) observed+=v[ll][j]^v[ll][k];
				double chiSqTest = ((observed - expected)*(observed - expected))/_d;
				if (chiSqTest > X2_critical_value) T++;
			}
		}
	}
	double C_m_2 = (m * (m - 1))/2, _n = n, _m = m;
	double E_T_H0 =  a1 * _n * C_m_2;
	double variance = E_T_H0 * (1 - a1);
	double Z_T = (T - E_T_H0) / sqrt(variance);
	return Z_T <= Z_critical_normal_distribution;
}

int ComputeSAC(int i,bitset<MAX_SIZE>(*cipher)(bitset<MAX_SIZE>,int n,int m),int n,int m,int l){
	int t = 0;double expected=l/2;
	for (int r = 0; r < l; r++) {
		bitset<MAX_SIZE> X_r_i = D[r];
		X_r_i.flip(i);
		bitset<MAX_SIZE> Y_r_i = cipher(X_r_i,n,m);
		bitset<MAX_SIZE> V_r_i = Y_r_i ^ output[r];
		for (int k = 0; k < m; k++)
			bitCount[i][k] += V_r_i.test(k);
	}
	for (int j = 0; j < m; j++) {
		double observed = bitCount[i][j];	//cant 1s en la columna j
		double chiSqTest = ((observed - expected) * (observed - expected)) / expected;
		if (chiSqTest > X2_critical_value) t++;
	}
	return t;
}

bool SAC_stream_cipher___Parallel(bitset<MAX_SIZE>(*cipher)(bitset<MAX_SIZE>,int n,int m),double a1,int n,int m,int l,bitset<MAX_SIZE>D[]){
	double T = 0;
	parallelFor(n, [&](int start, int end){ 
	    for(int i = start; i < end; ++i)
	        rejections[i] = ComputeSAC(i,cipher,n,m,l);
	});
	T = accumulate(rejections, rejections + n, T);
	cout << "Numero de rechazos T = " << T << endl;
	double _n = n, _m = m;
	cout << "valor alpha_1 = " << a1 << endl;
	double E_T_H0 =  a1 * _n * _m;
	cout << "valor esperado E(T|H0) = " << E_T_H0 << endl;
	double variance = E_T_H0 * (1 - a1);
	cout << "varianza = " << variance << endl;
	double Z_T = (T - E_T_H0) / sqrt(variance);
	cout << "valor Z_T = " << Z_T << endl;
	return Z_T <= Z_critical_normal_distribution;
}

bool SAC_stream_cipher_Sequential(bitset<MAX_SIZE>(*cipher)(bitset<MAX_SIZE>,int n,int m),double a1,int n,int m,int l,bitset<MAX_SIZE>D[]){
	double T = 0, expected=l/2;
	for (int i = 0; i < n; i++) {
		for (int r = 0; r < l; r++) {
			bitset<MAX_SIZE> X_r_i = D[r];
			X_r_i.flip(i);
			bitset<MAX_SIZE> Y_r_i = cipher(X_r_i,n,m);
			bitset<MAX_SIZE> V_r_i = Y_r_i ^ output[r];
			for (int k = 0; k < m; k++)
				bitCount[i][k] += V_r_i.test(k);
		}
		for (int j = 0; j < m; j++) {
			double observed = bitCount[i][j];	//cant 1s en la columna j
			double chiSqTest = ((observed - expected) * (observed - expected)) / expected;
			if (chiSqTest > X2_critical_value) T++;
		}
	}
	double _n = n, _m = m;
	double E_T_H0 =  a1 * _n * _m;
	double variance = E_T_H0 * (1 - a1);
	double Z_T = (T - E_T_H0) / sqrt(variance);
	return Z_T <= Z_critical_normal_distribution;
}

void make_data_test(int n,int l){	
	set<string> vectors; int ll = 0;
	while (vectors.size() < l){
		string s = "";
		for (int j = 0; j < n; j++){
			int rd = rng();
			if (rd % 2) s+="1";
			else s+="0";
		}
		vectors.insert(s);
	}
	for (auto &&v : vectors)D[ll++]=bitset<MAX_SIZE>(v);
}

void pre_compute(bitset<MAX_SIZE>(*cipher)(bitset<MAX_SIZE>,int n,int m),int n,int m,int l){
	make_data_test(n,l);//build D
	for (int i = 0; i < l; i++)output[i]=cipher(D[i],n,m);
}

int main() {
	
	int n, m, l, alg_option, cipher_opt;

	cout << "Seleccione la opcion deseada" << endl;
	cout << "1 Ejecutar algoritmo BIC" << endl;
	cout << "2 Ejecutar algoritmo SAC" << endl;
	cin >> alg_option;

	cout << "Seleccione la opcion deseada" << endl;
	cout << "1 Evaluar cifrador RC4" << endl;
	cout << "2 Evaluar cifrador HC256" << endl;
	cin >> cipher_opt;

	if(cipher_opt==1){
		cout << "Entrar parametro n en intervalo [8,2048] divisible por 8" << endl;
		cin >> n;
		if(n<8 || n>2048 || n%8!=0){
			cout << "Entrada incorrecta" <<endl;
			return 0;
		}
		cout << "Entrar parametro m divisible por 8" << endl;
		cin >> m;
		if(m > MAX_SIZE){
			cout << "Maximo valor m excedido" <<endl;
			return 0;
		}
		if(m%8!=0){
			cout << "Entrada incorrecta" <<endl;
			return 0;
		}
		cout << "Entrar parametro L" << endl;
		cin >> l;
		if(l > MAX_L){
			cout << "Maximo valor L excedido" <<endl;
			return 0;
		}
	}else if (cipher_opt==2)
	{
		n = 256;
		cout << "Entrar parametro m en intervalo [32,512] divisible por 32" << endl;
		cin >> m;
		if(m<32 || m>512 || m%32!=0){
			cout << "Entrada incorrecta" <<endl;
			return 0;
		}
		cout << "Entrar parametro L" << endl;
		cin >> l;
		if(l > MAX_L){
			cout << "Maximo valor L excedido" <<endl;
			return 0;
		}
	}else
	{
		cout << "Entrada incorrecta" << endl;
	}	

	int pow = 1;
	for (int exp = 1; exp <= n; exp++) {
		pow*=2;
		if(pow>=l) break;
	}
	if(pow < l) {
		cout << "2^n = "<< pow << " < L" << endl;
		cout << "Entrada incorrecta" <<endl;
		return 0;
	}

	cout << "parametros n=" << n << " m=" << m << " L=" << l << endl;

	if(cipher_opt==1){//Build D[r], output[r] sets
		pre_compute(&stream_cipher_RC4,n,m,l);
	} else pre_compute(&stream_cipher_HC256,n,m,l);	
	cout << "entrada Xr, Yr computada" << endl;

	if(alg_option==1){//BIC
		bool TEST_RESULT;
		auto start = high_resolution_clock::now();
		if(cipher_opt==1) TEST_RESULT = BIC_stream_cipher___Parallel(&stream_cipher_RC4,0.01,n,m,l,D);
		else TEST_RESULT = BIC_stream_cipher___Parallel(&stream_cipher_HC256,0.01,n,m,l,D);
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<seconds>(stop - start);
		cout << "Tiempo de computo en segundos: " << duration.count() << " seconds" << "\n";
		auto duration2 = duration_cast<minutes>(stop - start);
		cout << "Tiempo de computo en minutos: " << duration2.count() << " minutos" << "\n";
		if(TEST_RESULT)	cout << "El cifrador satisface el BIC" << endl;
		else cout << "El cifrador no satisface el BIC" << endl;
	}else {
		bool TEST_RESULT;
		auto start = high_resolution_clock::now();		
		if(cipher_opt==1) TEST_RESULT = SAC_stream_cipher___Parallel(&stream_cipher_RC4,0.01,n,m,l,D);
		else TEST_RESULT = SAC_stream_cipher___Parallel(&stream_cipher_HC256,0.01,n,m,l,D);
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<seconds>(stop - start);
		cout << "Tiempo de computo en segundos: " << duration.count() << " seconds" << "\n";
		auto duration2 = duration_cast<minutes>(stop - start);
		cout << "Tiempo de computo en minutos: " << duration2.count() << " minutos" << "\n";
		if(TEST_RESULT)	cout << "El cifrador satisface el SAC" << endl;
		else cout << "El cifrador no satisface el SAC" << endl;
	}
	return 0;
}
