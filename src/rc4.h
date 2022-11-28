#include <stdio.h>
#include <stdlib.h>

void rc4_keygen(int key_len, unsigned char key[]) {
    for(int i = 0; i < key_len; i++) {
        key[i] = (unsigned char) (256.0*rand()/(RAND_MAX+1.0));
    }
}

void rc4_expandkey(int key_len, unsigned char key[], unsigned char K[]) {
    for(int i = 0; i < 256; i++) {
        K[i] = key[i % key_len];
    }
}

void swap(unsigned char *s, unsigned int i, unsigned int j) {
    unsigned char temp = s[i];
    s[i] = s[j];
    s[j] = temp;
}

void rc4_ksa(unsigned char S[], unsigned char K[]) {
    for (int i = 0; i < 256; i++) {
        S[i] = i;
    }
    for (int j = 0, i = 0; i < 256; i++) {
        j = (j + K[i] + S[i]) & 0xFF;
        swap(S, i, j);
    }
}

unsigned char rc4_prga(int &i, int &j, unsigned char S[]) {
    i = (i + 1) & 0xFF;
    j = (j + S[i]) & 0xFF;
    swap(S, i, j);
    return(S[(S[i] + S[j]) & 0xFF]);
}
