int somma(int a, int b){
    return a+b;
}

void somma_ref(int &a, int&b){
    a+=b; //modifica il valore di a
    return;
}

int somma_ref_const(const int& a, const int &b){
    return a+b;
}