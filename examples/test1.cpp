
Vector_t<double> tank6Eq(Vector_t<double> x, Vector_t<double> u){
    // parameter
    k = 0.015;
    k2 = 0.01;
    g = 9.81; 

    // differential equation
    Vector_t<double> result;
    result[0]=u[0]+0.1+k2*(4-x[5])-k*sqrt(2*g)*sqrt(x[0]);
    result[1]=k*sqrt(2*g)*(sqrt(x[0])-sqrt(x[1]));
    result[2]=k*sqrt(2*g)*(sqrt(x[1])-sqrt(x[2]));
    result[3]=k*sqrt(2*g)*(sqrt(x[2])-sqrt(x[3]));
    result[4]=k*sqrt(2*g)*(sqrt(x[3])-sqrt(x[4]));
    result[5]=k*sqrt(2*g)*(sqrt(x[4])-sqrt(x[5]));
    return result;
}

double (*func_type)(Vector_t<double> x, Vector_t<double> u);
// parameter
k = 0.015;
k2 = 0.01;
g = 9.81; 
double tank6Eq_f0(Vector_t<double> x, Vector_t<double> u){
    return u[0]+0.1+k2*(4-x[5])-k*sqrt(2*g)*sqrt(x[0]);
}

double tank6Eq_f1(Vector_t<double> x, Vector_t<double> u){
    return k*sqrt(2*g)*(sqrt(x[0])-sqrt(x[1]));
}

double tank6Eq_f2(Vector_t<double> x, Vector_t<double> u){
    return k*sqrt(2*g)*(sqrt(x[1])-sqrt(x[2]));
}

double tank6Eq_f3(Vector_t<double> x, Vector_t<double> u){
    return k*sqrt(2*g)*(sqrt(x[2])-sqrt(x[3]));
}

double tank6Eq_f4(Vector_t<double> x, Vector_t<double> u){
    return k*sqrt(2*g)*(sqrt(x[3])-sqrt(x[4]));
}

double tank6Eq_f5(Vector_t<double> x, Vector_t<double> u){
    return k*sqrt(2*g)*(sqrt(x[4])-sqrt(x[5]));
}

int main(){
    func_type tank6Eq_f[6];
    // func_type tank6Eq;
}