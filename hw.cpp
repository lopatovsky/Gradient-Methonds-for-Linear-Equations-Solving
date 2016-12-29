#include <bits/stdc++.h>
#include <unistd.h>

using namespace std;
typedef pair<int,int> ii;
typedef vector<int> vi;
typedef vector<ii> vii;
#define MP make_pair
#define PB push_back
#define ff first
#define ss second
#define TR(it,c) for( typeof(c.begin()) it = c.begin(); it != c.end(); ++it )
#define TRR(it,c) for( typeof(c.rbegin()) it = c.rbegin(); it != c.rend(); ++it
#define REP(i,a) for (int i = 0; i < (a); i++)
#define FOR(i,a,b) for (int i = (a); i <= (b); i++)
#define FORD(i,a,b) for (int i = (a); i >= (b); i--)

#define DRI(a) int a; scanf("%d", &a);
#define DRII(a, b) int a, b; scanf("%d %d", &a, &b);
#define DRIII(a, b, c) int a, b, c; scanf("%d %d %d", &a, &b, &c);
#define RI(a) scanf("%d", &a);
#define RII(a, b) scanf("%d %d", &a, &b);
#define RIII(a, b, c) scanf("%d %d %d", &a, &b, &c);
#define MM(arr, num) memset((arr), (num), sizeof((arr)))
#define DEB(x) cerr << ">>> " << (#x) << " -> " << (x) << endl;
#define DEBA(x,n) cerr << (#x) << " "; deba((x),(n));
void deba( double * a, int n){ cerr << "| "; REP(i,n) cerr << a[i] << " "; cerr << "|" << endl;}


inline bool EQ(double a, double b) { return fabs(a-b) < 1e-9; }

const int INF = 1<<30;
typedef long long ll;
typedef unsigned long long ull;
//#define double long double
/*******************************************************/

const double EPS = 0.5;

const int M = 2600000;
const int N = 200000;

//const char FILE_A[] = "mat_cct_7650.txt";
//const char FILE_B[] = "prava_str_cct_7650.txt";

const char FILE_A[] = "mat_ps_8200.txt";
const char FILE_B[] = "prava_str_ps_8200.txt";

//const char FILE_A[] = "smallA.txt";
//const char FILE_B[] = "smallB.txt";

//const char FILE_A[] = "A.txt";
//const char FILE_B[] = "B.txt";


int n,m;


int row[M];
int col[M];
double b[N];
double sol[N];
double val[M];
//double test[] = {0,1,0};

double x[2][N]; //todo try [N][2]
double r[2][N];
double s[2][N];

double ark[N];

const int N2D = 9000;
double mat[N2D][N2D];

void fill_in_2D(){
    if(n > N2D) { cout << "To big matrix" << endl; return; }
    REP(i,n)REP(j,n) mat[i][j] = 0;
    REP(i, m) { mat[row[i]][col[i]] = mat[col[i]][row[i]] = val[i]; }
}

void gauss(){

    fill_in_2D();

    FOR(i, 1, n-1){
        FOR( j, 0, i-1){
            double a = - mat[j][j];
            double d =   mat[i][j];
            FOR(k, j, n-1) mat[i][k] = a * mat[i][k] + d * mat[j][k];
            b[i] = a * b[i] + d * b[j];
        }
    }

    REP(i,n){ REP(j,n) cout << mat[i][j] << " "; cout << endl; }
    REP(i,n) cout << b[i] << "-"; cout << endl;

    FORD(i,n-1,0){
        double v = b[i];
        FORD(j,n-1,i+1) v -= sol[j] * mat[i][j];
        sol[i] = v / mat[i][i];
    }

    REP(i,n) cout << sol[i] << " "; cout << endl;


}


inline void get_random_vec( double * r){
    REP(i,n) r[i] = rand() / (double)RAND_MAX;
    //DEBA(r,n);
}

inline void vec_sub( double * r, double * a, double * b ){
    REP(i,n) r[i] = a[i] - b[i];
}
inline void vec_add( double * r, double * a, double * b ){
    REP(i,n) r[i] = a[i] + b[i];
}

// a + c * b  -> c is const
inline void vec_add_mul( double * r, double * a, double * b, double alpha_b ){
    REP(i,n) r[i] = a[i] + alpha_b * b[i];
}
// a - c * b  -> c is const
inline void vec_sub_mul( double * r, double * a, double * b, double alpha_b ){
    REP(i,n) r[i] = a[i] - alpha_b * b[i];
}

inline void mat_vec_mul(  double * r, double * b){
    REP(i,n){
        r[i] = 0;
        REP(j,n) r[i] += mat[i][j] * b[j];
    }
}
inline double vec_mul( double * a, double * b ){
    double ret = 0;
    REP(i,n) ret += a[i] * b[i];
    return ret;
}

inline double vec_size(double * a){
    double ret = 0;
    REP(i,n) ret+= a[i]*a[i];
    return sqrt(ret);
}

inline double copy_vec( double * to, double * from){
    REP(i,n) to[i] = from[i];
}




void descent_2D(){

    fill_in_2D();
    get_random_vec( x[0] );
    //x[0][0] = 0;
    //x[0][1] = 1;

    mat_vec_mul( r[0] , x[0] ); //r[0] is Ax0 here.

    //DEBA(r[0],n);
   // DEBA(b,n);

    vec_sub(  r[0],  b,  r[0] );

    //DEBA(r[0],n);

    int k = 0, k1 = 1;
    double alpha, rk_rk, rk_ark;

    while(1){

        mat_vec_mul( ark, r[k] );

        //DEBA(ark,n)
        //DEBA(r[k],n)

        rk_rk = vec_mul( r[k], r[k] );
        rk_ark = vec_mul( r[k], ark );
        if(rk_ark == 0){ cout << "Zero division" << endl;  exit(0);}
        alpha =  rk_rk / rk_ark;


        vec_add_mul( x[k1], x[k], r[k], alpha );

        //DEBA(x[k1],n)

        vec_sub_mul( r[k1], r[k], ark,  alpha );

        DEB( vec_size(r[k1] ) );

        /*usleep(100000);
        vec_sub( s[0], test, x[k1]);
        DEB( vec_size(s[0]));
        */

        if( vec_size(r[k1] ) < EPS ) { break; }



        k  ^= 1;
        k1 ^= 1;
    }


    DEBA( x[k1] ,n );
    mat_vec_mul( x[k], x[k1] );
    vec_sub( r[0], b, x[k] );

    cout << "error: " <<  vec_size(r[0]) << endl;

}


void gradient_2D(){

    fill_in_2D();
    get_random_vec( x[0] );

    mat_vec_mul( r[0] , x[0] ); //r[0] is Ax0 here.
    vec_sub(  r[0],  b,  r[0] );
    copy_vec( s[0], r[0] );

    //DEBA(r[0],n);

    int k = 0, k1 = 1;
    double alpha, beta, rk1_rk1, rk_rk, sk_ask;

    while(1){

        mat_vec_mul( ark, s[k] );  // ark is here ask

        rk_rk = vec_mul( r[k], r[k] );
        sk_ask = vec_mul( s[k], ark );
        if(sk_ask == 0){ cout << "Zero division" << endl;  break;}
        alpha =  rk_rk / sk_ask;

        vec_add_mul( x[k1], x[k], s[k], alpha );
        vec_sub_mul( r[k1], r[k], ark,  alpha );


        DEB( vec_size(r[k1] ) );

        if( vec_size(r[k1] ) < EPS ) { break; }
        /*new part*/

        rk1_rk1 = vec_mul( r[k1], r[k1] );
        if(rk_rk == 0){ cout << "Zero division" << endl;  break;}
        beta = rk1_rk1 / rk_rk;

        DEB(beta)

        vec_add_mul( s[k1], r[k1], s[k], beta );
        /* */

        k  ^= 1;
        k1 ^= 1;
    }


    DEBA( x[k1] ,n );
    mat_vec_mul( x[k], x[k1] );
    vec_sub( r[0], b, x[k] );

    cout << "error: " <<  vec_size(r[0]) << endl;


}



int main() {

  srand( time(NULL) );

  FILE * fa = fopen(FILE_A,"r+");
  FILE * fb = fopen(FILE_B,"r+");

  if(fa == NULL ) cout << "wrong file a" << endl;
  if(fb == NULL ) cout << "wrong file b" << endl;

  fscanf(fa,"%d%d", &n,&m);
  REP(i,m){
    fscanf(fa,"%d %d %lf", &row[i], &col[i], &val[i] );
  }
  int n2; fscanf(fb,"%d", &n2);
  if(n2 != n) {cout << "chyba" << endl; return 1;}

  REP(i,n){
    fscanf(fb,"%lf", &b[i]);
  }

  //gauss();
  cout << "------" << endl;
  //descent_2D();
  cout << "------" << endl;
  gradient_2D();

  return 0;
}
