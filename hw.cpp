#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <vector>

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
void deba(int * a, int n){ cerr << "| "; REP(i,n) cerr << a[i] << " "; cerr << "|" << endl;}


inline bool EQ(double a, double b) { return fabs(a-b) < 1e-9; }

const int INF = 1<<30;
typedef long long ll;
typedef unsigned long long ull;
/*******************************************************/

const int M = 2600000;
const int N = 200000;

const char FILE_A[] = "smallA.txt";
const char FILE_B[] = "smallB.txt";

//const char FILE_A[] = "A.txt";
//const char FILE_B[] = "B.txt";


int n,m;


int l[M];
int r[M];
double b[N];
double sol[N];
double val[M];

double x[2][N]; //todo try [N][2]
double r[2][N];

const int N2D = 9000;
double mat[N2D][N2D];

void fill_in_2D(){
    if(n > N2D) { cout << "To big matrix" << endl; return; }
    REP(i,n)REP(j,n) mat[i][j] = 0;
    REP(i, m) mat[l[i]][r[i]] = val[i];
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
    REP(i,n) cout << b[i] << " "; cout << endl;

    FORD(i,n-1,0){
        double v = b[i];
        FORD(j,n-1,i+1) v -= sol[j] * mat[i][j];
        sol[i] = v / mat[i][i];
    }

    REP(i,n) cout << sol[i] << " "; cout << endl;


}


void pv( int [] v, int num){
    REP(i,num) cout << v[i] << " "; cout << endl;
}


inline void get_random_x(){
    REP(i,n) x[i] = rand() / (double)RAND_MAX;
    pv(x,n);
}



void descent_2D(){

    fill_in_2D();







}



int main() {

  srand( time(NULL) );

  FILE * fa = fopen(FILE_A,"r+");
  FILE * fb = fopen(FILE_B,"r+");

  fscanf(fa,"%d%d", &n,&m);
  REP(i,m){
    fscanf(fa,"%d %d %lf", &l[i], &r[i], &val[i] );
  }
  int n2; fscanf(fb,"%d", &n2);
  if(n2 != n) {cout << "chyba" << endl; return 1;}

  REP(i,n){
    fscanf(fb,"%lf", &b[i]);
  }

  gauss();


  return 0;
}
