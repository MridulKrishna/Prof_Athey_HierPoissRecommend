#ifndef ENV_HH
#define ENV_HH

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include <unistd.h>
#include <stdint.h>
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <map>
#include <list>
#include <vector>
#include<sys/stat.h>
#include "matrix.hh"
#include "log.hh"

typedef uint16_t yval_t;

typedef D2Array<yval_t> AdjMatrix;
typedef D2Array<double> Matrix;
typedef D3Array<double> D3;
typedef D2Array<KV> MatrixKV;
typedef D1Array<KV> KVArray;
typedef D1Array<KVI> KVIArray;

typedef std::map<uint64_t, yval_t> RatingMap;
typedef std::map<uint64_t, uint64_t> IDMap;
typedef std::map<uint32_t, uint32_t> FreqMap;
typedef std::map<string, uint32_t> FreqStrMap;
typedef std::map<string, uint32_t> StrMap;
typedef std::map<uint32_t, string> StrMapInv;

typedef D1Array<std::vector<uint32_t> *> SparseMatrix;
typedef D1Array<RatingMap *> SparseMatrixR;
typedef std::vector<Rating> RatingList;
typedef std::map<uint32_t, bool> UserMap;
typedef std::map<uint32_t, bool> MovieMap;
typedef std::map<uint32_t, bool> BoolMap;
typedef std::map<uint32_t, double> DoubleMap;
typedef std::map<uint32_t, Array *> ArrayMap;
typedef std::map<uint32_t, uint32_t> ValMap;
typedef std::map<uint32_t, vector<uint32_t> > MapVec;
typedef MapVec SparseMatrix2;
typedef std::map<Rating, bool> SampleMap;
typedef std::map<Rating, int> CountMap;
typedef std::map<Rating, double> ValueMap;
typedef std::map<uint32_t, string> StrMapInv;

class Env {
public:
  typedef enum { NETFLIX, MOVIELENS, MENDELEY, ECHONEST, NYT } Dataset;
  typedef enum { CREATE_TRAIN_TEST_SETS, TRAINING } Mode;
  Env(uint32_t N, uint32_t M, uint32_t K, uint32_t UC, uint32_t IC, string fname, string outfname, uint32_t rfreq, double rseed,
      uint32_t max_iterations, double Na, double Nap, double Nbp, double Nc, double Ncp, double Ndp, double Ne, double Nf,double pOffset, int scale, double scaleFactor, bool nLfirst, bool nOfirst, bool nSession, bool nFitpriors);
  ~Env() { fclose(_plogf); }
  
  static string prefix;
  
  static Logger::Level level;
  
  Dataset dataset;
  uint32_t n;  // users
  uint32_t m;  // movies
  uint32_t nTrain; // users in train set
  uint32_t mTrain; // movies in train set
  uint32_t k;  // factors
  uint32_t uc;  // Number of user characteristics
  uint32_t ic;  // Number of item characteristics
  uint32_t t;
  uint32_t mini_batch_size;
  double offset;
  
  // Hyperparameters
  double a, ap, bp, c, cp, dp, e, f;
  
  double alpha;
  double tau0;
  double tau1;
  double heldout_ratio;
  double validation_ratio;
  int reportfreq;
  double epsilon;
  double logepsilon;
  bool nolambda;
  bool strid;
  bool logl;
  uint32_t max_iterations;
  double seed;
  bool save_state_now;
  string datfname;
  string outfname;
  string label;
  bool nmi;
  string ground_truth_fname;
  bool model_load;
  string model_location;
  bool gen_heldout;
  uint32_t online_iterations;
  double meanchangethresh;
  bool batch;
  Mode mode;
  bool binary_data;
  bool bias;
  bool hier;
  bool explore;
  bool vb;
  bool nmf;
  bool nmfload;
  bool lda;
  bool vwlda;
  bool write_training;
  uint32_t rating_threshold;
  bool graphchi;
  bool wals;
  double wals_l;
  uint32_t wals_C;
  bool als;
  bool chinmf;
  bool climf;
  
  bool mle_item;
  bool mle_user;
  bool canny;
  bool ctr;
  
  int scale;
  double scaleFactor;
 
  bool lfirst;
  bool ofirst;
  
  bool session;
  
  bool fitpriors;
  
  static const int ONES = 1;
  static const int MEAN = 2;
  static const int STD = 3;
  
  template<class T> static void plog(string s, const T &v);
  static string file_str(string fname);
  static string outfile_str(string fname);
  
private:
  static FILE *_plogf;
};


template<class T> inline void
Env::plog(string s, const T &v)
{
  fprintf(_plogf, "%s: %s\n", s.c_str(), v.s().c_str());
  fflush(_plogf);
}

template<> inline void
Env::plog(string s, const double &v)
{
  fprintf(_plogf, "%s: %.9f\n", s.c_str(), v);
  fflush(_plogf);
}

template<> inline void
Env::plog(string s, const string &v)
{
  fprintf(_plogf, "%s: %s\n", s.c_str(), v.c_str());
  fflush(_plogf);
}

template<> inline void
Env::plog(string s, const bool &v)
{
  fprintf(_plogf, "%s: %s\n", s.c_str(), v ? "True": "False");
  fflush(_plogf);
}

template<> inline void
Env::plog(string s, const int &v)
{
  fprintf(_plogf, "%s: %d\n", s.c_str(), v);
  fflush(_plogf);
}

template<> inline void
Env::plog(string s, const unsigned &v)
{
  fprintf(_plogf, "%s: %d\n", s.c_str(), v);
  fflush(_plogf);
}

template<> inline void
Env::plog(string s, const short unsigned int &v)
{
  fprintf(_plogf, "%s: %d\n", s.c_str(), v);
  fflush(_plogf);
}

template<> inline void
Env::plog(string s, const uint64_t &v)
{
  fprintf(_plogf, "%s: %" PRIu64 "\n", s.c_str(), v);
  fflush(_plogf);
}

#ifdef __APPLE__
template<> inline void
Env::plog(string s, const long unsigned int &v)
{
  fprintf(_plogf, "%s: %lu\n", s.c_str(), v);
  fflush(_plogf);
}
#endif

inline string
Env::file_str(string fname)
{
  string s = prefix + fname;
  return s;
}

inline string
Env::outfile_str(string fname)
{
  string s = prefix + fname;
  return s;
}

inline
Env::Env(uint32_t N, uint32_t M, uint32_t K, uint32_t UC, uint32_t IC, string Nfname, string Noutfname, uint32_t rfreq,double rseed,
         uint32_t max_iterations, double Na, double Nap, double Nbp, double Nc, double Ncp, double Ndp, double Ne, double Nf,double pOffset, int nScale, double nScaleFactor, bool nLfirst, bool nOfirst, bool nSession, bool nFitpriors)
: n(N),
m(M),
k(K),
uc(UC),
ic(IC),
t(2),
mini_batch_size(1000),
a(Na), ap(Nap), bp(Nbp), c(Nc), cp(Ncp), dp(Ndp), e(Ne), f(Nf),
tau0(0),
tau1(0),
heldout_ratio(0.2),
validation_ratio(0.01),
reportfreq(rfreq),
epsilon(0.001),
logepsilon(log(epsilon)),
nolambda(true),
max_iterations(max_iterations),
seed(rseed),
save_state_now(false),
datfname(Nfname),
outfname(Noutfname),
online_iterations(1),
meanchangethresh(0.001),
mode(TRAINING),
offset(pOffset),
scale(nScale),
scaleFactor(nScaleFactor),
lfirst(nLfirst),
ofirst(nOfirst),
session(nSession),
fitpriors(nFitpriors)
{
  ostringstream sa;
  sa << "n" << n << "-";
  sa << "m" << m << "-";
  sa << "k" << k << "-";
  sa << "uc" << uc << "-";
  sa << "ic" << ic;
  if (label != "")
    sa << "-" << label;
  else if (datfname.length() > 3) {
    string q = datfname.substr(0,2);
    if (isalpha(q[0]))
      sa << "-" << q;
  }
  if (scale != MEAN) {
    if (scale == STD) {
      sa << "-std";
    }
    else if (scale == ONES) {
      sa << "-ones";
    }
  }
  
  if (scaleFactor != 1) {
    sa << "-scfact" << scaleFactor;
  }
  
  if (offset != 1) {
    sa << "-offset" << offset;
  }
  
  if (lfirst) {
    sa << "-lfirst";
  }
  
  if (ofirst) {
    sa << "-ofirst";
  }
  
  if (session) {
    sa << "-session";
  }
  
  if (fitpriors) {
    sa << "-fpriors";
  }

  prefix = sa.str();
  level = Logger::TEST;
  
  fprintf(stdout, "+ Creating directory %s\n", prefix.c_str());
  fflush(stdout);
  
  struct stat buffer;
  
  int out = stat((outfname+"/"+prefix).c_str(), &buffer);
  
  if ( out != 0 ) {
    cout << "+ Creating directory " << outfname << "/" << prefix << endl;
    cout << "HI 1 2 3 ... " << (outfname+"/"+prefix).c_str() << " " << endl;
    assert(!mkdir((outfname+"/"+prefix).c_str(), S_IRUSR | S_IWUSR | S_IXUSR));
  }
  
  assert (Logger::initialize(prefix, "infer.log",
                             true, level) >= 0);

  string name = outfname+"/"+prefix +"/param.txt";
//  cout << name << endl;
  _plogf = fopen(name.c_str(), "w");
  if (!_plogf)  {
    printf("cannot open param file:%s\n",  strerror(errno));
    exit(-1);
  }
  hier = true; //To compute log-likelihood for the HPF algorithm
 
  plog("n", n);
  plog("k", k);
  plog("t", t);
  plog("test_ratio", heldout_ratio);
  plog("validation_ratio", validation_ratio);
  plog("seed", seed);
  plog("a", a);
  plog("ap", ap);
  plog("bp", bp);
  plog("c", c);
  plog("cp", cp);
  plog("dp", dp);
  plog("e", e);
  plog("f", f);
  plog("reportfreq", reportfreq);
  plog("vb", vb);
  plog("bias", bias);
  plog("hier", hier);
  plog("nmf", nmf);
  plog("lda", lda);
  plog("wals_l", wals_l);
  plog("wals_C", wals_C);
  plog("mle_user", mle_user);
  plog("mle_item", mle_item);
  
  //string ndatfname = file_str("/network.dat");
  //unlink(ndatfname.c_str());
  //assert (symlink(datfname.c_str(), ndatfname.c_str()) >= 0);
  //unlink(file_str("/mutual.txt").c_str());
}

/*
 src: http://www.delorie.com/gnu/docs/glibc/libc_428.html
 Subtract the `struct timeval' values X and Y,
 storing the result in RESULT.
 Return 1 if the difference is negative, otherwise 0.
 */
inline int
timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y)
{
  /* Perform the carry for the later subtraction by updating y. */
  if (x->tv_usec < y->tv_usec) {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }
  if (x->tv_usec - y->tv_usec > 1000000) {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;
    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }
  
  /* Compute the time remaining to wait.
   tv_usec is certainly positive. */
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;
  
  /* Return 1 if result is negative. */
  return x->tv_sec < y->tv_sec;
}

inline void
timeval_add (struct timeval *result, const struct timeval *x)
{
  result->tv_sec  += x->tv_sec;
  result->tv_usec += x->tv_usec;
  
  if (result->tv_usec >= 1000000) {
    result->tv_sec++;
    result->tv_usec -= 1000000;
  }
}

#endif
