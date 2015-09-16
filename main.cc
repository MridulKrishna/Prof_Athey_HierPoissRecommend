#include "env.hh"
#include "hgaprec.hh"
#include "ratings.hh"
#include "log.hh"

#include <stdlib.h>
#include <string>
#include <sstream>
#include <signal.h>
#include <iostream>
#include <chrono>

using namespace std::chrono;

string Env::prefix = "";
Logger::Level Env::level = Logger::DEBG;
FILE *Env::_plogf = NULL;
void usage();
void test();

Env *env_global = NULL;
volatile sig_atomic_t sig_handler_active = 0;

void
term_handler(int sig)
{
  if (env_global) {
    printf("Got signal. Saving model state.\n");
    fflush(stdout);
    env_global->save_state_now = 1;
  } else {
    signal(sig, SIG_DFL);
    raise(sig);
  }
}

int mainArray(int argc, char **argv) {
  Array arr(5);
  arr.get(0) = 3;
  arr.get(1) = 2;
  arr.get(2) = 4;
  arr.get(3) = 5;
  arr.get(4) = 8;
  
  arr.print();
  
  Array arr2(1);
  
  arr2.copy_from(arr.subarray(2,2));
  
  arr2.print();
  
  Matrix mat(2,3);
  
  mat.get(0,0) = 1;
  mat.get(0,1) = 2;
  mat.get(0,2) = 3;
  mat.get(1,0) = 4;
  mat.get(1,1) = 5;
  mat.get(1,2) = 6;
  
  mat.print();
  
  Array cs(3);
  mat.colsum(cs);
  cs.print();
  
  Array rs(2);
  mat.rowsum(rs);
  rs.print();
  
  
  return 0;
  
}

high_resolution_clock::time_point now() {
  return high_resolution_clock::now();
}

void timing(high_resolution_clock::time_point tIni, high_resolution_clock::time_point tFin) {
  duration<double> time_span = duration_cast<duration<double> >(tFin - tIni);
  cout << time_span.count() << endl;
}

int main(int argc, char **argv) {
  signal(SIGTERM, term_handler);
  if (argc <= 1) {
    printf("gaprec -dir <netflix-dataset-dir> -n <users>" \
           "-m <movies> -k <dims> -label <out-dir-tag>\n");
    exit(0);
  }
  
  double offset = 1;
  
  string fname;
  string outfname = "";
  
  uint32_t n = 0, m = 0;
  uint32_t k = 0;
  uint32_t uc = 0, ic = 0;
  uint32_t rfreq = 10;
  uint32_t max_iterations = 1000;
  double rand_seed = 0;
  
  double a=0.3, ap=1.5, bp=0.3, c=0.3, cp=1.5, dp=0.3, e=0.3, f=0.3;
  
  int scale = Env::MEAN;
  double scaleFactor = 1;
  
  uint32_t i = 0;
  
  bool lfirst = false;
  bool ofirst = false;
  
  // Parse parameters
  while (i <= argc - 1) {
    if (strcmp(argv[i], "-dir") == 0) {
      fname = string(argv[++i]);
      fprintf(stdout, "+ dir = %s\n", fname.c_str());
    } else if (strcmp(argv[i], "-n") == 0) {
      n = atoi(argv[++i]);
      fprintf(stdout, "+ n = %d\n", n);
    } else if (strcmp(argv[i], "-m") == 0) {
      m = atoi(argv[++i]);
      fprintf(stdout, "+ m = %d\n", m);
    } else if (strcmp(argv[i], "-k") == 0) {
      k = atoi(argv[++i]);
      fprintf(stdout, "+ k = %d\n", k);
    } else if (strcmp(argv[i], "-uc") == 0) {
      uc = atoi(argv[++i]);
      fprintf(stdout, "+ uc = %d\n", uc);
    } else if (strcmp(argv[i], "-ic") == 0) {
      ic = atoi(argv[++i]);
      fprintf(stdout, "+ ic = %d\n", ic);
    } else if (strcmp(argv[i], "-rfreq") == 0) {
      rfreq = atoi(argv[++i]);
      fprintf(stdout, "+ rfreq = %d\n", rfreq);
    } else if (strcmp(argv[i], "-max-iterations") == 0) {
      max_iterations = atoi(argv[++i]);
      fprintf(stdout, "+ max iterations %d\n", max_iterations);
    } else if (strcmp(argv[i], "-seed") == 0) {
      rand_seed = atof(argv[++i]);
      fprintf(stdout, "+ random seed set to %.5f\n", rand_seed);
    } else if (strcmp(argv[i], "-a") == 0) {
      a = atof(argv[++i]);
    } else if (strcmp(argv[i], "-ap") == 0) {
      ap = atof(argv[++i]);
    } else if (strcmp(argv[i], "-bp") == 0) {
      bp = atof(argv[++i]);
    } else if (strcmp(argv[i], "-c") == 0) {
      c = atof(argv[++i]);
    } else if (strcmp(argv[i], "-cp") == 0) {
      cp = atof(argv[++i]);
    } else if (strcmp(argv[i], "-dp") == 0) {
      dp = atof(argv[++i]);
    } else if (strcmp(argv[i], "-e") == 0) {
      e = atof(argv[++i]);
    } else if (strcmp(argv[i], "-f") == 0) {
      f = atof(argv[++i]);
    } else if (strcmp(argv[i], "-nooffset") == 0) {
      offset = 0;
    } else if (strcmp(argv[i], "-offset") == 0) {
      offset = atof(argv[++i]);
    } else if (strcmp(argv[i], "-outdir") == 0) {
      outfname = string(argv[++i]);
      fprintf(stdout, "+ output dir = %s\n", outfname.c_str());
    } else if (strcmp(argv[i], "-std") == 0) {
      scale = Env::STD;
    } else if (strcmp(argv[i], "-ones") == 0) {
      scale = Env::ONES;
    } else if (strcmp(argv[i], "-scfact") == 0) {
      scaleFactor = atof(argv[++i]);
    } else if (strcmp(argv[i], "-lfirst") == 0) {
      lfirst = 1;
    } else if (strcmp(argv[i], "-ofirst") == 0) {
      ofirst = 1;
    } else if (i > 0) {
      fprintf(stdout,  "error: unknown option %s\n", argv[i]);
      assert(0);
    }
    ++i;
  };
  
  if ( outfname.compare("") == 0 ) {
    outfname = fname;
  }
    
  // Initializes the environment: variables to run the code
  Env env(n, m, k, uc, ic, fname, outfname, rfreq, rand_seed, max_iterations, a, ap, bp, c, cp, dp, e, f, offset, scale, scaleFactor, lfirst, ofirst);
  env_global = &env;
  
  // Reads the input files
  Ratings ratings(env);
  if (ratings.read(fname.c_str()) < 0) {
    fprintf(stderr, "error reading dataset from dir %s; quitting\n",
            fname.c_str());
    return -1;
  }

  ratings.readValidationAndTest(fname.c_str());
  ratings.readObserved(fname.c_str());
  
  bp = sqrt(ap*cp*a*c/(ap-1)/(cp-1)*(k+uc+ic)*n*m/ratings.totRating);
  dp = bp;

  env.bp = bp;
  env.dp = dp;
  
//  ratings._userObsScale.print();
//  ratings._itemObsScale.print();
  
  cout << "Constructing hgaprec" << endl;
  HGAPRec hgaprec(env, ratings);
  
  typedef high_resolution_clock::time_point moment;
  
  moment t1 = now();
  
  cout << "Running vb_hier()" << endl;
  hgaprec.vb_hier();
  
  moment t2 = now();
  
  timing(t1,t2);
  
  return 0; 
}