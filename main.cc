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

int mainRepeatedScale(int argc, char **argv) {
  double scaleFs[6];
  int J = 6;
  scaleFs[0] = 1;
  scaleFs[1] = 2;
  scaleFs[2] = 3;
  scaleFs[3] = 5;
  scaleFs[4] = 8;
  scaleFs[5] = 10;
  
  int offset = 1;
  
  string fname;
  string outfname = "";
  
  uint32_t n = 0, m = 0;
  uint32_t k = 0;
  uint32_t uc = 0, ic = 0;
  string ground_truth_fname;
  uint32_t rfreq = 10;
  string label;
  bool logl = false;
  uint32_t max_iterations = 1000;
  bool nmi = false;
  bool strid = false;
  double rand_seed = 0;
  
  bool test = false;
  bool batch = true;
  bool online = false;
  bool gen_heldout = false;
  
  bool model_load = false;
  string model_location = "";
  
  bool hol_load = false;
  string hol_location = "";
  bool vb = true;
  
  bool pred_accuracy = false;
  bool gt_accuracy = false;
  bool p = false;
  double a=0.3, ap=0.3, bp=0.3, c=0.3, cp=0.3, dp=0.3, e=0.3, f=0.3;
  Env::Dataset dataset = Env::MENDELEY;
  bool binary_data = false;
  bool bias = false;
  bool hier = false;
  bool explore = false;
  bool gen_ranking_for_users = false;
  bool rmse = false;
  bool nmf = false;
  bool nmfload = false;
  bool vwload = false;
  bool lda = false;
  bool vwlda = false;
  bool msr = false;
  bool write_training = false;
  uint32_t rating_threshold = 1;
  bool chi = false;
  bool wals = false;
  double wals_l = 0.1;
  uint32_t wals_C = 10;
  
  bool als = false;
  bool chinmf = false;
  bool climf = false;
  
  bool mle_item = false;
  bool mle_user = false;
  bool canny = false;
  bool ctr = false;
  
  uint32_t i = 0;
  int scale = Env::MEAN;
  double scaleFactor = 1;
  
  // Parse parameters
  while (i <= argc - 1) {
    if (strcmp(argv[i], "-dir") == 0) {
      fname = string(argv[++i]);
      fprintf(stdout, "+ dir = %s\n", fname.c_str());
    } else if (strcmp(argv[i], "-n") == 0) {
      n = atoi(argv[++i]);
      fprintf(stdout, "+ n = %d\n", n);
    } else if (strcmp(argv[i], "-p") == 0) {
      p = true;
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
    } else if (strcmp(argv[i], "-nmi") == 0) {
      ground_truth_fname = string(argv[++i]);
      fprintf(stdout, "+ ground truth fname = %s\n",
              ground_truth_fname.c_str());
      nmi = true;
    } else if (strcmp(argv[i], "-rfreq") == 0) {
      rfreq = atoi(argv[++i]);
      fprintf(stdout, "+ rfreq = %d\n", rfreq);
    } else if (strcmp(argv[i], "-strid") == 0) {
      strid = true;
      fprintf(stdout, "+ strid mode\n");
    } else if (strcmp(argv[i], "-label") == 0) {
      label = string(argv[++i]);
    } else if (strcmp(argv[i], "-logl") == 0) {
      logl = true;
      fprintf(stdout, "+ logl mode\n");
    } else if (strcmp(argv[i], "-max-iterations") == 0) {
      max_iterations = atoi(argv[++i]);
      fprintf(stdout, "+ max iterations %d\n", max_iterations);
    } else if (strcmp(argv[i], "-seed") == 0) {
      rand_seed = atof(argv[++i]);
      fprintf(stdout, "+ random seed set to %.5f\n", rand_seed);
    } else if (strcmp(argv[i], "-load") == 0) {
      model_load = true;
      model_location = string(argv[++i]);
      fprintf(stdout, "+ loading theta from %s\n", model_location.c_str());
    } else if (strcmp(argv[i], "-test") == 0) {
      test = true;
      fprintf(stdout, "+ test mode\n");
    } else if (strcmp(argv[i], "-batch") == 0) {
      batch = true;
      fprintf(stdout, "+ batch inference\n");
    } else if (strcmp(argv[i], "-online") == 0) {
      batch = false;
      fprintf(stdout, "+ online inference\n");
    } else if (strcmp(argv[i], "-gen-heldout") == 0) {
      gen_heldout = true;
      fprintf(stdout, "+ generate held-out files from dataset\n");
    } else if (strcmp(argv[i], "-pred-accuracy") == 0) {
      pred_accuracy = true;
      fprintf(stdout, "+ compute predictive accuracy\n");
    } else if (strcmp(argv[i], "-gt-accuracy") == 0) {
      gt_accuracy = true;
      fprintf(stdout, "+ compute  accuracy to ground truth\n");
    } else if (strcmp(argv[i], "-netflix") == 0) {
      dataset = Env::NETFLIX;
    } else if (strcmp(argv[i], "-mendeley") == 0) {
      dataset = Env::MENDELEY;
    } else if (strcmp(argv[i], "-movielens") == 0) {
      dataset = Env::MOVIELENS;
    } else if (strcmp(argv[i], "-echonest") == 0) {
      dataset = Env::ECHONEST;
    } else if (strcmp(argv[i], "-nyt") == 0) {
      dataset = Env::NYT;
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
    } else if (strcmp(argv[i], "-binary-data") == 0) {
      binary_data = true;
    } else if (strcmp(argv[i], "-bias") == 0) {
      bias = true;
    } else if (strcmp(argv[i], "-hier") == 0) {
      hier = true;
    } else if (strcmp(argv[i], "-mle-user") == 0) {
      mle_user = true;
    } else if (strcmp(argv[i], "-mle-item") == 0) {
      mle_item = true;
    } else if (strcmp(argv[i], "-canny") == 0) {
      canny = true;
    } else if (strcmp(argv[i], "-gen-ranking") == 0) {
      gen_ranking_for_users = true;
    } else if (strcmp(argv[i], "-rmse") == 0) {
      rmse = true;
    } else if (strcmp(argv[i], "-novb") == 0) {
      vb = false;
    } else if (strcmp(argv[i], "-msr") == 0) {
      msr = true;
    } else if (strcmp(argv[i], "-nmf") == 0) {
      nmf = true;
    } else if (strcmp(argv[i], "-nmfload") == 0) {
      nmfload = true;
    } else if (strcmp(argv[i], "-vwload") == 0) {
      vwload = true;
    } else if (strcmp(argv[i], "-lda") == 0) {
      lda = true;
    } else if (strcmp(argv[i], "-vwlda") == 0) {
      vwlda = true;
    } else if (strcmp(argv[i], "-write-training") == 0) {
      write_training = true;
    } else if (strcmp(argv[i], "-chi") == 0) {
      chi = true;
    } else if (strcmp(argv[i], "-chinmf") == 0) {
      chinmf = true;
    } else if (strcmp(argv[i], "-als") == 0) {
      als = true;
    } else if (strcmp(argv[i], "-wals") == 0) {
      wals = true;
    } else if (strcmp(argv[i], "-wals_l") == 0) {
      wals_l = atof(argv[++i]);
    } else if (strcmp(argv[i], "-wals_C") == 0) {
      wals_C = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-climf") == 0) {
      climf = true;
    } else if (strcmp(argv[i], "-ctr") == 0) {
      ctr = true;
    } else if (strcmp(argv[i], "-rating-threshold") == 0) {
      rating_threshold = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-nooffset") == 0) {
      offset = 0;
    } else if (strcmp(argv[i], "-outdir") == 0) {
      outfname = string(argv[++i]);
      fprintf(stdout, "+ output dir = %s\n", outfname.c_str());
    } else if (strcmp(argv[i], "-std") == 0) {
      scale = Env::STD;
    } else if (strcmp(argv[i], "-ones") == 0) {
      scale = Env::ONES;
    } else if (strcmp(argv[i], "-scfact") == 0) {
      scaleFactor = atof(argv[++i]);
    } else if (i > 0) {
      fprintf(stdout,  "error: unknown option %s\n", argv[i]);
      assert(0);
    }
    ++i;
  };
  
  scale = Env::MEAN;
  
  for ( int j = 0; j<J; ++j) {
    scaleFactor = scaleFs[j];
    if ( outfname.compare("") == 0 ) {
      outfname = fname;
    }
    
    // Initializes the environment: variables to run the code
    Env env(n, m, k, uc, ic, fname, outfname, nmi, ground_truth_fname, rfreq,
            strid, label, logl, rand_seed, max_iterations,
            model_load, model_location,
            gen_heldout, a, ap, bp, c, cp, dp, e, f, dataset,
            batch, binary_data, bias, hier,
            explore, vb, nmf, nmfload, lda, vwlda,
            write_training, rating_threshold,
            chi, wals, wals_l, wals_C,
            als, chinmf, climf,
            mle_item, mle_user, canny, ctr, offset, scale, scaleFactor);
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
    
    cout << "Constructing hgaprec" << endl;
    HGAPRec hgaprec(env, ratings);
    
    typedef high_resolution_clock::time_point moment;
    
    moment t1 = now();
    
    cout << "Running vb_hier()" << endl;
    hgaprec.vb_hier();
    cout << "\n" << endl;
    moment t2 = now();
    
    timing(t1,t2);
  }
  
  scale = Env::STD;
  
  for ( int j = 0; j<J; ++j) {
    double scaleFactor = scaleFs[j];
    if ( outfname.compare("") == 0 ) {
      outfname = fname;
    }
    
    // Initializes the environment: variables to run the code
    Env env(n, m, k, uc, ic, fname, outfname, nmi, ground_truth_fname, rfreq,
            strid, label, logl, rand_seed, max_iterations,
            model_load, model_location,
            gen_heldout, a, ap, bp, c, cp, dp, e, f, dataset,
            batch, binary_data, bias, hier,
            explore, vb, nmf, nmfload, lda, vwlda,
            write_training, rating_threshold,
            chi, wals, wals_l, wals_C,
            als, chinmf, climf,
            mle_item, mle_user, canny, ctr, offset, scale, scaleFactor);
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
    
    cout << "Constructing hgaprec" << endl;
    HGAPRec hgaprec(env, ratings);
    
    typedef high_resolution_clock::time_point moment;
    
    moment t1 = now();
    
    cout << "Running vb_hier()" << endl;
    hgaprec.vb_hier();
    cout << "\n" << endl;
    moment t2 = now();
    
    timing(t1,t2);
  }
  
  scale = Env::ONES;
  
  for ( int j = 0; j<J; ++j) {
    double scaleFactor = scaleFs[j];
    if ( outfname.compare("") == 0 ) {
      outfname = fname;
    }
    
    // Initializes the environment: variables to run the code
    Env env(n, m, k, uc, ic, fname, outfname, nmi, ground_truth_fname, rfreq,
            strid, label, logl, rand_seed, max_iterations,
            model_load, model_location,
            gen_heldout, a, ap, bp, c, cp, dp, e, f, dataset,
            batch, binary_data, bias, hier,
            explore, vb, nmf, nmfload, lda, vwlda,
            write_training, rating_threshold,
            chi, wals, wals_l, wals_C,
            als, chinmf, climf,
            mle_item, mle_user, canny, ctr, offset, scale, scaleFactor);
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
    
    cout << "Constructing hgaprec" << endl;
    HGAPRec hgaprec(env, ratings);
    
    typedef high_resolution_clock::time_point moment;
    
    moment t1 = now();
    
    cout << "Running vb_hier()" << endl;
    hgaprec.vb_hier();
    cout << "\n" << endl;
    moment t2 = now();
    
    timing(t1,t2);
  }
  
  return 0;
}

int mainOnce(int argc, char **argv) {
  signal(SIGTERM, term_handler);
  if (argc <= 1) {
    printf("gaprec -dir <netflix-dataset-dir> -n <users>" \
           "-m <movies> -k <dims> -label <out-dir-tag>\n");
    exit(0);
  }
  
  int offset = 1;
  
  string fname;
  string outfname = "";
  
  uint32_t n = 0, m = 0;
  uint32_t k = 0;
  uint32_t uc = 0, ic = 0;
  string ground_truth_fname;
  uint32_t rfreq = 10;
  string label;
  bool logl = false;
  uint32_t max_iterations = 1000;
  bool nmi = false;
  bool strid = false;
  double rand_seed = 0;
  
  bool test = false;
  bool batch = true;
  bool online = false;
  bool gen_heldout = false;
  
  bool model_load = false;
  string model_location = "";
  
  bool hol_load = false;
  string hol_location = "";
  bool vb = true;
  
  bool pred_accuracy = false;
  bool gt_accuracy = false;
  bool p = false;
  double a=0.3, ap=0.3, bp=0.3, c=0.3, cp=0.3, dp=0.3, e=0.3, f=0.3;
  Env::Dataset dataset = Env::MENDELEY;
  bool binary_data = false;
  bool bias = false;
  bool hier = false;
  bool explore = false;
  bool gen_ranking_for_users = false;
  bool rmse = false;
  bool nmf = false;
  bool nmfload = false;
  bool vwload = false;
  bool lda = false;
  bool vwlda = false;
  bool msr = false;
  bool write_training = false;
  uint32_t rating_threshold = 1;
  bool chi = false;
  bool wals = false;
  double wals_l = 0.1;
  uint32_t wals_C = 10;
  
  bool als = false;
  bool chinmf = false;
  bool climf = false;
  
  bool mle_item = false;
  bool mle_user = false;
  bool canny = false;
  bool ctr = false;
  
  int scale = Env::MEAN;
  double scaleFactor = 1;
  
  uint32_t i = 0;
  
  // Parse parameters
  while (i <= argc - 1) {
    if (strcmp(argv[i], "-dir") == 0) {
      fname = string(argv[++i]);
      fprintf(stdout, "+ dir = %s\n", fname.c_str());
    } else if (strcmp(argv[i], "-n") == 0) {
      n = atoi(argv[++i]);
      fprintf(stdout, "+ n = %d\n", n);
    } else if (strcmp(argv[i], "-p") == 0) {
      p = true;
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
    } else if (strcmp(argv[i], "-nmi") == 0) {
      ground_truth_fname = string(argv[++i]);
      fprintf(stdout, "+ ground truth fname = %s\n",
              ground_truth_fname.c_str());
      nmi = true;
    } else if (strcmp(argv[i], "-rfreq") == 0) {
      rfreq = atoi(argv[++i]);
      fprintf(stdout, "+ rfreq = %d\n", rfreq);
    } else if (strcmp(argv[i], "-strid") == 0) {
      strid = true;
      fprintf(stdout, "+ strid mode\n");
    } else if (strcmp(argv[i], "-label") == 0) {
      label = string(argv[++i]);
    } else if (strcmp(argv[i], "-logl") == 0) {
      logl = true;
      fprintf(stdout, "+ logl mode\n");
    } else if (strcmp(argv[i], "-max-iterations") == 0) {
      max_iterations = atoi(argv[++i]);
      fprintf(stdout, "+ max iterations %d\n", max_iterations);
    } else if (strcmp(argv[i], "-seed") == 0) {
      rand_seed = atof(argv[++i]);
      fprintf(stdout, "+ random seed set to %.5f\n", rand_seed);
    } else if (strcmp(argv[i], "-load") == 0) {
      model_load = true;
      model_location = string(argv[++i]);
      fprintf(stdout, "+ loading theta from %s\n", model_location.c_str());
    } else if (strcmp(argv[i], "-test") == 0) {
      test = true;
      fprintf(stdout, "+ test mode\n");
    } else if (strcmp(argv[i], "-batch") == 0) {
      batch = true;
      fprintf(stdout, "+ batch inference\n");
    } else if (strcmp(argv[i], "-online") == 0) {
      batch = false;
      fprintf(stdout, "+ online inference\n");
    } else if (strcmp(argv[i], "-gen-heldout") == 0) {
      gen_heldout = true;
      fprintf(stdout, "+ generate held-out files from dataset\n");
    } else if (strcmp(argv[i], "-pred-accuracy") == 0) {
      pred_accuracy = true;
      fprintf(stdout, "+ compute predictive accuracy\n");
    } else if (strcmp(argv[i], "-gt-accuracy") == 0) {
      gt_accuracy = true;
      fprintf(stdout, "+ compute  accuracy to ground truth\n");
    } else if (strcmp(argv[i], "-netflix") == 0) {
      dataset = Env::NETFLIX;
    } else if (strcmp(argv[i], "-mendeley") == 0) {
      dataset = Env::MENDELEY;
    } else if (strcmp(argv[i], "-movielens") == 0) {
      dataset = Env::MOVIELENS;
    } else if (strcmp(argv[i], "-echonest") == 0) {
      dataset = Env::ECHONEST;
    } else if (strcmp(argv[i], "-nyt") == 0) {
      dataset = Env::NYT;
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
    } else if (strcmp(argv[i], "-binary-data") == 0) {
      binary_data = true;
    } else if (strcmp(argv[i], "-bias") == 0) {
      bias = true;
    } else if (strcmp(argv[i], "-hier") == 0) {
      hier = true;
    } else if (strcmp(argv[i], "-mle-user") == 0) {
      mle_user = true;
    } else if (strcmp(argv[i], "-mle-item") == 0) {
      mle_item = true;
    } else if (strcmp(argv[i], "-canny") == 0) {
      canny = true;
    } else if (strcmp(argv[i], "-gen-ranking") == 0) {
      gen_ranking_for_users = true;
    } else if (strcmp(argv[i], "-rmse") == 0) {
      rmse = true;
    } else if (strcmp(argv[i], "-novb") == 0) {
      vb = false;
    } else if (strcmp(argv[i], "-msr") == 0) {
      msr = true;
    } else if (strcmp(argv[i], "-nmf") == 0) {
      nmf = true;
    } else if (strcmp(argv[i], "-nmfload") == 0) {
      nmfload = true;
    } else if (strcmp(argv[i], "-vwload") == 0) {
      vwload = true;
    } else if (strcmp(argv[i], "-lda") == 0) {
      lda = true;
    } else if (strcmp(argv[i], "-vwlda") == 0) {
      vwlda = true;
    } else if (strcmp(argv[i], "-write-training") == 0) {
      write_training = true;
    } else if (strcmp(argv[i], "-chi") == 0) {
      chi = true;
    } else if (strcmp(argv[i], "-chinmf") == 0) {
      chinmf = true;
    } else if (strcmp(argv[i], "-als") == 0) {
      als = true;
    } else if (strcmp(argv[i], "-wals") == 0) {
      wals = true;
    } else if (strcmp(argv[i], "-wals_l") == 0) {
      wals_l = atof(argv[++i]);
    } else if (strcmp(argv[i], "-wals_C") == 0) {
      wals_C = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-climf") == 0) {
      climf = true;
    } else if (strcmp(argv[i], "-ctr") == 0) {
      ctr = true;
    } else if (strcmp(argv[i], "-rating-threshold") == 0) {
      rating_threshold = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-nooffset") == 0) {
      offset = 0;
    } else if (strcmp(argv[i], "-outdir") == 0) {
      outfname = string(argv[++i]);
      fprintf(stdout, "+ output dir = %s\n", outfname.c_str());
    } else if (strcmp(argv[i], "-std") == 0) {
      scale = Env::STD;
    } else if (strcmp(argv[i], "-ones") == 0) {
      scale = Env::ONES;
    } else if (strcmp(argv[i], "-scfact") == 0) {
      scaleFactor = atof(argv[++i]);
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
  Env env(n, m, k, uc, ic, fname, outfname, nmi, ground_truth_fname, rfreq,
          strid, label, logl, rand_seed, max_iterations,
          model_load, model_location,
          gen_heldout, a, ap, bp, c, cp, dp, e, f, dataset,
          batch, binary_data, bias, hier,
          explore, vb, nmf, nmfload, lda, vwlda,
          write_training, rating_threshold,
          chi, wals, wals_l, wals_C,
          als, chinmf, climf,
          mle_item, mle_user, canny, ctr, offset, scale, scaleFactor);
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
  
  cout << "Constructing hgaprec" << endl;
  HGAPRec hgaprec(env, ratings);
  
  typedef high_resolution_clock::time_point moment;
  
  moment t1 = now();
  
  cout << "Running vb_hier()" << endl;
  hgaprec.vb_hier();
  cout << "\n" << endl;
  moment t2 = now();
  
  timing(t1,t2);
  
  return 0; 
}

int main(int argc, char **argv) {
  mainOnce(argc, argv);
//  mainRepeatedScale(argc, argv);
}
/*
int
mainOri(int argc, char **argv)
{
  signal(SIGTERM, term_handler);
  if (argc <= 1) {
    printf("gaprec -dir <netflix-dataset-dir> -n <users>" \
	   "-m <movies> -k <dims> -label <out-dir-tag>\n");
    exit(0);
  }

  string fname;
  uint32_t n = 0, m = 0;
  uint32_t k = 0;
  string ground_truth_fname;
  uint32_t rfreq = 10;
  string label;
  bool logl = false;
  uint32_t max_iterations = 1000;
  bool nmi = false;
  bool strid = false;
  double rand_seed = 0;

  bool test = false;
  bool batch = true;
  bool online = false;
  bool gen_heldout = false;

  bool model_load = false;
  string model_location = "";

  bool hol_load = false;
  string hol_location = "";
  bool vb = true;
  
  bool pred_accuracy = false;
  bool gt_accuracy = false;
  bool p = false;
  double a = 0.3, b = 0.3, c = 0.3, d = 0.3;
  Env::Dataset dataset = Env::MENDELEY;
  bool binary_data = false;
  bool bias = false;
  bool hier = false;
  bool explore = false;
  bool gen_ranking_for_users = false;
  bool rmse = false;
  bool nmf = false;
  bool nmfload = false;
  bool vwload = false;
  bool lda = false;
  bool vwlda = false;
  bool msr = false;
  bool write_training = false;
  uint32_t rating_threshold = 1;
  bool chi = false;
  bool wals = false;
  double wals_l = 0.1;
  uint32_t wals_C = 10;

  bool als = false;
  bool chinmf = false;
  bool climf = false;

  bool mle_item = false;
  bool mle_user = false;
  bool canny = false;
  bool ctr = false;

  uint32_t i = 0;
    
  // Parse parameters
  while (i <= argc - 1) {
    if (strcmp(argv[i], "-dir") == 0) {
      fname = string(argv[++i]);
      fprintf(stdout, "+ dir = %s\n", fname.c_str());
    } else if (strcmp(argv[i], "-n") == 0) {
      n = atoi(argv[++i]);
      fprintf(stdout, "+ n = %d\n", n);
    } else if (strcmp(argv[i], "-p") == 0) {
      p = true;
    } else if (strcmp(argv[i], "-m") == 0) {
      m = atoi(argv[++i]);
      fprintf(stdout, "+ m = %d\n", m);
    } else if (strcmp(argv[i], "-k") == 0) {
      k = atoi(argv[++i]);
      fprintf(stdout, "+ k = %d\n", k);
    } else if (strcmp(argv[i], "-nmi") == 0) {
      ground_truth_fname = string(argv[++i]);
      fprintf(stdout, "+ ground truth fname = %s\n",
	      ground_truth_fname.c_str());
      nmi = true;
    } else if (strcmp(argv[i], "-rfreq") == 0) {
      rfreq = atoi(argv[++i]);
      fprintf(stdout, "+ rfreq = %d\n", rfreq);
    } else if (strcmp(argv[i], "-strid") == 0) {
      strid = true;
      fprintf(stdout, "+ strid mode\n");
    } else if (strcmp(argv[i], "-label") == 0) {
      label = string(argv[++i]);
    } else if (strcmp(argv[i], "-logl") == 0) {
      logl = true;
      fprintf(stdout, "+ logl mode\n");
    } else if (strcmp(argv[i], "-max-iterations") == 0) {
      max_iterations = atoi(argv[++i]);
      fprintf(stdout, "+ max iterations %d\n", max_iterations);
    } else if (strcmp(argv[i], "-seed") == 0) {
      rand_seed = atof(argv[++i]);
      fprintf(stdout, "+ random seed set to %.5f\n", rand_seed);
    } else if (strcmp(argv[i], "-load") == 0) {
      model_load = true;
      model_location = string(argv[++i]);
      fprintf(stdout, "+ loading theta from %s\n", model_location.c_str());
    } else if (strcmp(argv[i], "-test") == 0) {
      test = true;
      fprintf(stdout, "+ test mode\n");
    } else if (strcmp(argv[i], "-batch") == 0) {
      batch = true;
      fprintf(stdout, "+ batch inference\n");
    } else if (strcmp(argv[i], "-online") == 0) {
      batch = false;
      fprintf(stdout, "+ online inference\n");
    } else if (strcmp(argv[i], "-gen-heldout") == 0) {
      gen_heldout = true;
      fprintf(stdout, "+ generate held-out files from dataset\n");
    } else if (strcmp(argv[i], "-pred-accuracy") == 0) {
      pred_accuracy = true;
      fprintf(stdout, "+ compute predictive accuracy\n");
    } else if (strcmp(argv[i], "-gt-accuracy") == 0) {
      gt_accuracy = true;
      fprintf(stdout, "+ compute  accuracy to ground truth\n");
    } else if (strcmp(argv[i], "-netflix") == 0) {
      dataset = Env::NETFLIX;
    } else if (strcmp(argv[i], "-mendeley") == 0) {
      dataset = Env::MENDELEY;
    } else if (strcmp(argv[i], "-movielens") == 0) {
      dataset = Env::MOVIELENS;
    } else if (strcmp(argv[i], "-echonest") == 0) {
      dataset = Env::ECHONEST;
    } else if (strcmp(argv[i], "-nyt") == 0) {
      dataset = Env::NYT;
    } else if (strcmp(argv[i], "-a") == 0) {
      a = atof(argv[++i]);
    } else if (strcmp(argv[i], "-b") == 0) {
      b = atof(argv[++i]);
    } else if (strcmp(argv[i], "-c") == 0) {
      c = atof(argv[++i]);
    } else if (strcmp(argv[i], "-d") == 0) {
      d = atof(argv[++i]);
    } else if (strcmp(argv[i], "-binary-data") == 0) {
      binary_data = true;
    } else if (strcmp(argv[i], "-bias") == 0) {
      bias = true;
    } else if (strcmp(argv[i], "-hier") == 0) {
      hier = true;
    } else if (strcmp(argv[i], "-mle-user") == 0) {
      mle_user = true;
    } else if (strcmp(argv[i], "-mle-item") == 0) {
      mle_item = true;
    } else if (strcmp(argv[i], "-canny") == 0) {
      canny = true;
    } else if (strcmp(argv[i], "-gen-ranking") == 0) {
      gen_ranking_for_users = true;
    } else if (strcmp(argv[i], "-rmse") == 0) {
      rmse = true;
    } else if (strcmp(argv[i], "-novb") == 0) {
      vb = false;
    } else if (strcmp(argv[i], "-msr") == 0) {
      msr = true;
    } else if (strcmp(argv[i], "-nmf") == 0) {
      nmf = true;
    } else if (strcmp(argv[i], "-nmfload") == 0) {
      nmfload = true;
    } else if (strcmp(argv[i], "-vwload") == 0) {
      vwload = true;
    } else if (strcmp(argv[i], "-lda") == 0) {
      lda = true;
    } else if (strcmp(argv[i], "-vwlda") == 0) {
      vwlda = true;
    } else if (strcmp(argv[i], "-write-training") == 0) {
      write_training = true;
    } else if (strcmp(argv[i], "-chi") == 0) {
      chi = true;
    } else if (strcmp(argv[i], "-chinmf") == 0) {
      chinmf = true;
    } else if (strcmp(argv[i], "-als") == 0) {
      als = true;
    } else if (strcmp(argv[i], "-wals") == 0) {
      wals = true;
    } else if (strcmp(argv[i], "-wals_l") == 0) {
      wals_l = atof(argv[++i]);
    } else if (strcmp(argv[i], "-wals_C") == 0) {
      wals_C = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-climf") == 0) {
      climf = true;
    } else if (strcmp(argv[i], "-ctr") == 0) {
      ctr = true;
    } else if (strcmp(argv[i], "-rating-threshold") == 0) {
      rating_threshold = atoi(argv[++i]);
    } else if (i > 0) {
      fprintf(stdout,  "error: unknown option %s\n", argv[i]);
      assert(0);
    } 
    ++i;
  };
  
  // Initializes the environment: variables to run the code
  Env env(n, m, k, fname, nmi, ground_truth_fname, rfreq, 
	  strid, label, logl, rand_seed, max_iterations, 
	  model_load, model_location, 
	  gen_heldout, a, b, c, d, dataset, 
	  batch, binary_data, bias, hier, 
	  explore, vb, nmf, nmfload, lda, vwlda, 
	  write_training, rating_threshold, 
	  chi, wals, wals_l, wals_C,
	  als, chinmf, climf, 
	  mle_item, mle_user, canny, ctr);
  env_global = &env;
  
  // Reads the input files
  Ratings ratings(env);
  if (ratings.read(fname.c_str()) < 0) {
    fprintf(stderr, "error reading dataset from dir %s; quitting\n", 
	    fname.c_str());
    return -1;
  }

  if (rmse) {
    cout << "rmse" << endl;
    HGAPRec hgaprec(env, ratings);
    hgaprec.compute_rmse();
    exit(0);
  }
  
  if (chi) {
        cout << "chi" << endl;
    HGAPRec hgaprec(env, ratings);
    hgaprec.write_chi_training_matrix(wals_C);
    if (env.chinmf)
      hgaprec.run_chi_nmf();
    else if (env.als)
      hgaprec.run_chi_als();
    else if (env.wals)
      hgaprec.run_chi_wals(wals_l);
    else if (env.climf)
      hgaprec.run_chi_climf();
    exit(0);
  }
  
  if (test) {
        cout << "test" << endl;
    HGAPRec hgaprec(env, ratings);
    hgaprec.test();
    exit(0);
  }

  if (msr) {
        cout << "msr" << endl;
    HGAPRec hgaprec(env, ratings);
    hgaprec.gen_msr_csv();
    exit(0);
  }

  if (write_training) {
        cout << "write_training" << endl;
    HGAPRec hgaprec(env, ratings);
    if (lda)
      hgaprec.write_lda_training_matrix();
    else if (nmf)
      hgaprec.write_nmf_training_matrix();
    else if (vwlda)
      hgaprec.write_vwlda_training_matrix();
    exit(0);
  }

  if (nmfload) {
        cout << "nmfload" << endl;
    HGAPRec hgaprec(env, ratings);
    hgaprec.load_nmf_beta_and_theta();
    exit(0);
  }

  if (vwload) {
        cout << "vwload" << endl;
    HGAPRec hgaprec(env, ratings);
    hgaprec.load_vwlda_beta_and_theta();
    exit(0);
  }

  if (nmf) {
        cout << "nmf" << endl;
    HGAPRec hgaprec(env, ratings);
#ifdef NMFLIB    
    hgaprec.nmf();
#else
    printf("run configure with --enable-nmflib to use the -nmf option!\n");
    fflush(stdout);
#endif
    //hgaprec.write_training_matrix();
    //hgaprec.load_nmf_beta_and_theta();
    exit(0);
  }

  if (lda) {
        cout << "lda" << endl;
    HGAPRec hgaprec(env, ratings);
    hgaprec.load_lda_beta_and_theta();
    exit(0);
  }
  
  if (vwlda) {
        cout << "vwla" << endl;
    HGAPRec hgaprec(env, ratings);
    hgaprec.write_vwlda_training_matrix();
    hgaprec.run_vwlda();
    exit(0);
  }

  // if (ctr) {
  //   HGAPRec hgaprec(env, ratings);
  //   hgaprec.load_ctr_beta_and_theta();
  //   exit(0);
  // }
  
  if (gen_ranking_for_users) {
        cout << "gen_ranking_for_users" << endl;
    HGAPRec hgaprec(env, ratings);
    hgaprec.gen_ranking_for_users(true);
    exit(0);
  }
  
  // Runs the main algorithm
  if (batch) {
    HGAPRec hgaprec(env, ratings);
    
    if (bias && !hier)
      hgaprec.vb_bias();
    else if (hier)
      hgaprec.vb_hier();
    else if (mle_item)
      hgaprec.vb_mle_item();
    else if (mle_user)
      hgaprec.vb_mle_user();
    else if (canny)
      hgaprec.vb_canny();
    else
      hgaprec.vb();
  } else {
    printf("Quitting. Online inference not implemented.\n");
    fflush(stdout);
    exit(0);
  }
}*/
