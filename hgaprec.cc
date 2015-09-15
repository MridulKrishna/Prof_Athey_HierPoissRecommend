#include "hgaprec.hh"
#include "env.hh"
#include <iostream>
#include <iomanip>

#ifdef HAVE_NMFLIB
#include "./nmflib/include/common.h"
#include "./nmflib/include/nmfdriver.h"
#include <math.h>
#endif

// The constructor saves the environment, ratings, and parameters. It also runs the constructors for the random variables.
HGAPRec::HGAPRec(Env &env, Ratings &ratings)
  : _env(env), _ratings(ratings),
_n(env.n), _m(env.m), _uc(env.uc),_ic(env.ic),_k(env.k),_offset(env.offset),
//_validation_map(),
//_test_map(),
//_validation_users_of_movie(),
//_leave_one_out(),
_iter(0),
_start_time(time(0)),
_theta("theta", 0.3, 0.3, _n,_k,&_r),
_beta("beta", 0.3, 0.3, _m,_k,&_r),
_thetabias("thetabias", 0.3, 0.3, _n, 1, &_r),
_betabias("betabias", 0.3, 0.3, _m, 1, &_r),
_htheta("htheta", env.a, env.bp, _n, _k, &_r),
_hbeta("hbeta", env.c, env.dp, _m, _k, &_r),
_hsigma("hsigma", env.e, env.e*env.bp/(env.c*env.a), _ratings._itemObsScale , _n, _ic, &_r),
_hrho("hrho", env.f, env.f*env.dp/(env.c*env.a), _ratings._userObsScale , _m, _uc, &_r),
_thetarate("thetarate", env.ap, env.ap/env.bp, _n, &_r),
_betarate("betarate", env.cp, env.cp/env.dp, _m, &_r),
//phi(_k+_ic+_uc),
_theta_mle(_n, _k),
_beta_mle(_m, _k),
_old_theta_mle(_n, _k),
_old_beta_mle(_m, _k),
_lda_gamma(NULL), _lda_beta(NULL),
_nmf_theta(NULL), _nmf_beta(NULL),
_ctr_theta(NULL), _ctr_beta(NULL),
_prev_h(.0), _nh(.0),
_save_ranking_file(false),
_use_rate_as_score(true),
_topN_by_user(100),
_maxval(0), _minval(65536)
{
//  cout << env.dp << " " << env.bp << endl;
//  cout << "Offset: " << _offset << endl;
  // Initializes the random number generator
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  _r = gsl_rng_alloc(T);
  if (_env.seed)
    gsl_rng_set(_r, _env.seed);
  Env::plog("infer n:", _n);

  // Creates various output files it will use later
  
  string name = _env.outfname+"/"+_env.prefix +"/heldout.txt";
  _hf = fopen(name.c_str(), "w");
  if (!_hf)  {
    printf("cannot open heldout file:%s\n",  strerror(errno));
    exit(-1);
  }
  name = _env.outfname+"/"+_env.prefix +"/validation.txt";
  _vf = fopen(name.c_str(), "w");
  if (!_vf)  {
    printf("cannot open validation file:%s\n",  strerror(errno));
    exit(-1);
  }
  name = _env.outfname+"/"+_env.prefix +"/test.txt";
  _tf = fopen(name.c_str(), "w");
  if (!_tf)  {
    printf("cannot open test file:%s\n",  strerror(errno));
    exit(-1);
  }
  name = _env.outfname+"/"+_env.prefix +"/logl.txt";
  cout << name << endl;
  _af = fopen(name.c_str(), "w");
  if (!_af)  {
    printf("cannot open logl file:%s\n",  strerror(errno));
    exit(-1);
  }
  name = _env.outfname+"/"+_env.prefix +"/precision.txt";
  _pf = fopen(name.c_str(), "w");
  if (!_pf)  {
    printf("cannot open precision file:%s\n",  strerror(errno));
    exit(-1);
  }
  name = _env.outfname+"/"+_env.prefix +"/ndcg.txt";
  _df = fopen(name.c_str(), "w");
  if (!_df)  {
    printf("cannot open ndcg file:%s\n",  strerror(errno));
    exit(-1);
  }
  name = _env.outfname+"/"+_env.prefix +"/rmse.txt";
  _rf = fopen(name.c_str(), "w");
  if (!_rf)  {
    printf("cannot open rmse file:%s\n",  strerror(errno));
    exit(-1);
  }  

//  // Loads the other two files
//  if (!_env.write_training) {
////    cout << "aqui" << endl;
//    load_validation_and_test_sets();
//  }

  if (!_env.hier) {
    Env::plog("theta shape:", _theta.sprior());
    Env::plog("theta rate:", _theta.rprior());
    Env::plog("beta shape:", _beta.sprior());
    Env::plog("beta rate:", _beta.rprior());
  } else {
//    Env::plog("htheta shape:", _htheta.sprior());
//    Env::plog("htheta rate:", _htheta.rprior());
//
//    Env::plog("hbeta shape:", _hbeta.sprior());
//    Env::plog("hbeta rate:", _hbeta.rprior());
//
//    Env::plog("thetarate shape:", _thetarate.sprior());
//    Env::plog("thetarate rate:", _thetarate.rprior());
//
//    Env::plog("betarate shape:", _betarate.sprior());
//    Env::plog("betarate rate:", _thetarate.rprior());
  }
}

HGAPRec::~HGAPRec()
{
  fclose(_hf);
  fclose(_vf);
  fclose(_af);
  fclose(_pf);
  fclose(_tf);
  fclose(_rf);
}

// Deprecated. Now loads these in class Ratings

//void
//HGAPRec::load_validation_and_test_sets()
//{
//  char buf[4096];
//  sprintf(buf, "%s/validation.tsv", _env.datfname.c_str());
//  FILE *validf = fopen(buf, "r");
//  assert(validf);
//  if (_env.dataset == Env::NYT)
//    _ratings.read_nyt_train(validf, &_validation_map);
//  else
//    // Saves the ratings from the validation file in validation_map
//    _ratings.read_generic(validf, &_validation_map);
//  fclose(validf);
//
//  // Loop with iterator on the elements saved in _validation_map
//  // Builds a histogram with the number of users per movie in _validation_users_of_movie
//  for (CountMap::const_iterator i = _validation_map.begin();
//       i != _validation_map.end(); ++i) {
//    const Rating &r = i->first;
//    const Rating r2 = r;
//    _validation_users_of_movie[r.second]++;
//  }
//
//  sprintf(buf, "%s/test.tsv", _env.datfname.c_str());
//  FILE *testf = fopen(buf, "r");
//  assert(testf);
//  if (_env.dataset == Env::NYT)
//    _ratings.read_nyt_train(testf, &_test_map);
//  else
//    // Saves the ratings from the validation file in test_map
//    _ratings.read_generic(testf, &_test_map);
//  fclose(testf);
//  
//  // XXX: keeps one heldout test item for each user
//  // assumes leave-one-out
//  // JCC: Loop with iterator on the elements saved in _validation_map
//  // Builds a histogram with the number of users per movie in _validation_users_of_movie
//  for (CountMap::const_iterator i = _test_map.begin();
//       i != _test_map.end(); ++i) {
//    const Rating &r = i->first;
//    _leave_one_out[r.first] = r.second;
//    debug("adding %d -> %d to leave one out", r.first, r.second);
//  }
//
//  printf("+ loaded validation and test sets from %s\n", _env.datfname.c_str());
//  fflush(stdout);
//  Env::plog("test ratings", _test_map.size());
//  Env::plog("validation ratings", _validation_map.size());
//}

void
HGAPRec::initialize() {
  // Initializes the xi and eta parameters
  _thetarate.initialize2(0,_offset);
//  _thetarate.initialize2(_k*_env.a+_ic*_env.e,_offset);
  _thetarate.compute_expectations();
//  _thetarate.shape_curr().print();
//    _thetarate.rate_curr().print();
  
  _betarate.initialize2(0,_offset);
//  _betarate.initialize2(_k*_env.c+_uc*_env.f,_offset);
  _betarate.compute_expectations();
//  _betarate.shape_curr().print();
//  _betarate.rate_curr().print();
//  _betarate.expected_logv().print();

  
  // Initializes the beta and theta vectors
  _hbeta.initialize(_offset);
  _hbeta.initialize_exp(_offset);
//      _hbeta.shape_curr().print();
//        _hbeta.rate_curr().print();
  
  _htheta.initialize(_offset);
  _htheta.initialize_exp(_offset);
//  _htheta.shape_curr().print();
//  _htheta.rate_curr().print();
  
  _hsigma.initialize(_offset);
  _hsigma.initialize_exp(_offset);
//  _hsigma.shape_curr().print();
//  _hsigma.rate_curr().print();
  
  _hrho.initialize(_offset);
  _hrho.initialize_exp(_offset);
//  _hrho.shape_curr().print();
//  _hrho.rate_curr().print();
}

// Calculates the vector of probabilites for the multinomial distribution and saves it in argument phi
void
HGAPRec::get_phi(GPBase<Matrix> &theta, uint32_t ai, GPBase<Matrix> &beta, uint32_t bi, Array &phi)
{
  // Checks that the sizes of beta, theta, and phi agree
  assert (phi.size() == theta.k() &&
          phi.size() == beta.k());
  // Checks that the position does not exceed the array dimensions
  assert (ai < theta.n() && bi < beta.n());
  
  // Gets the expected log of theta and beta
  const double  **elogtheta = theta.expected_logv().const_data();
  const double  **elogbeta = beta.expected_logv().const_data();
  
  // Makes phi a zero vector
  phi.zero();
  
  // Adds each one of the elements of phi
  for (uint32_t k = 0; k < _k; ++k)
    phi[k] = elogtheta[ai][k] + elogbeta[bi][k];
  
  // Normalizes phi so it adds up to one
  phi.lognormalize();
}

// Calculates the vector of probabilites for the multinomial distribution and saves it in argument phi
void
HGAPRec::get_phi(GPBase<Matrix> &theta, uint32_t u, GPBase<Matrix> &beta, uint32_t i, GPBase<Matrix> &sigma, GPBase<Matrix> &rho, uint32_t ic, uint32_t uc, Array &phi)
{
  // Checks that the sizes of beta, theta, and phi agree
  assert (phi.size() == theta.k()+sigma.k()+rho.k() &&
	  phi.size() == beta.k()+sigma.k()+rho.k());
  // Checks that the position does not exceed the array dimensions
  assert (u < theta.n() && i < beta.n());
  
  // Gets the expected log of theta, beta, sigma, and rho
  const double  **elogtheta = theta.expected_logv().const_data();
  const double  **elogbeta = beta.expected_logv().const_data();
  const double  **elogsigma = sigma.expected_logv().const_data();
  const double  **elogrho = rho.expected_logv().const_data();
  
  const Matrix* userChar = &_ratings._userObs;
  const Matrix* itemChar = &_ratings._itemObs;
  
  // Makes phi a zero vector
  phi.zero();
  
  // Adds each one of the elements of phi
  for (uint32_t k = 0; k < _k; ++k) {
    phi[k] = elogtheta[u][k] + elogbeta[i][k];
//    cout << "phik " << phi[k] << endl;
  }
  for (uint32_t l = 0; l < _ic; ++l) {
    phi[_k+l] = elogsigma[u][l] + log(itemChar->get(i,l));
//    cout << "phil " << phi[_k+l] << endl;
  }
  for (uint32_t m = 0; m < _uc; ++m) {
    phi[_k+_ic+m] =  log(userChar->get(u,m)) + elogrho[i][m];
//    cout << "phim " << phi[_k+_ic+m] << endl;
//    cout << _k+_ic+m << " " << log(userChar->get(u,m)) << " " << elogrho[i][m] << endl;
  }
  
  // Normalizes phi so it adds up to one
  phi.lognormalize();
}

// Calculates the vector of probabilites for the multinomial distribution and saves it in argument phi. Takes into account that the item and user observables have to be scaled down by popularity and activity.
void
HGAPRec::get_phi(GPMatrix &theta, uint32_t u, GPMatrix &beta, uint32_t i, GPMatrix &sigma, GPMatrix &rho, GPArray &xi, GPArray &eta, uint32_t ic, uint32_t uc, Array &phi)
{
  // Checks that the sizes of beta, theta, and phi agree
  assert (phi.size() == theta.k()+sigma.k()+rho.k() &&
          phi.size() == beta.k()+sigma.k()+rho.k());
  // Checks that the position does not exceed the array dimensions
  assert (u < theta.n() && i < beta.n());
  
  // Gets the expected log of theta, beta, sigma, and rho
  const double  **elogtheta = theta.expected_logv().const_data();
  const double  **elogbeta = beta.expected_logv().const_data();
  const double  **elogsigma = sigma.expected_logv().const_data();
  const double  **elogrho = rho.expected_logv().const_data();
  const double  *elogxi = xi.expected_logv().const_data();
  const double  *elogeta = eta.expected_logv().const_data();
  
  const Matrix* userChar = &_ratings._userObs;
  const Matrix* itemChar = &_ratings._itemObs;
  
  // Makes phi a zero vector
  phi.zero();
  
//  _hsigma.shape_curr().print();
//  _hsigma.rate_curr().print();
  
  // Adds each one of the elements of phi
  for (uint32_t k = 0; k < _k; ++k) {
    phi[k] = elogtheta[u][k] + elogbeta[i][k];
//        cout << "phik " << phi[k] << endl;
  }
  for (uint32_t l = 0; l < _ic; ++l) {
    phi[_k+l] = elogsigma[u][l] -elogeta[i] + log(itemChar->get(i,l));
//        cout << "phil " << phi[_k+l] << endl;
  }
  for (uint32_t m = 0; m < _uc; ++m) {
    phi[_k+_ic+m] =  elogrho[i][m] -elogxi[u] + log(userChar->get(u,m));
//        cout << "phim " << phi[_k+_ic+m] << endl;
//        cout << _k+_ic+m << " " << log(userChar->get(u,m)) << " " << elogrho[i][m] << " " << elogxi[u] << endl;
  }
  
  // Normalizes phi so it adds up to one
  phi.lognormalize();
}

// Calculates the vector of probabilites for the multinomial distribution and saves it in argument phi. Only takes into account factors that are purely latent.
void
HGAPRec::get_phi(GPMatrix &theta, uint32_t u, GPMatrix &beta, uint32_t i, Array &phi)
{
  // Checks that the sizes of beta, theta, and phi agree
  assert (phi.size() == theta.k() &&
          phi.size() == beta.k() );
  // Checks that the position does not exceed the array dimensions
  assert (u < theta.n() && i < beta.n());
  
  // Gets the expected log of theta, beta, sigma, and rho
  const double  **elogtheta = theta.expected_logv().const_data();
  const double  **elogbeta = beta.expected_logv().const_data();
  
  // Makes phi a zero vector
  phi.zero();
  
  //  _hsigma.shape_curr().print();
  //  _hsigma.rate_curr().print();
  
  // Adds each one of the elements of phi
  for (uint32_t k = 0; k < _k; ++k) {
    phi[k] = elogtheta[u][k] + elogbeta[i][k];
    //        cout << "phik " << phi[k] << endl;
  }
  
  // Normalizes phi so it adds up to one
  phi.lognormalize();
}

void
HGAPRec::get_phi(GPBase<Matrix> &a, uint32_t ai, 
		 GPBase<Matrix> &b, uint32_t bi, 
		 double biasa, double biasb,
		 Array &phi)
{
  assert (phi.size() == a.k() + 2 &&
	  phi.size() == b.k() + 2);
  assert (ai < a.n() && bi < b.n());
  const double  **eloga = a.expected_logv().const_data();
  const double  **elogb = b.expected_logv().const_data();
  phi.zero();
  for (uint32_t k = 0; k < _k; ++k)
    phi[k] = eloga[ai][k] + elogb[bi][k];
  phi[_k] = biasa;
  phi[_k+1] = biasb;
  phi.lognormalize();
}

void
HGAPRec::get_phi(Matrix &a, uint32_t ai, 
		 GPBase<Matrix> &b, uint32_t bi, 
		 Array &phi)
{
  assert (phi.size() == a.n() &&
	  phi.size() == b.k());
  assert (ai < a.m() && bi < b.n());
  const double  **elogb = b.expected_logv().const_data();
  const double  **ad = a.const_data();
  phi.zero();
  for (uint32_t k = 0; k < _k; ++k)
    phi[k] = log(ad[ai][k]) + elogb[bi][k];
  phi.lognormalize();
}


void
HGAPRec::get_phi(GPBase<Matrix> &a, uint32_t ai, 
		 Matrix &b, uint32_t bi, 
		 Array &phi)
{
  assert (phi.size() == a.k() &&
	  phi.size() == b.n());
  const double  **eloga = a.expected_logv().const_data();
  const double  **bd = b.const_data();
  phi.zero();
  for (uint32_t k = 0; k < _k; ++k)
    phi[k] = log(bd[bi][k]) + eloga[ai][k];
  phi.lognormalize();
}


void
HGAPRec::write_lda_training_matrix()
{
  uint32_t x = 0;
  FILE *f = fopen(Env::file_str("/ldatrain.tsv").c_str(), "w");
  lerr("n = %d", _n);
  for (uint32_t n = 0; n < _n; ++n) {
    const vector<uint32_t> *movies = _ratings.get_movies(n);
    IDMap::const_iterator it = _ratings.seq2user().find(n);
    if (!movies || movies->size() == 0) {
      lerr("0 movies for user %d (%d)", n, it->second);
      continue;
    }
    
    uint32_t c = 0;
    for (uint32_t m = 0; m < _m; ++m)
      if (_ratings.r(n,m) == 0)
	c++;
    
    if (c == _m) {// null row
      lerr("null row found for user %d (%d)", n, it->second);
      continue;
    }

    x++;
    fprintf(f, "%d ", movies->size());
    
    for (uint32_t j = 0; j < movies->size(); ++j) {
      uint32_t m = (*movies)[j];
      yval_t y = _ratings.r(n,m);

      fprintf(f, " %d:%d", m, y);
    }
    fprintf(f, "\n");
  }
  fclose(f);
  lerr("wrote %d lines", x);
}

void
HGAPRec::write_chi_training_matrix(double wals_C)
{
  ostringstream sa_t, sa_v;
  uint32_t x = 0;
  FILE *f = fopen(Env::file_str("/chifull.tsv").c_str(), "w");
  FILE *g = fopen(Env::file_str("/chitrain.tsv").c_str(), "w");
  FILE *h = fopen(Env::file_str("/chivalidation.tsv").c_str(), "w");

  // matrix market format
  fprintf(f, "%%%%MatrixMarket matrix coordinate real general\n");
  //fprintf(f, "%% Generated 27-Feb-2012\n");
  fprintf(g, "%%%%MatrixMarket matrix coordinate real general\n");
  //fprintf(g, "%% Generated 27-Feb-2012\n");
  fprintf(h, "%%%%MatrixMarket matrix coordinate real general\n");
  //fprintf(h, "%% Generated 27-Feb-2012\n");


  BoolMap nusers_t, nusers_v, nusers;
  BoolMap nitems_t, nitems_v, nitems;
  uint32_t nratings_t = 0, nratings_v = 0;
  for (uint32_t n = 0; n < _n; ++n) {
    const vector<uint32_t> *movies = _ratings.get_movies(n);
    IDMap::const_iterator it = _ratings.seq2user().find(n);
    if (!movies || movies->size() == 0) {
      lerr("0 movies for user %d (%d)", n, it->second);
      continue;
    }
    uint32_t c = 0;
    for (uint32_t m = 0; m < _m; ++m)
      if (_ratings.r(n,m) == 0)
	c++;
    
    if (c == _m) {// null row
      lerr("null row found for user %d (%d)", n, it->second);
      continue;
    }

    nusers[n] = true;
    nusers_t[n] = true;
    for (uint32_t j = 0; j < movies->size(); ++j)  {
      uint32_t m = (*movies)[j];
      yval_t y = _ratings.r(n,m);
      nitems[m] = true;
      nitems_t[m] = true;
      nratings_t++;
      //
      // note: p_ui is binary
      // while r_ui determines the confidence c_ui as (1 + wals_C * y)
      //
      if (_env.wals) {
	sa_t << n+1 << " " << m+1 << " " << (int)(1 + y * wals_C) << " ";
	sa_t << (int)(y > 0 ? 1: 0) << "\n";
      } else if (_env.chinmf || _env.als || _env.climf)
	sa_t << n+1 << " " << m+1 << " " << (int)y << "\n";

      if (y > _maxval)
	_maxval = y;
      if (y < _minval)
	_minval = y;
    }
  }

  CountMap *mp = &_ratings._validation_map;
  for (CountMap::const_iterator i = mp->begin();
       i != mp->end(); ++i) {
    const Rating &e = i->first;
    uint32_t n = e.first;
    uint32_t m = e.second;
    yval_t r = i->second;
    nusers[n]=true;
    nusers_v[n]=true;
    nitems[m]=true;
    nitems_v[m]=true;
    nratings_v++;
    if (_env.wals) {
      sa_v << n+1 << " " << m+1 << " " << (int)(1 + r * wals_C) << " ";
      sa_v << (int)(r > 0 ? 1: 0) << "\n";
    } else if (_env.chinmf || _env.als || _env.climf)
      sa_v << n+1 << " " << m+1 << " " << (int)r << "\n";

    if (r > _maxval)
      _maxval = r;
    if (r < _minval)
      _minval = r;
  }

  if (_minval == _maxval)
    _maxval += 1;

  fprintf(f, "%d\t%d\t%d\n", nusers.size(), nitems.size(), nratings_t+nratings_v);
  fprintf(g, "%d\t%d\t%d\n", nusers_t.size(), nitems_t.size(), nratings_t);
  fprintf(h, "%d\t%d\t%d\n", nusers_v.size(), nitems_v.size(), nratings_v);
  
  fprintf(f, "%s", sa_t.str().c_str());
  fprintf(f, "%s", sa_v.str().c_str());

  fprintf(g, "%s", sa_t.str().c_str());
  fprintf(h, "%s", sa_v.str().c_str());

  fclose(h);
  fclose(g);
  fclose(f);
}

void
HGAPRec::load_ctr_beta_and_theta()
{
  char buf[4096];
  _ctr_theta = new Matrix(_n, _k);
  _ctr_beta = new Matrix(_m, _k);

  lerr("loading final-U.dat");
  _ctr_theta->load("final-U.dat", 0, false, 0);
  lerr("loading final-V.dat");
  _ctr_beta->load("final-V.dat", 0, false, 0);

  FILE *g = fopen("user_map.dat", "r");
  uint32_t a, b;
  assert(g);
  while(!feof(g)) {
    if (fscanf(g, "%d,%d\n", &a, &b) < 0){ 
      printf("error: unexpected lines in file\n");
      fclose(g);
      exit(-1);
    }
    _ctr_user_to_idx[a] = b;
  }
  fclose(g);


  g = fopen("item_map.dat", "r");
  assert(g);
  while(!feof(g)) {
    if (fscanf(g, "%d,%d\n", &a, &b) < 0){ 
      printf("error: unexpected lines in file\n");
      fclose(g);
      exit(-1);
    }
    _ctr_item_to_idx[a] = b;
  }
  fclose(g);

  FILE *f = fopen(Env::file_str("/user-map.csv").c_str(), "w");
  for (uint32_t n = 0; n < _n; ++n) {
    IDMap::const_iterator it = _ratings.seq2user().find(n);
    assert (it != _ratings.seq2user().end());
    fprintf(f, "%d,%d\n", it->second,_ctr_user_to_idx[it->second]);
  }
  fclose(f);

  f = fopen(Env::file_str("/item-map.csv").c_str(), "w");
  for (uint32_t n = 0; n < _m; ++n) {
    IDMap::const_iterator it = _ratings.seq2movie().find(n);
    assert (it != _ratings.seq2movie().end());
    fprintf(f, "%d,%d\n", it->second,_ctr_item_to_idx[it->second]);
  }
  fclose(f);
  
  lerr("saving ctr_theta.tsv");
  _ctr_theta->save(Env::file_str("/ctr_theta.tsv").c_str(), _ratings.seq2user());
  lerr("saving ctr_beta.tsv");
  _ctr_beta->save(Env::file_str("/ctr_beta.tsv").c_str(), _ratings.seq2movie());
}

void
HGAPRec::load_chi_beta_and_theta()
{
  char buf[4096];
  _chi_gamma = new Matrix(_n, _k);
  _chi_beta = new Matrix(_m, _k);

  if (_env.bias) {
    _chi_ubias = new Matrix(_n, 1);
    _chi_vbias = new Matrix(_m, 1);
    _chi_global = new Matrix(1, 1);
    lerr("u bias");
    _chi_ubias->mm_load_rowmajor(Env::file_str("/chitrain.tsv_U_bias.mm").c_str());
    lerr("v bias");
    _chi_vbias->mm_load_rowmajor(Env::file_str("/chitrain.tsv_V_bias.mm").c_str());
    lerr("global");
    _chi_global->mm_load_rowmajor(Env::file_str("/chitrain.tsv_global_mean.mm").c_str());
  }

  if (_env.wals || _env.als || _env.climf) {
    lerr("loading gamma from %s", Env::file_str("/chitrain.tsv_U.mm").c_str());
    _chi_gamma->mm_load_rowmajor(Env::file_str("/chitrain.tsv_U.mm").c_str());

    lerr("loading beta from %s", Env::file_str("/chitrain.tsv_V.mm").c_str());
    _chi_beta->mm_load_rowmajor(Env::file_str("/chitrain.tsv_V.mm").c_str());
  } else if (_env.chinmf) {
    lerr("loading gamma from %s", Env::file_str("/chifull.tsv_U.mm").c_str());
    _chi_gamma->mm_load_rowmajor(Env::file_str("/chifull.tsv_U.mm").c_str());
    
    lerr("loading beta from %s", Env::file_str("/chifull.tsv_V.mm").c_str());
    _chi_beta->mm_load_rowmajor(Env::file_str("/chifull.tsv_V.mm").c_str());
  } 

  _chi_beta->save(Env::file_str("/chi_beta.tsv").c_str(), _ratings.seq2movie());
  _chi_gamma->save(Env::file_str("/chi_gamma.tsv").c_str(), _ratings.seq2user());

  sprintf(buf, "%s/test_users.tsv", _env.datfname.c_str());
  FILE *f = fopen(buf, "r");
  if (!f) { 
    lerr("cannot open %s", buf);
    return;
  }
  _sampled_users.clear();
  _ratings.read_test_users(f, &_sampled_users);
  fclose(f);
  compute_precision(true);
  printf("DONE writing ranking.tsv in output directory\n");
  //compute_rmse();
  fflush(stdout);
  
  delete _chi_gamma;
  delete _chi_beta;
}


void
HGAPRec::run_chi_als()
{
  char chicmd[4096];
  sprintf(chicmd, "cd %s; GRAPHCHI_ROOT=/scratch/pgopalan/graphchi-cpp " \
	  "/scratch/pgopalan/graphchi-cpp/toolkits/collaborative_filtering/als " \
	  "--training=%s --validation=%s --lambda=0.01 "			\
	  "--minval=1 --maxval=5 --max_iter=100 --quiet=1 --D=%d > %s 2>&1",
	  Env::file_str("").c_str(),
	  "chitrain.tsv", "chivalidation.tsv", _k, "chi.out");
  lerr("%s", chicmd);
  if (system(chicmd) < 0) {
    lerr("could not execute %s (%s)", chicmd, strerror(errno));
    exit(-1);
  }
  load_chi_beta_and_theta();
}

void
HGAPRec::run_chi_wals(double wals_l)
{
  char chicmd[4096];
  sprintf(chicmd, "cd %s; GRAPHCHI_ROOT=/scratch/pgopalan/graphchi-cpp " \
	  "/scratch/pgopalan/graphchi-cpp/toolkits/collaborative_filtering/wals " \
	  "--training=%s --validation=%s --lambda=0.01 "		\
	  "--implicitratingtype=1 --implicitratingweight=1 --implicitratingpercentage=1 " \
	  "--minval=0 --maxval=1  --quiet=1 --D=%d --max_iter=100 > %s 2>&1",
	  Env::file_str("").c_str(),
	  "chitrain.tsv", "chivalidation.tsv", wals_l, _k, "chi.out");
  lerr("%s", chicmd);
  if (system(chicmd) < 0) {
    lerr("could not execute %s (%s)", chicmd, strerror(errno));
    exit(-1);
  }
  load_chi_beta_and_theta();
}

void
HGAPRec::run_chi_climf()
{
  char chicmd[4096];
  sprintf(chicmd, "cd %s; GRAPHCHI_ROOT=/scratch/pgopalan/graphchi-cpp " \
	  "/scratch/pgopalan/graphchi-cpp/toolkits/collaborative_filtering/climf " \
	  "--training=%s --validation=%s --binary_relevance_thresh=4 "	\
	  "--quiet=1 --sgd_gamma=1e-6 --D=%d --max_iter=500 --sgd_step_dec=0.9999 "\
	  "--sgd_lambda=1e-6 > %s 2>&1",
	  Env::file_str("").c_str(),
	  "chitrain.tsv", "chivalidation.tsv", _k, "chi.out");
  lerr("%s", chicmd);
  if (system(chicmd) < 0) {
    lerr("could not execute %s (%s)", chicmd, strerror(errno));
    exit(-1);
  }
  load_chi_beta_and_theta();
}

void
HGAPRec::run_chi_sgd()
{
  char chicmd[4096];
  sprintf(chicmd, "cd %s; GRAPHCHI_ROOT=/scratch/pgopalan/graphchi-cpp "\
	  "/scratch/pgopalan/graphchi-cpp/toolkits/collaborative_filtering/biassgd " \
	  "--training=%s --validation=%s --biassgd_lambda=1e-3 --biassgd_gamma=1e-3 " \
	  "--minval=1 --maxval=5 --max_iter=100 --quiet=1 --D=%d > %s 2>&1",
	  Env::file_str("").c_str(),
	  "chitrain.tsv", "chivalidation.tsv", _k, "chi.out");
  lerr("%s", chicmd);
  if (system(chicmd) < 0) {
    lerr("could not execute %s (%s)", chicmd, strerror(errno));
    exit(-1);
  }
  load_chi_beta_and_theta();
}

void
HGAPRec::run_chi_pmf()
{
  char chicmd[4096];
  sprintf(chicmd, "cd %s; GRAPHCHI_ROOT=/scratch/pgopalan/graphchi-cpp "\
	  "/scratch/pgopalan/graphchi-cpp/toolkits/collaborative_filtering/pmf "\
	  "--training=%s --minval=1 --maxval=5 --max_iter=100 --pmf_burn_in=5 "\
	  "--allow_zeros=1 --quiet=1 --D=%d --matrixmarket=true --pmf_additional_output=1 "\
	  "> %s 2>&1",
	  Env::file_str("").c_str(),
	  "chitrain.tsv", _k, "chi.out");
  lerr("%s", chicmd);
  if (system(chicmd) < 0) {
    lerr("could not execute %s (%s)", chicmd, strerror(errno));
    exit(-1);
  }
  load_chi_beta_and_theta();
}

void
HGAPRec::run_chi_nmf()
{
  char chicmd[4096];
  sprintf(chicmd, "cd %s; GRAPHCHI_ROOT=/scratch/pgopalan/graphchi-cpp "\
	  "/scratch/pgopalan/graphchi-cpp/toolkits/collaborative_filtering/nmf "\
	  "--training=%s --minval=%d --maxval=%d --max_iter=500 --quiet=1 "\
	  "--D=%d> %s 2>&1",
	  Env::file_str("").c_str(),
	  "chifull.tsv", _minval, _maxval,
	  _k, "chi.out");
  lerr("%s", chicmd);
  if (system(chicmd) < 0) {
    lerr("could not execute %s (%s)", chicmd, strerror(errno));
    exit(-1);
  }
  load_chi_beta_and_theta();
}


void
HGAPRec::run_vwlda()
{
  char vwcmd[4096]; 
  sprintf(vwcmd, "/scratch/pgopalan/vowpal_wabbit/vowpalwabbit/vw  "	\
	  "--lda %d  --lda_alpha %f --lda_rho %f --lda_D %d "		\
	  "--minibatch %d --power_t %f --initial_t %d "			\
	  "%s -b %d -p %s --readable_model %s > %s 2>&1",
	  _k, (double)1.0/_k, 
	  (double)1.0/_k, 
	  _n, 
	  256, 0.5, 1,
	  Env::file_str("/ldatrain.tsv").c_str(), 
	  (int)(log2(_m)+1),
	  Env::file_str("/gamma.tsv").c_str(), 
	  Env::file_str("/beta.tsv").c_str(),
	  Env::file_str("/vwlda.out").c_str());
  lerr("%s", vwcmd);
  if (system(vwcmd) < 0) {
    lerr("could not execute %s (%s)", vwcmd, strerror(errno));
    exit(-1);
  }
  load_vwlda_beta_and_theta();
}

void
HGAPRec::write_vwlda_training_matrix()
{
  uint32_t x = 0;
  FILE *f = fopen(Env::file_str("/ldatrain.tsv").c_str(), "w");
  for (uint32_t n = 0; n < _n; ++n) {
    const vector<uint32_t> *movies = _ratings.get_movies(n);
    IDMap::const_iterator it = _ratings.seq2user().find(n);
    if (!movies || movies->size() == 0) {
      lerr("0 movies for user %d (%d)", n, it->second);
      continue;
    }
    
    uint32_t c = 0;
    for (uint32_t m = 0; m < _m; ++m)
      if (_ratings.r(n,m) == 0)
	c++;
    
    if (c == _m) {// null row
      lerr("null row found for user %d (%d)", n, it->second);
      continue;
    }

    x++;
    fprintf(f, "|");
    
    for (uint32_t j = 0; j < movies->size(); ++j) {
      uint32_t m = (*movies)[j];
      yval_t y = _ratings.r(n,m);

      fprintf(f, " %d:%d", m, y);
    }
    fprintf(f, "\n");
  }
  fclose(f);
  lerr("wrote %d lines", x);
}


void
HGAPRec::write_nmf_training_matrix()
{
  uint32_t nrows = 0;
  for (uint32_t n = 0; n < _n; ++n) {
    const vector<uint32_t> *movies = _ratings.get_movies(n);
    IDMap::const_iterator it = _ratings.seq2user().find(n);
    if (movies->size() == 0) {
      lerr("0 movies for user %d (%d)", n, it->second);
      continue;
    }
    
    uint32_t c = 0;
    for (uint32_t m = 0; m < _m; ++m)
      if (_ratings.r(n,m) == 0)
	c++;

    if (c == _m) {// null row
      lerr("null row found for user %d (%d)", n, it->second);
      continue;
    }
    nrows++;
  }

  FILE *f = fopen(Env::file_str("/trainm.tsv").c_str(), "w");
  // only for libNMF
  fprintf(f, "%d\n", nrows);
  fprintf(f, "%d\n", _m);
  for (uint32_t n = 0; n < _n; ++n) {
    const vector<uint32_t> *movies = _ratings.get_movies(n);
    IDMap::const_iterator it = _ratings.seq2user().find(n);
    if (movies->size() == 0) {
      lerr("0 movies for user %d (%d)", n, it->second);
      continue;
    }
    
    uint32_t c = 0;
    for (uint32_t m = 0; m < _m; ++m)
      if (_ratings.r(n,m) == 0)
	c++;

    if (c == _m) {// null row
      lerr("null row found for user %d (%d)", n, it->second);
      continue;
    }

    for (uint32_t m = 0; m < _m; ++m)
      fprintf(f, "%d\t", _ratings.r(n,m));
    fprintf(f, "\n");
  }
  fclose(f);
}

void
HGAPRec::load_lda_beta_and_theta()
{
  char buf[4096];
  _lda_gamma = new Matrix(_n, _k);
  _lda_beta = new Matrix(_k, _m);

  _lda_gamma->load("gamma.tsv", 0);
  _lda_beta->load("beta.tsv", 0);
  
  _lda_gamma->normalize1();
  _lda_beta->exp1();

  IDMap m;
  _lda_gamma->save(Env::file_str("/lda_gamma.tsv").c_str(), m);
  _lda_beta->save(Env::file_str("/lda_beta.tsv").c_str(), m);
  _lda_beta->save_transpose(Env::file_str("/lda_beta_transpose.tsv").c_str(), m);

  sprintf(buf, "%s/test_users.tsv", _env.datfname.c_str());
  FILE *f = fopen(buf, "r");
  if (!f) { 
    lerr("cannot open %s", buf);
    return;
  }
  _sampled_users.clear();
  _ratings.read_test_users(f, &_sampled_users);
  fclose(f);
  compute_precision(true);
  printf("DONE writing ranking.tsv in output directory\n");
  fflush(stdout);
  
  delete _lda_gamma;
  delete _lda_beta;
}

void
HGAPRec::load_vwlda_beta_and_theta()
{
  char buf[4096];
  _lda_gamma = new Matrix(_n, _k);
  _lda_beta = new Matrix(_k, _m);

  _lda_beta->load(Env::file_str("/beta.tsv").c_str(), 1, true, 11);
  _lda_gamma->load(Env::file_str("/gamma.tsv").c_str(), 0);
  
  _lda_gamma->normalize1();
  _lda_beta->normalize1();
  //_lda_beta->exp1();

  _lda_gamma->save(Env::file_str("/lda_gamma.tsv").c_str(), _ratings.seq2user());
  _lda_beta->save(Env::file_str("/lda_beta.tsv").c_str(), _ratings.seq2movie());
  _lda_beta->save_transpose(Env::file_str("/lda_beta_transpose.tsv").c_str(), _ratings.seq2movie());

  sprintf(buf, "%s/test_users.tsv", _env.datfname.c_str());
  FILE *f = fopen(buf, "r");
  if (!f) { 
    lerr("cannot open %s", buf);
    return;
  }
  _sampled_users.clear();
  _ratings.read_test_users(f, &_sampled_users);
  fclose(f);
  compute_precision(true);
  printf("DONE writing ranking.tsv in output directory\n");
  fflush(stdout);
  
  delete _lda_gamma;
  delete _lda_beta;
}

void
HGAPRec::load_nmf_beta_and_theta()
{
  _nmf_theta = new Matrix(_n, _k);
  _nmf_beta = new Matrix(_m, _k);

  char buf[4096];
  _nmf_theta->nmf_load(Env::file_str("/theta.tsv").c_str(), false);
  _nmf_beta->nmf_load(Env::file_str("/beta.tsv").c_str(), true);

  IDMap m;
  _nmf_theta->save(Env::file_str("/nmf_theta.tsv").c_str(), m);
  _nmf_beta->save(Env::file_str("/nmf_beta.tsv").c_str(), m);

  sprintf(buf, "%s/test_users.tsv", _env.datfname.c_str());
  FILE *f = fopen(buf, "r");
  if (!f) { 
    lerr("cannot open %s", buf);
    return;
  }
  _sampled_users.clear();
  _ratings.read_test_users(f, &_sampled_users);
  fclose(f);
  compute_precision(true);
  printf("DONE writing ranking.tsv in output directory\n");
  fflush(stdout);

  delete _nmf_theta;
  delete _nmf_beta;
}

/*
void
HGAPRec::load_nmf_beta_and_theta()
{
  char buf[4096];
  _nmf_theta.load("theta.tsv", 0);
  _nmf_beta.load("beta.tsv", 0);

  IDMap m;
  _nmf_theta.save(Env::file_str("/nmf_theta.tsv").c_str(), m);
  _nmf_beta.save(Env::file_str("/nmf_beta.tsv").c_str(), m);

  sprintf(buf, "%s/test_users.tsv", _env.datfname.c_str());
  FILE *f = fopen(buf, "r");
  if (!f) { 
    lerr("cannot open %s", buf);
    return;
  }
  _sampled_users.clear();
  _ratings.read_test_users(f, &_sampled_users);
  fclose(f);
  compute_precision(true);
  printf("DONE writing ranking.tsv in output directory\n");
  fflush(stdout);
}
*/

#ifdef HAVE_NMFLIB
void
HGAPRec::nmf()
{
  write_nmf_training_matrix();
  string s = Env::file_str("/trainm.tsv");
  string thetas = Env::file_str("/theta.tsv");
  string betas = Env::file_str("/beta.tsv");

  options_t opts;
  opts.rep = 1;
  opts.init = nndsvd;
  opts.min_init = 0;
  opts.max_init = 1;
  opts.w_out = thetas.c_str();
  opts.h_out = betas.c_str();
  opts.TolX = 2.0E-02;
  opts.TolFun = 2.0E-02;
  opts.nndsvd_maxiter = -1; //if set to -1 - default value will be set in generateMatrix
  opts.nndsvd_blocksize = 1;
  opts.nndsvd_tol = 2E-08;
  opts.nndsvd_ncv = -1;		

  nmfDriver(s.c_str(), _k, 100, NULL, NULL, mu, &opts);
  load_nmf_beta_and_theta();
}
#endif

// Main method for the hierarchical model (the one in the paper)
void
HGAPRec::vb_hier()
{
  // Initial values of the parameters, according to the posterior plus a random shock

  initialize();
//  _betarate.shape_curr().print();
//  _betarate.rate_curr().print();
//  _thetarate.shape_curr().print();
//  _thetarate.rate_curr().print();
//  _hbeta.shape_curr().print();
//  _hbeta.rate_curr().print();
//  _htheta.shape_curr().print();
//  _htheta.rate_curr().print();
//  _hsigma.shape_curr().print();
//  _hsigma.rate_curr().print();
//  _hrho.shape_curr().print();
//  _hrho.rate_curr().print();
  
//  _ratings._userObsScale.print();
//  _ratings._itemObsScale.print();
  
  cout << "Initialized" << endl;
//  cout << _env.reportfreq << endl;
  
  // Matrices that save the means of variables
  Matrix betaMeans(_k,_env.max_iterations+1);
  Matrix thetaMeans(_k,_env.max_iterations+1);
  Matrix sigmaMeans(_ic,_env.max_iterations+1);
  Matrix rhoMeans(_uc,_env.max_iterations+1);
  
  Array xiMeans(_env.max_iterations+1);
  Array etaMeans(_env.max_iterations+1);
  
  Array betaMean(_k);
  Array thetaMean(_k);
  Array sigmaMean(_ic);
  Array rhoMean(_uc);
  
  _hbeta.expected_means(betaMean);
  _htheta.expected_means(thetaMean);
  _hsigma.expected_means(sigmaMean);
  _hrho.expected_means(rhoMean);
  
  double xiMean = _betarate.expected_mean();
  double etaMean = _thetarate.expected_mean();
  
  betaMeans.set_col(0, betaMean);
  thetaMeans.set_col(0, thetaMean);
  sigmaMeans.set_col(0, sigmaMean);
  rhoMeans.set_col(0, rhoMean);
  
  xiMeans.set(0,xiMean);
  etaMeans.set(0,etaMean);
  
//  cout << "Eta" << endl;
//  _betarate.shape_curr().print();
//  _betarate.rate_curr().print();
//  cout << "Xi" << endl;
//  _thetarate.shape_curr().print();
//  _thetarate.rate_curr().print();
  
  //lerr("htheta = %s", _htheta.rate_next().s().c_str());

  ///////////
// Original code
//  // What is bias?
//  uint32_t x;
//  if (_env.bias)
//    x = _k+2;
//  else
//    x = _k;
  ///////////
  
  uint32_t x = _k+_uc+_ic;
//
  // Constructs the array for the parameters of the multinomial distribution
  Array phi(x);
  
  bool stop = false;
  
  while (!stop) {
    // Stop if the max number of iterations is reached
    if (_iter > _env.max_iterations) {
      exit(0);
    }
    
    // Loop over users
    for (uint32_t n = 0; n < _n; ++n) {
//      cout << "Loop over users " << n << endl;
      
      // Gets the matrix of items for each user and stores it in movies
      const vector<uint32_t> *movies = _ratings.get_movies(n);
      
      // Loop over each user's items
      for (uint32_t j = 0; movies && j < movies->size(); ++j) {
//        cout << "Loop over items " << j <<endl;
        // Gets the code of the movie
        uint32_t m = (*movies)[j];
        
//        //////// JCC: Tests
//        IDMap seq2m = _ratings.seq2movie();
//        IDMap seq2u = _ratings.seq2user();
//        
//        uint32_t username = seq2u.at(n);
//        uint32_t moviename = seq2m.at(j);
//        cout << "User: " << n << " " << username << endl;
//        cout << "Movie: " << j << " " << m << " " << moviename << endl;
//        ////////
        
        // Get the movie rating
        yval_t y = _ratings.r(n,m);

//        //////// JCC: Tests
//        uint32_t y2 = _ratings.r(n,m);
//        cout << _ratings.r(n,m) << endl;
//        cout << "Rating: " << y << endl;
//        cout << "Rating: " << y2 << endl;
//        ////////
        
        // Finds phi from the current parameters of hbeta, htheta, hsigma, and hrho (the equation in step 1 of the algorithm in the paper)
        get_phi(_htheta, n, _hbeta, m, _hsigma, _hrho, _thetarate, _betarate, _ic, _uc, phi);
        
//        get_phi(_htheta, n, _hbeta, m, _hsigma, _hrho, _ic, _uc, phi);
        
//        cout << phi.sum(0) << endl;
//        if (_iter < 2) {
//          cout << "Phi: " << endl;
//          phi.print();
//        }
        
//        _phi.set_row( r, phi);
      
        // Makes phi sum up to y to get y_{ui} phi_{uik}
        if (y > 1) {
          phi.scale(y);
        }
        

//        if (_iter == 0 && n == 0 && j == 1)
//        cout << n << " " << j << endl;
//          phi.print();
        
//        cout << phi.sum(0) << endl;
//        phi.print();
        
        // Defines the subarrays of phi for latent variables, user observables, and item observables and updates the next shape parameter of theta and beta (gamma and kappa) by adding y_{ui} phi_{uik} to the nth row of gamma and the mth row of kappa (the first equation in steps 2 and 3 of the algorithm in the paper)
        Array phik(_k);
        Array phil(_ic);
        Array phim(_uc);
        
        if (_k>0) {
          phik.copy_from(phi.subarray(0,_k-1));
          _htheta.update_shape_next1(n, phik);
          _hbeta.update_shape_next1(m, phik);
        }
        
        if (_ic > 0) {
          phil.copy_from(phi.subarray(_k,_k+_ic-1));
          _hsigma.update_shape_next1(n, phil);
        }
        
        if ( _uc > 0) {
          phim.copy_from(phi.subarray(_k+_ic,_k+_ic+_uc-1));
          _hrho.update_shape_next1(m, phim);
        }
        
        if (_env.bias) {
          _thetabias.update_shape_next3(n, 0, phi[_k]);
          _betabias.update_shape_next3(m, 0, phi[_k+1]);
        }
//        phi.print();
      }
    } // End of loop over users/movies
    
    debug("htheta = %s", _htheta.expected_v().s().c_str());
    debug("hbeta = %s", _hbeta.expected_v().s().c_str());
    
    //----------------------------------
    // Updates for user parameters
    //----------------------------------
    
    // If there are latent characteristics...
    
    if (_k>0) {
      Array betarowsum(_k);
      // Saves the sums over items of expected values for each factor ( the second part of \gamma^{rte}_{uk})
      _hbeta.sum_rows(betarowsum);
      
      // Sets the prior rate based on expectations with current parameters, i.e.,\frac{\kappa^{shp}}{\kappa^{rte}}
      _htheta.set_prior_rate(_thetarate.expected_v(),
                             _thetarate.expected_logv());
      debug("adding %s to theta rate", _thetarate.expected_v().s().c_str());
      debug("betarowsum %s", betarowsum.s().c_str());
      
      // Adds the previous sum (betarowsum) to \frac{\kappa^{shp}}{\kappa^{shp}} in the next rate
      _htheta.update_rate_next(betarowsum);
      // Swaps the current and the next values for the parameters
      _htheta.swap();
      
      // Computes expectations and log expectations based on the new parameters
      _htheta.compute_expectations();
    }
    
    // If there are observed item characteristics...
    if (_ic > 0) {
      // Computes the sums of user characteristics (second term for \mu_{ul}^{rte}
      Array itemSum(_ic);
      _ratings._itemObs.weighted_colsum(_betarate.expected_inv(),itemSum);
//      _ratings._itemObs.colsum(itemSum);
    
      // Sets the prior rate based on weighted expectations with current parameters, i.e.,x_l * e/(ca) * \frac{\kappa^{shp}}{\kappa^{rte}}
      _hsigma.set_prior_rate_scaled(_thetarate.expected_v(),_env.e/(_env.c*_env.a), _ratings._itemObsScale);
      
      // Adds the previous sum (itemSum) to \frac{\kappa^{shp}}{\kappa^{shp}} in the next rate
      _hsigma.update_rate_next(itemSum);
      // Swaps the current and the next values for the parameters
      _hsigma.swap();
      
      // Computes expectations and log expectations based on the new parameters
      _hsigma.compute_expectations();
    }
    
//    cout << "Theta: " << endl;
//    _htheta.shape_curr().print();
//    _htheta.rate_curr().print();
//
//    cout << "Sigma: " << endl;
//    _hsigma.shape_curr().print();
//    _hsigma.rate_curr().print();
    
    //----------------------------------
    // Updates for item parameters
    //----------------------------------
    
    // If there are latent variables...
    if (_k>0) {
      Array thetarowsum(_k);
      // Saves the sums of expected values over users for each factor (the second part of \lambda^{rte}_{ik})
      _htheta.sum_rows(thetarowsum);
      
      // Sets the prior rate based on expectations with current parameters, i.e., \frac{\tau^{shp}}{\tau^{rte}}
      _hbeta.set_prior_rate(_betarate.expected_v(),
                            _betarate.expected_logv());
      
      // Adds the previous sum to \frac{\tau^{shp}}{\tau^{shp}} in the next rate
      _hbeta.update_rate_next(thetarowsum);
      // Swaps the current and the next values for the parameters
      _hbeta.swap();
      
      // Computes expectations and log expectations based on the new parameters
      _hbeta.compute_expectations();
    }
    
    // If there are observed user characteristics...
    if (_ic > 0) {
      // Computes the sums of user characteristics (second term for \mu_{ul}^{rte}
      Array userSum(_uc);
      
//      _thetarate.shape_curr().print();
//      _thetarate.rate_curr().print();
//      _thetarate.expected_inv().print();
      _ratings._userObs.weighted_colsum(_thetarate.expected_inv(),userSum);
//      _ratings._userObs.colsum(userSum);
      
      // Sets the prior rate based on expectations with current parameters, i.e.,\frac{\tau^{shp}}{\tau^{rte}}
      _hrho.set_prior_rate_scaled(_betarate.expected_v(),_env.f/(_env.c*_env.a),_ratings._userObsScale);
      
      // Adds the previous sum to \frac{\kappa^{shp}}{\kappa^{shp}} in the next rate
      _hrho.update_rate_next(userSum);
      // Swaps the current and the next values for the parameters
      _hrho.swap();
      
      // Computes expectations and log expectations based on the new parameters
//      cout << "rho " << _iter << endl;
      _hrho.compute_expectations();
    }
    
//    cout << "Beta: " << endl;
//    _hbeta.shape_curr().print();
//    _hbeta.rate_curr().print();
//    
//    cout << "Rho: " << endl;
//    _hrho.shape_curr().print();
//    _hrho.rate_curr().print();
    
    //-------------------------------------
    // Update the rate parameters for the popularity and activity parameters
    //-------------------------------------

    // Adds Ka+Le to the shape parameter of xi
    _thetarate.update_shape_next(_k*_env.a+_ic*_env.e);
    
    // If there are latent variables...
    if (_k>0) {
      Array thetacolsum(_n);
      // Computes the second term of \kappa^{rte}_{uk})
      _htheta.sum_cols(thetacolsum);
      
      // Adds the new term to the rate parameter of xi
      _thetarate.update_rate_next(thetacolsum);
      
//      thetacolsum.print();
    }
    
    // With unobserved item characteristics...
    if (_ic > 0) {
      Array sigmacolsum(_n);
      // Computes the third term of \kappa^{rte}_{uk})
      _hsigma.sum_cols_weight(_ratings._itemObsScale,sigmacolsum);
//      _ratings._itemObsScale.print();
      // Adds the new term to the rate parameter of xi
//      _thetarate.update_rate_next(sigmacolsum.scale(_env.e/(_env.c*_env.a)));
      double scale = _env.e/(_env.c*_env.a);
      _thetarate.update_rate_next(sigmacolsum,scale);
//      sigmacolsum.print();
    }
    
    debug("thetacolsum = %s", thetacolsum.s().c_str());
    
    // Swaps the current and the next values for the parameters
    _thetarate.swap();
    // Computes the expectations with the (new) current values
    _thetarate.compute_expectations();
    
    
    // Adds Kc+Mf to the shape parameter of eta
    _betarate.update_shape_next(_k*_env.c+_uc*_env.f);
    
    // If there are latent variables...
    if (_k>0) {
      
      Array betacolsum(_m);
      
      // Computes the second term of \tau^{rte}_{uk})
      _hbeta.sum_cols(betacolsum);
      // Adds the new term to the rate parameter of eta
      _betarate.update_rate_next(betacolsum);
    }
    
    // With unobserved user characteristics...
    if (_uc > 0) {
      Array rhocolsum(_m);
      // Computes the third term of \tau^{rte}_{uk})
      _hrho.sum_cols_weight(_ratings._userObsScale,rhocolsum);
      // Adds the new term to the rate parameter of eta
      _betarate.update_rate_next(rhocolsum.scale(_env.f/(_env.c*_env.a)));
    }
    
    debug("betacolsum = %s", betacolsum.s().c_str());
    
    // Swaps the current and the next values for the parameters
    _betarate.swap();
    // Computes the expectations with the (new) current values
    _betarate.compute_expectations();
    
//    if (_iter == 59 || _iter == 0) {
//      cout << "Beta shape: " << endl;
//      _hbeta.shape_curr().print();
//      cout << "Beta rate: " << endl;
//      _hbeta.rate_curr().print();
//      cout << "Theta shape: " << endl;
//      _htheta.shape_curr().print();
//      cout << "Theta rate: " << endl;
//      _htheta.rate_curr().print();
//      
//      cout << "Eta shape: " << endl;
//      _betarate.shape_curr().print();
//      cout << "Eta rate: " << endl;
//      _betarate.rate_curr().print();
//      cout << "Xi shape: " << endl;
//      _thetarate.shape_curr().print();
//      cout << "Xi rate: " << endl;
//      _thetarate.rate_curr().print();
//    }
    
    printf("iteration %d\n", _iter);
    
    Array betaMean(_k);
    Array thetaMean(_k);
    Array sigmaMean(_ic);
    Array rhoMean(_uc);
    
    _hbeta.expected_means(betaMean);
    _htheta.expected_means(thetaMean);
    _hsigma.expected_means(sigmaMean);
    _hrho.expected_means(rhoMean);
    
    double xiMean = _betarate.expected_mean();
    double etaMean = _thetarate.expected_mean();
    
    betaMeans.set_col(_iter+1, betaMean);
    thetaMeans.set_col(_iter+1, thetaMean);
    sigmaMeans.set_col(_iter+1, sigmaMean);
    rhoMeans.set_col(_iter+1, rhoMean);
    
    xiMeans.set(_iter+1,xiMean);
    etaMeans.set(_iter+1,etaMean);
    
//    _betarate.expected_v().print();
//    _thetarate.expected_v().print();
    
    fflush(stdout);
    if (_iter % _env.reportfreq == 0) {
      
      compute_likelihood(false);
      stop = compute_likelihood(true);
      //compute_rmse();
      save_model();
      // Computes and saves number of relevant recommendations among best ranked items
      compute_precision(false);
      // Computes and saves average ranking of items in test set
      compute_itemrank(false);
      //gen_ranking_for_users(false);
      if (_env.logl)
        logl();
      
      string nameBeta = string("/betaMeans.tsv");
      string nameTheta = string("/thetaMeans.tsv");
      string nameSigma = string("/sigmaMeans.tsv");
      string nameRho = string("/rhoMeans.tsv");
      string nameXi = string("/xiMeans.tsv");
      string nameEta = string("/etaMeans.tsv");
      
      betaMeans.save(_env.outfname+"/"+Env::outfile_str(nameBeta));
      thetaMeans.save(_env.outfname+"/"+Env::outfile_str(nameTheta));
      sigmaMeans.save(_env.outfname+"/"+Env::outfile_str(nameSigma));
      rhoMeans.save(_env.outfname+"/"+Env::outfile_str(nameRho));
      
      xiMeans.save(_env.outfname+"/"+Env::outfile_str(nameXi));
      etaMeans.save(_env.outfname+"/"+Env::outfile_str(nameEta));
    }
    
    if (stop) {
      do_on_stop();
//      _thetarate.expected_inv().print();
    } else {
      // Saves the matrices in files
      if (_env.save_state_now) {
        lerr("Saving state at iteration %d duration %d secs", _iter, duration());
        do_on_stop();
      }
      _iter++;

    }
//    if ( _iter == 55) {
//      _betarate.shape_curr().print();
//      _betarate.rate_curr().print();
//    }
  }
}

// Main method for the hierarchical model (the one in the paper)
void
HGAPRec::vb_hier_latents_first()
{
  // Initial values of the parameters, according to the posterior plus a random shock
  
  initialize();
  //  _betarate.shape_curr().print();
  //  _betarate.rate_curr().print();
  //  _thetarate.shape_curr().print();
  //  _thetarate.rate_curr().print();
  //  _hbeta.shape_curr().print();
  //  _hbeta.rate_curr().print();
  //  _htheta.shape_curr().print();
  //  _htheta.rate_curr().print();
  //  _hsigma.shape_curr().print();
  //  _hsigma.rate_curr().print();
  //  _hrho.shape_curr().print();
  //  _hrho.rate_curr().print();
  
  //  _ratings._userObsScale.print();
  //  _ratings._itemObsScale.print();
  
  cout << "Initialized" << endl;
  //  cout << _env.reportfreq << endl;
  
  // Matrices that save the means of variables
  Matrix betaMeans(_k,_env.max_iterations+1);
  Matrix thetaMeans(_k,_env.max_iterations+1);
  Matrix sigmaMeans(_ic,_env.max_iterations+1);
  Matrix rhoMeans(_uc,_env.max_iterations+1);
  
  Array xiMeans(_env.max_iterations+1);
  Array etaMeans(_env.max_iterations+1);
  
  Array betaMean(_k);
  Array thetaMean(_k);
  Array sigmaMean(_ic);
  Array rhoMean(_uc);
  
  _hbeta.expected_means(betaMean);
  _htheta.expected_means(thetaMean);
  _hsigma.expected_means(sigmaMean);
  _hrho.expected_means(rhoMean);
  
  double xiMean = _betarate.expected_mean();
  double etaMean = _thetarate.expected_mean();
  
  betaMeans.set_col(0, betaMean);
  thetaMeans.set_col(0, thetaMean);
  sigmaMeans.set_col(0, sigmaMean);
  rhoMeans.set_col(0, rhoMean);
  
  xiMeans.set(0,xiMean);
  etaMeans.set(0,etaMean);
  
  //  cout << "Eta" << endl;
  //  _betarate.shape_curr().print();
  //  _betarate.rate_curr().print();
  //  cout << "Xi" << endl;
  //  _thetarate.shape_curr().print();
  //  _thetarate.rate_curr().print();
  
  // Constructs the array for the parameters of the multinomial distribution
  Array phiLatents(_k);
  
  while (_iter < 100) {
    // Stop if the max number of iterations is reached
    if (_iter > _env.max_iterations) {
      exit(0);
    }
    
    // Loop over users
    for (uint32_t n = 0; n < _n; ++n) {
      //      cout << "Loop over users " << n << endl;
      
      // Gets the matrix of items for each user and stores it in movies
      const vector<uint32_t> *movies = _ratings.get_movies(n);
      
      // Loop over each user's items
      for (uint32_t j = 0; movies && j < movies->size(); ++j) {
        //        cout << "Loop over items " << j <<endl;
        // Gets the code of the movie
        uint32_t m = (*movies)[j];
        
        // Get the movie rating
        yval_t y = _ratings.r(n,m);
       
        // Finds phi from the current parameters of hbeta, htheta (the equation in step 1 of the algorithm in the paper)
        get_phi(_htheta, n, _hbeta, m, phiLatents);
        
        //        get_phi(_htheta, n, _hbeta, m, _hsigma, _hrho, _ic, _uc, phi);
        
        //        cout << phi.sum(0) << endl;
        //        if (_iter < 2) {
        //          cout << "Phi: " << endl;
        //          phi.print();
        //        }
        
        //        _phi.set_row( r, phi);
        
        // Makes phi sum up to y to get y_{ui} phi_{uik}
        if (y > 1) {
          phiLatents.scale(y);
        }
        
        
        //        if (_iter == 0 && n == 0 && j == 1)
        //        cout << n << " " << j << endl;
        //          phi.print();
        
        //        cout << phi.sum(0) << endl;
        //        phi.print();
        
        // Defines the subarrays of phi for latent variables, user observables, and item observables and updates the next shape parameter of theta and beta (gamma and kappa) by adding y_{ui} phi_{uik} to the nth row of gamma and the mth row of kappa (the first equation in steps 2 and 3 of the algorithm in the paper)
        _htheta.update_shape_next1(n, phiLatents);
        _hbeta.update_shape_next1(m, phiLatents);

                //        phi.print();
      }
    } // End of loop over users/movies
    
    debug("htheta = %s", _htheta.expected_v().s().c_str());
    debug("hbeta = %s", _hbeta.expected_v().s().c_str());
    
    //----------------------------------
    // Updates for user parameters
    //----------------------------------
    
    // If there are latent characteristics...
    
    if (_k>0) {
      Array betarowsum(_k);
      // Saves the sums over items of expected values for each factor ( the second part of \gamma^{rte}_{uk})
      _hbeta.sum_rows(betarowsum);
      
      // Sets the prior rate based on expectations with current parameters, i.e.,\frac{\kappa^{shp}}{\kappa^{rte}}
      _htheta.set_prior_rate(_thetarate.expected_v(),
                             _thetarate.expected_logv());
      debug("adding %s to theta rate", _thetarate.expected_v().s().c_str());
      debug("betarowsum %s", betarowsum.s().c_str());
      
      // Adds the previous sum (betarowsum) to \frac{\kappa^{shp}}{\kappa^{shp}} in the next rate
      _htheta.update_rate_next(betarowsum);
      // Swaps the current and the next values for the parameters
      _htheta.swap();
      
      // Computes expectations and log expectations based on the new parameters
      _htheta.compute_expectations();
    }
    printf("iteration %d\n", _iter);
    
    _iter++;
  }
  
  _iter = 0;
  
  bool stop = false;
  
  uint32_t x = _k+_uc+_ic;
  //
  // Constructs the array for the parameters of the multinomial distribution
  Array phi(x);
  
  while (!stop) {
    // Stop if the max number of iterations is reached
    if (_iter > _env.max_iterations) {
      exit(0);
    }
    
    // Loop over users
    for (uint32_t n = 0; n < _n; ++n) {
      //      cout << "Loop over users " << n << endl;
      
      // Gets the matrix of items for each user and stores it in movies
      const vector<uint32_t> *movies = _ratings.get_movies(n);
      
      // Loop over each user's items
      for (uint32_t j = 0; movies && j < movies->size(); ++j) {
        //        cout << "Loop over items " << j <<endl;
        // Gets the code of the movie
        uint32_t m = (*movies)[j];
        
        //        //////// JCC: Tests
        //        IDMap seq2m = _ratings.seq2movie();
        //        IDMap seq2u = _ratings.seq2user();
        //
        //        uint32_t username = seq2u.at(n);
        //        uint32_t moviename = seq2m.at(j);
        //        cout << "User: " << n << " " << username << endl;
        //        cout << "Movie: " << j << " " << m << " " << moviename << endl;
        //        ////////
        
        // Get the movie rating
        yval_t y = _ratings.r(n,m);
        
        //        //////// JCC: Tests
        //        uint32_t y2 = _ratings.r(n,m);
        //        cout << _ratings.r(n,m) << endl;
        //        cout << "Rating: " << y << endl;
        //        cout << "Rating: " << y2 << endl;
        //        ////////
        
        // Finds phi from the current parameters of hbeta, htheta, hsigma, and hrho (the equation in step 1 of the algorithm in the paper)
        get_phi(_htheta, n, _hbeta, m, _hsigma, _hrho, _thetarate, _betarate, _ic, _uc, phi);
        
        //        get_phi(_htheta, n, _hbeta, m, _hsigma, _hrho, _ic, _uc, phi);
        
        //        cout << phi.sum(0) << endl;
        //        if (_iter < 2) {
        //          cout << "Phi: " << endl;
        //          phi.print();
        //        }
        
        //        _phi.set_row( r, phi);
        
        // Makes phi sum up to y to get y_{ui} phi_{uik}
        if (y > 1) {
          phi.scale(y);
        }
        
        
        //        if (_iter == 0 && n == 0 && j == 1)
        //        cout << n << " " << j << endl;
        //          phi.print();
        
        //        cout << phi.sum(0) << endl;
        //        phi.print();
        
        // Defines the subarrays of phi for latent variables, user observables, and item observables and updates the next shape parameter of theta and beta (gamma and kappa) by adding y_{ui} phi_{uik} to the nth row of gamma and the mth row of kappa (the first equation in steps 2 and 3 of the algorithm in the paper)
        Array phik(_k);
        Array phil(_ic);
        Array phim(_uc);
        
        if (_k>0) {
          phik.copy_from(phi.subarray(0,_k-1));
          _htheta.update_shape_next1(n, phik);
          _hbeta.update_shape_next1(m, phik);
        }
        
        if (_ic > 0) {
          phil.copy_from(phi.subarray(_k,_k+_ic-1));
          _hsigma.update_shape_next1(n, phil);
        }
        
        if ( _uc > 0) {
          phim.copy_from(phi.subarray(_k+_ic,_k+_ic+_uc-1));
          _hrho.update_shape_next1(m, phim);
        }
        
        if (_env.bias) {
          _thetabias.update_shape_next3(n, 0, phi[_k]);
          _betabias.update_shape_next3(m, 0, phi[_k+1]);
        }
        //        phi.print();
      }
    } // End of loop over users/movies
    
    debug("htheta = %s", _htheta.expected_v().s().c_str());
    debug("hbeta = %s", _hbeta.expected_v().s().c_str());
    
    //----------------------------------
    // Updates for user parameters
    //----------------------------------
    
    // If there are latent characteristics...
    
    if (_k>0) {
      Array betarowsum(_k);
      // Saves the sums over items of expected values for each factor ( the second part of \gamma^{rte}_{uk})
      _hbeta.sum_rows(betarowsum);
      
      // Sets the prior rate based on expectations with current parameters, i.e.,\frac{\kappa^{shp}}{\kappa^{rte}}
      _htheta.set_prior_rate(_thetarate.expected_v(),
                             _thetarate.expected_logv());
      debug("adding %s to theta rate", _thetarate.expected_v().s().c_str());
      debug("betarowsum %s", betarowsum.s().c_str());
      
      // Adds the previous sum (betarowsum) to \frac{\kappa^{shp}}{\kappa^{shp}} in the next rate
      _htheta.update_rate_next(betarowsum);
      // Swaps the current and the next values for the parameters
      _htheta.swap();
      
      // Computes expectations and log expectations based on the new parameters
      _htheta.compute_expectations();
    }
    
    // If there are observed item characteristics...
    if (_ic > 0) {
      // Computes the sums of user characteristics (second term for \mu_{ul}^{rte}
      Array itemSum(_ic);
      _ratings._itemObs.weighted_colsum(_betarate.expected_inv(),itemSum);
      //      _ratings._itemObs.colsum(itemSum);
      
      // Sets the prior rate based on weighted expectations with current parameters, i.e.,x_l * e/(ca) * \frac{\kappa^{shp}}{\kappa^{rte}}
      _hsigma.set_prior_rate_scaled(_thetarate.expected_v(),_env.e/(_env.c*_env.a), _ratings._itemObsScale);
      
      // Adds the previous sum (itemSum) to \frac{\kappa^{shp}}{\kappa^{shp}} in the next rate
      _hsigma.update_rate_next(itemSum);
      // Swaps the current and the next values for the parameters
      _hsigma.swap();
      
      // Computes expectations and log expectations based on the new parameters
      _hsigma.compute_expectations();
    }
    
    //    cout << "Theta: " << endl;
    //    _htheta.shape_curr().print();
    //    _htheta.rate_curr().print();
    //
    //    cout << "Sigma: " << endl;
    //    _hsigma.shape_curr().print();
    //    _hsigma.rate_curr().print();
    
    //----------------------------------
    // Updates for item parameters
    //----------------------------------
    
    // If there are latent variables...
    if (_k>0) {
      Array thetarowsum(_k);
      // Saves the sums of expected values over users for each factor (the second part of \lambda^{rte}_{ik})
      _htheta.sum_rows(thetarowsum);
      
      // Sets the prior rate based on expectations with current parameters, i.e., \frac{\tau^{shp}}{\tau^{rte}}
      _hbeta.set_prior_rate(_betarate.expected_v(),
                            _betarate.expected_logv());
      
      // Adds the previous sum to \frac{\tau^{shp}}{\tau^{shp}} in the next rate
      _hbeta.update_rate_next(thetarowsum);
      // Swaps the current and the next values for the parameters
      _hbeta.swap();
      
      // Computes expectations and log expectations based on the new parameters
      _hbeta.compute_expectations();
    }
    
    // If there are observed user characteristics...
    if (_ic > 0) {
      // Computes the sums of user characteristics (second term for \mu_{ul}^{rte}
      Array userSum(_uc);
      
      //      _thetarate.shape_curr().print();
      //      _thetarate.rate_curr().print();
      //      _thetarate.expected_inv().print();
      _ratings._userObs.weighted_colsum(_thetarate.expected_inv(),userSum);
      //      _ratings._userObs.colsum(userSum);
      
      // Sets the prior rate based on expectations with current parameters, i.e.,\frac{\tau^{shp}}{\tau^{rte}}
      _hrho.set_prior_rate_scaled(_betarate.expected_v(),_env.f/(_env.c*_env.a),_ratings._userObsScale);
      
      // Adds the previous sum to \frac{\kappa^{shp}}{\kappa^{shp}} in the next rate
      _hrho.update_rate_next(userSum);
      // Swaps the current and the next values for the parameters
      _hrho.swap();
      
      // Computes expectations and log expectations based on the new parameters
      //      cout << "rho " << _iter << endl;
      _hrho.compute_expectations();
    }
    
    //    cout << "Beta: " << endl;
    //    _hbeta.shape_curr().print();
    //    _hbeta.rate_curr().print();
    //
    //    cout << "Rho: " << endl;
    //    _hrho.shape_curr().print();
    //    _hrho.rate_curr().print();
    
    //-------------------------------------
    // Update the rate parameters for the popularity and activity parameters
    //-------------------------------------
    
    // Adds Ka+Le to the shape parameter of xi
    _thetarate.update_shape_next(_k*_env.a+_ic*_env.e);
    
    // If there are latent variables...
    if (_k>0) {
      Array thetacolsum(_n);
      // Computes the second term of \kappa^{rte}_{uk})
      _htheta.sum_cols(thetacolsum);
      
      // Adds the new term to the rate parameter of xi
      _thetarate.update_rate_next(thetacolsum);
      
      //      thetacolsum.print();
    }
    
    // With unobserved item characteristics...
    if (_ic > 0) {
      Array sigmacolsum(_n);
      // Computes the third term of \kappa^{rte}_{uk})
      _hsigma.sum_cols_weight(_ratings._itemObsScale,sigmacolsum);
      //      _ratings._itemObsScale.print();
      // Adds the new term to the rate parameter of xi
      //      _thetarate.update_rate_next(sigmacolsum.scale(_env.e/(_env.c*_env.a)));
      double scale = _env.e/(_env.c*_env.a);
      _thetarate.update_rate_next(sigmacolsum,scale);
      //      sigmacolsum.print();
    }
    
    debug("thetacolsum = %s", thetacolsum.s().c_str());
    
    // Swaps the current and the next values for the parameters
    _thetarate.swap();
    // Computes the expectations with the (new) current values
    _thetarate.compute_expectations();
    
    
    // Adds Kc+Mf to the shape parameter of eta
    _betarate.update_shape_next(_k*_env.c+_uc*_env.f);
    
    // If there are latent variables...
    if (_k>0) {
      
      Array betacolsum(_m);
      
      // Computes the second term of \tau^{rte}_{uk})
      _hbeta.sum_cols(betacolsum);
      // Adds the new term to the rate parameter of eta
      _betarate.update_rate_next(betacolsum);
    }
    
    // With unobserved user characteristics...
    if (_uc > 0) {
      Array rhocolsum(_m);
      // Computes the third term of \tau^{rte}_{uk})
      _hrho.sum_cols_weight(_ratings._userObsScale,rhocolsum);
      // Adds the new term to the rate parameter of eta
      _betarate.update_rate_next(rhocolsum.scale(_env.f/(_env.c*_env.a)));
    }
    
    debug("betacolsum = %s", betacolsum.s().c_str());
    
    // Swaps the current and the next values for the parameters
    _betarate.swap();
    // Computes the expectations with the (new) current values
    _betarate.compute_expectations();
    
    //    if (_iter == 59 || _iter == 0) {
    //      cout << "Beta shape: " << endl;
    //      _hbeta.shape_curr().print();
    //      cout << "Beta rate: " << endl;
    //      _hbeta.rate_curr().print();
    //      cout << "Theta shape: " << endl;
    //      _htheta.shape_curr().print();
    //      cout << "Theta rate: " << endl;
    //      _htheta.rate_curr().print();
    //
    //      cout << "Eta shape: " << endl;
    //      _betarate.shape_curr().print();
    //      cout << "Eta rate: " << endl;
    //      _betarate.rate_curr().print();
    //      cout << "Xi shape: " << endl;
    //      _thetarate.shape_curr().print();
    //      cout << "Xi rate: " << endl;
    //      _thetarate.rate_curr().print();
    //    }
    
    printf("iteration %d\n", _iter);
    
    Array betaMean(_k);
    Array thetaMean(_k);
    Array sigmaMean(_ic);
    Array rhoMean(_uc);
    
    _hbeta.expected_means(betaMean);
    _htheta.expected_means(thetaMean);
    _hsigma.expected_means(sigmaMean);
    _hrho.expected_means(rhoMean);
    
    double xiMean = _betarate.expected_mean();
    double etaMean = _thetarate.expected_mean();
    
    betaMeans.set_col(_iter+1, betaMean);
    thetaMeans.set_col(_iter+1, thetaMean);
    sigmaMeans.set_col(_iter+1, sigmaMean);
    rhoMeans.set_col(_iter+1, rhoMean);
    
    xiMeans.set(_iter+1,xiMean);
    etaMeans.set(_iter+1,etaMean);
    
    //    _betarate.expected_v().print();
    //    _thetarate.expected_v().print();
    
    fflush(stdout);
    if (_iter % _env.reportfreq == 0) {
      
      compute_likelihood(false);
      stop = compute_likelihood(true);
      //compute_rmse();
      save_model();
      // Computes and saves number of relevant recommendations among best ranked items
      compute_precision(false);
      // Computes and saves average ranking of items in test set
      compute_itemrank(false);
      //gen_ranking_for_users(false);
      if (_env.logl)
        logl();
      
      string nameBeta = string("/betaMeans.tsv");
      string nameTheta = string("/thetaMeans.tsv");
      string nameSigma = string("/sigmaMeans.tsv");
      string nameRho = string("/rhoMeans.tsv");
      string nameXi = string("/xiMeans.tsv");
      string nameEta = string("/etaMeans.tsv");
      
      betaMeans.save(_env.outfname+"/"+Env::outfile_str(nameBeta));
      thetaMeans.save(_env.outfname+"/"+Env::outfile_str(nameTheta));
      sigmaMeans.save(_env.outfname+"/"+Env::outfile_str(nameSigma));
      rhoMeans.save(_env.outfname+"/"+Env::outfile_str(nameRho));
      
      xiMeans.save(_env.outfname+"/"+Env::outfile_str(nameXi));
      etaMeans.save(_env.outfname+"/"+Env::outfile_str(nameEta));
    }
    
    if (stop) {
      do_on_stop();
      //      _thetarate.expected_inv().print();
    } else {
      // Saves the matrices in files
      if (_env.save_state_now) {
        lerr("Saving state at iteration %d duration %d secs", _iter, duration());
        do_on_stop();
      }
      _iter++;
      
    }
    //    if ( _iter == 55) {
    //      _betarate.shape_curr().print();
    //      _betarate.rate_curr().print();
    //    }
  }
}

// Variation of the main method that first fits the latents, then the user observables, and finally the item observables. This is run repeatedly until convergence.
void
HGAPRec::vb_hier_cycles()
{
  // Initial values of the parameters, according to the posterior plus a random shock
  
  initialize();
  //  _betarate.shape_curr().print();
  //  _betarate.rate_curr().print();
  //  _thetarate.shape_curr().print();
  //  _thetarate.rate_curr().print();
  //  _hbeta.shape_curr().print();
  //  _hbeta.rate_curr().print();
  //  _htheta.shape_curr().print();
  //  _htheta.rate_curr().print();
  //  _hsigma.shape_curr().print();
  //  _hsigma.rate_curr().print();
  //  _hrho.shape_curr().print();
  //  _hrho.rate_curr().print();
  
  cout << "Initialized" << endl;
  //  cout << _env.reportfreq << endl;
  
  // Matrices that save the means of variables
  Matrix betaMeans(_k,_env.max_iterations+1);
  Matrix thetaMeans(_k,_env.max_iterations+1);
  Matrix sigmaMeans(_ic,_env.max_iterations+1);
  Matrix rhoMeans(_uc,_env.max_iterations+1);
  
  Array betaMean(_k);
  Array thetaMean(_k);
  Array sigmaMean(_ic);
  Array rhoMean(_uc);
  
  _hbeta.expected_means(betaMean);
  _htheta.expected_means(thetaMean);
  _hsigma.expected_means(sigmaMean);
  _hrho.expected_means(rhoMean);
  
  betaMeans.set_col(0, betaMean);
  thetaMeans.set_col(0, thetaMean);
  sigmaMeans.set_col(0, sigmaMean);
  rhoMeans.set_col(0, rhoMean);
  
  //lerr("htheta = %s", _htheta.rate_next().s().c_str());
  
  ///////////
  // Original code
  //  // What is bias?
  //  uint32_t x;
  //  if (_env.bias)
  //    x = _k+2;
  //  else
  //    x = _k;
  ///////////
  
  uint32_t x = _k+_uc+_ic;
  //
  // Constructs the array for the parameters of the multinomial distribution
  Array phi(x);
  
  bool stop = false;
  
  //  _ratings._itemObsScale.print();
  //  _ratings._userObsScale.print();
  
  while (!stop) {
    
    // Stop if the max number of iterations is reached
    if (_iter > _env.max_iterations) {
      exit(0);
    }
    
//    int mod50 = (_iter+1) % 100;
    int cycle = (_iter+1) % 2;
    
    // Loop over users
    for (uint32_t n = 0; n < _n; ++n) {
      //      cout << "Loop over users " << n << endl;
      
      // Gets the matrix of items for each user and stores it in movies
      const vector<uint32_t> *movies = _ratings.get_movies(n);
      
      // Loop over each user's items
      for (uint32_t j = 0; movies && j < movies->size(); ++j) {
        //        cout << "Loop over items " << j <<endl;
        // Gets the code of the movie
        uint32_t m = (*movies)[j];
        
        //        //////// JCC: Tests
        //        IDMap seq2m = _ratings.seq2movie();
        //        IDMap seq2u = _ratings.seq2user();
        //
        //        uint32_t username = seq2u.at(n);
        //        uint32_t moviename = seq2m.at(j);
        //        cout << "User: " << n << " " << username << endl;
        //        cout << "Movie: " << j << " " << m << " " << moviename << endl;
        //        ////////
        
        // Get the movie rating
        yval_t y = _ratings.r(n,m);
        
        //        //////// JCC: Tests
        //        uint32_t y2 = _ratings.r(n,m);
        //        cout << _ratings.r(n,m) << endl;
        //        cout << "Rating: " << y << endl;
        //        cout << "Rating: " << y2 << endl;
        //        ////////
        
        // Finds phi from the current parameters of hbeta, htheta, hsigma, and hrho (the equation in step 1 of the algorithm in the paper)
        get_phi(_htheta, n, _hbeta, m, _hsigma, _hrho, _thetarate, _betarate, _ic, _uc, phi);
        
        //        cout << phi.sum(0) << endl;
        //        if (_iter == 0) {
        //          cout << "Phi: " << endl;
        //          phi.print();
        //        }
        
        //        _phi.set_row( r, phi);
        
        // Makes phi sum up to y to get y_{ui} phi_{uik}
        if (y > 1) {
          phi.scale(y);
        }
        
        
        //        if (_iter == 0 && n == 0 && j == 1)
        //        cout << n << " " << j << endl;
        //          phi.print();
        
        //        cout << phi.sum(0) << endl;
        //        phi.print();
        
        // Defines the subarrays of phi for latent variables, user observables, and item observables and updates the next shape parameter of theta and beta (gamma and kappa) by adding y_{ui} phi_{uik} to the nth row of gamma and the mth row of kappa (the first equation in steps 2 and 3 of the algorithm in the paper)
        Array phik(_k);
        Array phil(_ic);
        Array phim(_uc);
        
        if (_k>0) {
          phik.copy_from(phi.subarray(0,_k-1));
          _htheta.update_shape_next1(n, phik);
          _hbeta.update_shape_next1(m, phik);
        }
        
        if (_ic > 0) {
          phil.copy_from(phi.subarray(_k,_k+_ic-1));
          _hsigma.update_shape_next1(n, phil);
        }
        
        if ( _uc > 0) {
          phim.copy_from(phi.subarray(_k+_ic,_k+_ic+_uc-1));
          _hrho.update_shape_next1(m, phim);
        }
        
        if (_env.bias) {
          _thetabias.update_shape_next3(n, 0, phi[_k]);
          _betabias.update_shape_next3(m, 0, phi[_k+1]);
        }
        //        phi.print();
      }
    } // End of loop over users/movies
    
    debug("htheta = %s", _htheta.expected_v().s().c_str());
    debug("hbeta = %s", _hbeta.expected_v().s().c_str());
    
    //----------------------------------
    // Updates for user parameters
    //----------------------------------
    
    // If there are latent characteristics...
    
    if (_k>0 && cycle == 0) {
      Array betarowsum(_k);
      // Saves the sums over items of expected values for each factor ( the second part of \gamma^{rte}_{uk})
      _hbeta.sum_rows(betarowsum);
      
      // Sets the prior rate based on expectations with current parameters, i.e.,\frac{\kappa^{shp}}{\kappa^{rte}}
      _htheta.set_prior_rate(_thetarate.expected_v(),
                             _thetarate.expected_logv());
      debug("adding %s to theta rate", _thetarate.expected_v().s().c_str());
      debug("betarowsum %s", betarowsum.s().c_str());
      
      // Adds the previous sum (betarowsum) to \frac{\kappa^{shp}}{\kappa^{shp}} in the next rate
      _htheta.update_rate_next(betarowsum);
      // Swaps the current and the next values for the parameters
      _htheta.swap();
      
      // Computes expectations and log expectations based on the new parameters
      _htheta.compute_expectations();
    }
    
    // If there are observed item characteristics...
    if (_ic > 0 && cycle == 1) {
      // Computes the sums of user characteristics (second term for \mu_{ul}^{rte}
      Array itemSum(_ic);
      _ratings._itemObs.colsum(itemSum);
      
      // Sets the prior rate based on expectations with current parameters, i.e.,\frac{\kappa^{shp}}{\kappa^{rte}}
      _hsigma.set_prior_rate_scaled(_thetarate.expected_v(),
                                    _thetarate.expected_logv(),_ratings._itemObsScale);
      
      // Adds the previous sum (itemSum) to \frac{\kappa^{shp}}{\kappa^{shp}} in the next rate
      _hsigma.update_rate_next(itemSum);
      // Swaps the current and the next values for the parameters
      _hsigma.swap();
      
      // Computes expectations and log expectations based on the new parameters
      _hsigma.compute_expectations();
    }
    
    //    cout << "Theta: " << endl;
    //    _htheta.shape_curr().print();
    //    _htheta.rate_curr().print();
    //
    //    cout << "Sigma: " << endl;
    //    _hsigma.shape_curr().print();
    //    _hsigma.rate_curr().print();
    
    //----------------------------------
    // Updates for item parameters
    //----------------------------------
    
    // If there are latent variables...
    if (_k>0 && cycle == 0 ) {
      Array thetarowsum(_k);
      // Saves the sums of expected values over users for each factor (the second part of \lambda^{rte}_{ik})
      _htheta.sum_rows(thetarowsum);
      
      // Sets the prior rate based on expectations with current parameters, i.e., \frac{\tau^{shp}}{\tau^{rte}}
      _hbeta.set_prior_rate(_betarate.expected_v(),
                            _betarate.expected_logv());
      
      // Adds the previous sum to \frac{\tau^{shp}}{\tau^{shp}} in the next rate
      _hbeta.update_rate_next(thetarowsum);
      // Swaps the current and the next values for the parameters
      _hbeta.swap();
      
      // Computes expectations and log expectations based on the new parameters
      _hbeta.compute_expectations();
    }
    
    // If there are observed user characteristics...
    if (_ic > 0 && cycle == 1) {
      // Computes the sums of user characteristics (second term for \mu_{ul}^{rte}
      Array userSum(_uc);
      _ratings._userObs.colsum(userSum);
      
      // Sets the prior rate based on expectations with current parameters, i.e.,\frac{\tau^{shp}}{\tau^{rte}}
      _hrho.set_prior_rate_scaled(_betarate.expected_v(),
                                  _betarate.expected_logv(),_ratings._userObsScale);
      
      // Adds the previous sum to \frac{\kappa^{shp}}{\kappa^{shp}} in the next rate
      _hrho.update_rate_next(userSum);
      // Swaps the current and the next values for the parameters
      _hrho.swap();
      
      // Computes expectations and log expectations based on the new parameters
      _hrho.compute_expectations();
    }
    
    //    cout << "Beta: " << endl;
    //    _hbeta.shape_curr().print();
    //    _hbeta.rate_curr().print();
    //
    //    cout << "Rho: " << endl;
    //    _hrho.shape_curr().print();
    //    _hrho.rate_curr().print();
    
    //-------------------------------------
    // Update the rate parameters for the popularity and activity parameters
    //-------------------------------------
    
    // If there are latent variables...
    if (_k>0 && cycle == 0) {
      Array thetacolsum(_n);
      // Computes the second term of \kappa^{rte}_{uk})
      _htheta.sum_cols(thetacolsum);
      
      // Adds (K+L)a to the shape parameter of xi
      _thetarate.update_shape_next((_k+_ic) * _thetarate.sprior());
      // Adds the new term to the rate parameter of xi
      _thetarate.update_rate_next(thetacolsum);
    }
    
    // With unobserved item characteristics...
    if (_ic > 0 && cycle == 1) {
      Array sigmacolsum(_n);
      // Computes the third term of \kappa^{rte}_{uk})
      _hsigma.sum_cols_weight(_ratings._itemObsScale,sigmacolsum);
      // Adds the new term to the rate parameter of xi
      _thetarate.update_rate_next(sigmacolsum);
    }
    debug("thetacolsum = %s", thetacolsum.s().c_str());
    
    // If there are latent variables...
    if (_k>0  && cycle == 0 ) {
      // Swaps the current and the next values for the parameters
      _thetarate.swap();
      // Computes the expectations with the (new) current values
      _thetarate.compute_expectations();
      
      //    _thetarate.rate_curr().print();
      
      Array betacolsum(_m);
      
      // Computes the second term of \tau^{rte}_{uk})
      _hbeta.sum_cols(betacolsum);
      // Adds (K+M)c to the shape parameter of eta
      _betarate.update_shape_next((_k+_uc) * _betarate.sprior());
      // Adds the new term to the rate parameter of eta
      _betarate.update_rate_next(betacolsum);
    }
    
    // With unobserved user characteristics...
    if (_uc > 0 && cycle == 1 ) {
      Array rhocolsum(_m);
      // Computes the third term of \tau^{rte}_{uk})
      _hrho.sum_cols_weight(_ratings._userObsScale,rhocolsum);
      // Adds the new term to the rate parameter of eta
      _betarate.update_rate_next(rhocolsum);
    }
    
    debug("betacolsum = %s", betacolsum.s().c_str());
    
    // Swaps the current and the next values for the parameters
    _betarate.swap();
    // Computes the expectations with the (new) current values
    _betarate.compute_expectations();
    
    printf("iteration %d\n", _iter);
    
    Array betaMean(_k);
    Array thetaMean(_k);
    Array sigmaMean(_ic);
    Array rhoMean(_uc);
    
    _hbeta.expected_means(betaMean);
    _htheta.expected_means(thetaMean);
    _hsigma.expected_means(sigmaMean);
    _hrho.expected_means(rhoMean);
    
    betaMeans.set_col(_iter+1, betaMean);
    thetaMeans.set_col(_iter+1, thetaMean);
    sigmaMeans.set_col(_iter+1, sigmaMean);
    rhoMeans.set_col(_iter+1, rhoMean);
    
    fflush(stdout);
    if (_iter % _env.reportfreq == 0) {
      compute_likelihood(false);
      stop = compute_likelihood(true);
      //compute_rmse();
      save_model();
      // Computes and saves number of relevant recommendations among best ranked items
      compute_precision(false);
      // Computes and saves average ranking of items in test set
      compute_itemrank(false);
      //gen_ranking_for_users(false);
      if (_env.logl)
        logl();
      
      string nameBeta = string("/betaMeans.tsv");
      string nameTheta = string("/thetaMeans.tsv");
      string nameSigma = string("/sigmaMeans.tsv");
      string nameRho = string("/rhoMeans.tsv");
      
      betaMeans.save(_env.outfname+"/"+Env::outfile_str(nameBeta));
      thetaMeans.save(_env.outfname+"/"+Env::outfile_str(nameTheta));
      sigmaMeans.save(_env.outfname+"/"+Env::outfile_str(nameSigma));
      rhoMeans.save(_env.outfname+"/"+Env::outfile_str(nameRho));
    }
    
    // Saves matrices in files
    if (_env.save_state_now) {
      lerr("Saving state at iteration %d duration %d secs", _iter, duration());
      do_on_stop();
    }
    _iter++;
    
//    if (stop) {
//      do_on_stop();
//    } else {
//      // Saves the matrices in files
//      if (_env.save_state_now) {
//        lerr("Saving state at iteration %d duration %d secs", _iter, duration());
//        do_on_stop();
//      }
//      _iter++;
//      
//    }
    
    //    if ( _iter == 55) {
    //      _betarate.shape_curr().print();
    //      _betarate.rate_curr().print();
    //    }
  }
}

// Variation of the main method that first fits the latents, then the user observables, and finally the item observables. This is run repeatedly until convergence.
void
HGAPRec::vb_hier_cycles2()
{
  // Initial values of the parameters, according to the posterior plus a random shock
  
  initialize();
  //  _betarate.shape_curr().print();
  //  _betarate.rate_curr().print();
  //  _thetarate.shape_curr().print();
  //  _thetarate.rate_curr().print();
  //  _hbeta.shape_curr().print();
  //  _hbeta.rate_curr().print();
  //  _htheta.shape_curr().print();
  //  _htheta.rate_curr().print();
  //  _hsigma.shape_curr().print();
  //  _hsigma.rate_curr().print();
  //  _hrho.shape_curr().print();
  //  _hrho.rate_curr().print();
  
  cout << "Initialized" << endl;
  //  cout << _env.reportfreq << endl;
  
  // Matrices that save the means of variables
  Matrix betaMeans(_k,_env.max_iterations+1);
  Matrix thetaMeans(_k,_env.max_iterations+1);
  Matrix sigmaMeans(_ic,_env.max_iterations+1);
  Matrix rhoMeans(_uc,_env.max_iterations+1);
  
  Array betaMean(_k);
  Array thetaMean(_k);
  Array sigmaMean(_ic);
  Array rhoMean(_uc);
  
  _hbeta.expected_means(betaMean);
  _htheta.expected_means(thetaMean);
  _hsigma.expected_means(sigmaMean);
  _hrho.expected_means(rhoMean);
  
  betaMeans.set_col(0, betaMean);
  thetaMeans.set_col(0, thetaMean);
  sigmaMeans.set_col(0, sigmaMean);
  rhoMeans.set_col(0, rhoMean);
  
  //lerr("htheta = %s", _htheta.rate_next().s().c_str());
  
  ///////////
  // Original code
  //  // What is bias?
  //  uint32_t x;
  //  if (_env.bias)
  //    x = _k+2;
  //  else
  //    x = _k;
  ///////////
  
  uint32_t x = _k+_uc+_ic;
  //
  // Constructs the array for the parameters of the multinomial distribution
  Array phi(x);
  
  bool stop = false;
  
  //  _ratings._itemObsScale.print();
  //  _ratings._userObsScale.print();
  
  while (!stop) {
    
    // Stop if the max number of iterations is reached
    if (_iter > _env.max_iterations) {
      exit(0);
    }
    
    //    int mod50 = (_iter+1) % 100;
    int cycle = (_iter+1) % 2;
    cycle = 1-cycle;
    
    // Loop over users
    for (uint32_t n = 0; n < _n; ++n) {
      //      cout << "Loop over users " << n << endl;
      
      // Gets the matrix of items for each user and stores it in movies
      const vector<uint32_t> *movies = _ratings.get_movies(n);
      
      // Loop over each user's items
      for (uint32_t j = 0; movies && j < movies->size(); ++j) {
        //        cout << "Loop over items " << j <<endl;
        // Gets the code of the movie
        uint32_t m = (*movies)[j];
        
        //        //////// JCC: Tests
        //        IDMap seq2m = _ratings.seq2movie();
        //        IDMap seq2u = _ratings.seq2user();
        //
        //        uint32_t username = seq2u.at(n);
        //        uint32_t moviename = seq2m.at(j);
        //        cout << "User: " << n << " " << username << endl;
        //        cout << "Movie: " << j << " " << m << " " << moviename << endl;
        //        ////////
        
        // Get the movie rating
        yval_t y = _ratings.r(n,m);
        
        //        //////// JCC: Tests
        //        uint32_t y2 = _ratings.r(n,m);
        //        cout << _ratings.r(n,m) << endl;
        //        cout << "Rating: " << y << endl;
        //        cout << "Rating: " << y2 << endl;
        //        ////////
        
        // Finds phi from the current parameters of hbeta, htheta, hsigma, and hrho (the equation in step 1 of the algorithm in the paper)
        get_phi(_htheta, n, _hbeta, m, _hsigma, _hrho, _ic, _uc, phi);
        
        //        cout << phi.sum(0) << endl;
        //        if (_iter == 0) {
        //          cout << "Phi: " << endl;
        //          phi.print();
        //        }
        
        //        _phi.set_row( r, phi);
        
        // Makes phi sum up to y to get y_{ui} phi_{uik}
        if (y > 1) {
          phi.scale(y);
        }
        
        
        //        if (_iter == 0 && n == 0 && j == 1)
        //        cout << n << " " << j << endl;
        //          phi.print();
        
        //        cout << phi.sum(0) << endl;
        //        phi.print();
        
        // Defines the subarrays of phi for latent variables, user observables, and item observables and updates the next shape parameter of theta and beta (gamma and kappa) by adding y_{ui} phi_{uik} to the nth row of gamma and the mth row of kappa (the first equation in steps 2 and 3 of the algorithm in the paper)
        Array phik(_k);
        Array phil(_ic);
        Array phim(_uc);
        
        if (_k>0) {
          phik.copy_from(phi.subarray(0,_k-1));
          _htheta.update_shape_next1(n, phik);
          _hbeta.update_shape_next1(m, phik);
        }
        
        if (_ic > 0) {
          phil.copy_from(phi.subarray(_k,_k+_ic-1));
          _hsigma.update_shape_next1(n, phil);
        }
        
        if ( _uc > 0) {
          phim.copy_from(phi.subarray(_k+_ic,_k+_ic+_uc-1));
          _hrho.update_shape_next1(m, phim);
        }
        
        if (_env.bias) {
          _thetabias.update_shape_next3(n, 0, phi[_k]);
          _betabias.update_shape_next3(m, 0, phi[_k+1]);
        }
        //        phi.print();
      }
    } // End of loop over users/movies
    
    debug("htheta = %s", _htheta.expected_v().s().c_str());
    debug("hbeta = %s", _hbeta.expected_v().s().c_str());
    
    //----------------------------------
    // Updates for user parameters
    //----------------------------------
    
    // If there are latent characteristics...
    
    if (_k>0 && cycle == 0) {
      Array betarowsum(_k);
      // Saves the sums over items of expected values for each factor ( the second part of \gamma^{rte}_{uk})
      _hbeta.sum_rows(betarowsum);
      
      // Sets the prior rate based on expectations with current parameters, i.e.,\frac{\kappa^{shp}}{\kappa^{rte}}
      _htheta.set_prior_rate(_thetarate.expected_v(),
                             _thetarate.expected_logv());
      debug("adding %s to theta rate", _thetarate.expected_v().s().c_str());
      debug("betarowsum %s", betarowsum.s().c_str());
      
      // Adds the previous sum (betarowsum) to \frac{\kappa^{shp}}{\kappa^{shp}} in the next rate
      _htheta.update_rate_next(betarowsum);
      // Swaps the current and the next values for the parameters
      _htheta.swap();
      
      // Computes expectations and log expectations based on the new parameters
      _htheta.compute_expectations();
    }
    
    // If there are observed item characteristics...
    if (_ic > 0 && cycle == 1) {
      // Computes the sums of user characteristics (second term for \mu_{ul}^{rte}
      Array itemSum(_ic);
      _ratings._itemObs.colsum(itemSum);
      
      // Sets the prior rate based on expectations with current parameters, i.e.,\frac{\kappa^{shp}}{\kappa^{rte}}
      _hsigma.set_prior_rate_scaled(_thetarate.expected_v(),
                                    _thetarate.expected_logv(),_ratings._itemObsScale);
      
      // Adds the previous sum (itemSum) to \frac{\kappa^{shp}}{\kappa^{shp}} in the next rate
      _hsigma.update_rate_next(itemSum);
      // Swaps the current and the next values for the parameters
      _hsigma.swap();
      
      // Computes expectations and log expectations based on the new parameters
      _hsigma.compute_expectations();
    }
    
    //    cout << "Theta: " << endl;
    //    _htheta.shape_curr().print();
    //    _htheta.rate_curr().print();
    //
    //    cout << "Sigma: " << endl;
    //    _hsigma.shape_curr().print();
    //    _hsigma.rate_curr().print();
    
    //----------------------------------
    // Updates for item parameters
    //----------------------------------
    
    // If there are latent variables...
    if (_k>0 && cycle == 0 ) {
      Array thetarowsum(_k);
      // Saves the sums of expected values over users for each factor (the second part of \lambda^{rte}_{ik})
      _htheta.sum_rows(thetarowsum);
      
      // Sets the prior rate based on expectations with current parameters, i.e., \frac{\tau^{shp}}{\tau^{rte}}
      _hbeta.set_prior_rate(_betarate.expected_v(),
                            _betarate.expected_logv());
      
      // Adds the previous sum to \frac{\tau^{shp}}{\tau^{shp}} in the next rate
      _hbeta.update_rate_next(thetarowsum);
      // Swaps the current and the next values for the parameters
      _hbeta.swap();
      
      // Computes expectations and log expectations based on the new parameters
      _hbeta.compute_expectations();
    }
    
    // If there are observed user characteristics...
    if (_ic > 0 && cycle == 1) {
      // Computes the sums of user characteristics (second term for \mu_{ul}^{rte}
      Array userSum(_uc);
      _ratings._userObs.colsum(userSum);
      
      // Sets the prior rate based on expectations with current parameters, i.e.,\frac{\tau^{shp}}{\tau^{rte}}
      _hrho.set_prior_rate_scaled(_betarate.expected_v(),
                                  _betarate.expected_logv(),_ratings._userObsScale);
      
      // Adds the previous sum to \frac{\kappa^{shp}}{\kappa^{shp}} in the next rate
      _hrho.update_rate_next(userSum);
      // Swaps the current and the next values for the parameters
      _hrho.swap();
      
      // Computes expectations and log expectations based on the new parameters
      _hrho.compute_expectations();
    }
    
    //    cout << "Beta: " << endl;
    //    _hbeta.shape_curr().print();
    //    _hbeta.rate_curr().print();
    //
    //    cout << "Rho: " << endl;
    //    _hrho.shape_curr().print();
    //    _hrho.rate_curr().print();
    
    //-------------------------------------
    // Update the rate parameters for the popularity and activity parameters
    //-------------------------------------
    
    // If there are latent variables...
    if (_k>0 && cycle == 0) {
      Array thetacolsum(_n);
      // Computes the second term of \kappa^{rte}_{uk})
      _htheta.sum_cols(thetacolsum);
      
      // Adds (K+L)a to the shape parameter of xi
      _thetarate.update_shape_next((_k+_ic) * _thetarate.sprior());
      // Adds the new term to the rate parameter of xi
      _thetarate.update_rate_next(thetacolsum);
    }
    
    // With unobserved item characteristics...
    if (_ic > 0 && cycle == 1) {
      Array sigmacolsum(_n);
      // Computes the third term of \kappa^{rte}_{uk})
      _hsigma.sum_cols_weight(_ratings._itemObsScale,sigmacolsum);
      // Adds the new term to the rate parameter of xi
      _thetarate.update_rate_next(sigmacolsum);
    }
    debug("thetacolsum = %s", thetacolsum.s().c_str());
    
    // If there are latent variables...
    if (_k>0  && cycle == 0 ) {
      // Swaps the current and the next values for the parameters
      _thetarate.swap();
      // Computes the expectations with the (new) current values
      _thetarate.compute_expectations();
      
      //    _thetarate.rate_curr().print();
      
      Array betacolsum(_m);
      
      // Computes the second term of \tau^{rte}_{uk})
      _hbeta.sum_cols(betacolsum);
      // Adds (K+M)c to the shape parameter of eta
      _betarate.update_shape_next((_k+_uc) * _betarate.sprior());
      // Adds the new term to the rate parameter of eta
      _betarate.update_rate_next(betacolsum);
    }
    
    // With unobserved user characteristics...
    if (_uc > 0 && cycle == 1 ) {
      Array rhocolsum(_m);
      // Computes the third term of \tau^{rte}_{uk})
      _hrho.sum_cols_weight(_ratings._userObsScale,rhocolsum);
      // Adds the new term to the rate parameter of eta
      _betarate.update_rate_next(rhocolsum);
    }
    
    debug("betacolsum = %s", betacolsum.s().c_str());
    
    // Swaps the current and the next values for the parameters
    _betarate.swap();
    // Computes the expectations with the (new) current values
    _betarate.compute_expectations();
    
    printf("iteration %d\n", _iter);
    
    Array betaMean(_k);
    Array thetaMean(_k);
    Array sigmaMean(_ic);
    Array rhoMean(_uc);
    
    _hbeta.expected_means(betaMean);
    _htheta.expected_means(thetaMean);
    _hsigma.expected_means(sigmaMean);
    _hrho.expected_means(rhoMean);
    
    betaMeans.set_col(_iter+1, betaMean);
    thetaMeans.set_col(_iter+1, thetaMean);
    sigmaMeans.set_col(_iter+1, sigmaMean);
    rhoMeans.set_col(_iter+1, rhoMean);
    
    fflush(stdout);
    if (_iter % _env.reportfreq == 0) {
      compute_likelihood(false);
      stop = compute_likelihood(true);
      //compute_rmse();
      save_model();
      // Computes and saves number of relevant recommendations among best ranked items
      compute_precision(false);
      // Computes and saves average ranking of items in test set
      compute_itemrank(false);
      //gen_ranking_for_users(false);
      if (_env.logl)
        logl();
      
      string nameBeta = string("/betaMeans.tsv");
      string nameTheta = string("/thetaMeans.tsv");
      string nameSigma = string("/sigmaMeans.tsv");
      string nameRho = string("/rhoMeans.tsv");
      
      betaMeans.save(_env.outfname+"/"+Env::outfile_str(nameBeta));
      thetaMeans.save(_env.outfname+"/"+Env::outfile_str(nameTheta));
      sigmaMeans.save(_env.outfname+"/"+Env::outfile_str(nameSigma));
      rhoMeans.save(_env.outfname+"/"+Env::outfile_str(nameRho));
    }
    
    // Saves matrices in files
    if (_env.save_state_now) {
      lerr("Saving state at iteration %d duration %d secs", _iter, duration());
      do_on_stop();
    }
    _iter++;
    
    //    if (stop) {
    //      do_on_stop();
    //    } else {
    //      // Saves the matrices in files
    //      if (_env.save_state_now) {
    //        lerr("Saving state at iteration %d duration %d secs", _iter, duration());
    //        do_on_stop();
    //      }
    //      _iter++;
    //
    //    }
    
    //    if ( _iter == 55) {
    //      _betarate.shape_curr().print();
    //      _betarate.rate_curr().print();
    //    }
  }
}

// Calculates log likelihood. Validation tells whether it should be calculated for the validation or test set. If validation, also check the stopping criterion.
bool
HGAPRec::compute_likelihood(bool validationLikelihood)
{
//  testLogLikelihood = new Matrix(7998,5);
//  cout << validation << endl;
  // k: index of ratings in validation/training
  uint32_t k = 0, kzeros = 0, kones = 0;
  // s: stores the sum of log likelihoods
  double s = .0, szeros = 0, sones = 0;
  
  // Saves either _validation_map or _test_map in mp (depending on whether the parameter validation is true)
  CountMap *mp = NULL;
  FILE *ff = NULL;
  if (validationLikelihood) {
    mp = &_ratings._validation_map;
    ff = _vf;
  } else {
    mp = &_ratings._test_map;
    ff = _tf;
  }

//  cout << mp->size() << endl;
  
  // Loop over ratings
  for (CountMap::const_iterator i = mp->begin();
       i != mp->end(); ++i) {
    // Extracts the indices for the user and item
    const Rating &e = i->first;
    uint32_t n = e.first;
    uint32_t m = e.second;

    // Extracts the rating
    yval_t r = i->second;
    
    // Finds the log likelihood and adds it
//    cout << k << endl;
    double u = _env.hier ? rating_likelihood_hier(n,m,r) : rating_likelihood(n,m,r);
  
//    double rate;
//    double likelihood;
//    rating_likelihood_hier_return(n,m,r,rate,likelihood);
    
//    cout << "r: " << r << " rate: " << rate << " likelihood: " << likelihood << endl;
//    IDMap userMap = _ratings.seq2user();
//    IDMap itemMap = _ratings.seq2movie();
//    testLogLikelihood->set(k,0,r);
//    testLogLikelihood->set(k,1,rate);
//    testLogLikelihood->set(k,2,userMap[n]);
//    testLogLikelihood->set(k,3,itemMap[m]);
//    testLogLikelihood->set(k,4,likelihood);

    s += u;
    k += 1;
    
//        cout << k << endl;
  }
  
//  if ( validationLikelihood) {
//    string name = _env.outfname+"/"+Env::outfile_str("/likelihoodsValidation.tsv");
////    cout << name << endl;
////    testLogLikelihood->save(name.c_str());
//  } else {
//    string name = _env.outfname+"/"+Env::outfile_str("/likelihoodsTest.tsv");
//    cout << name << endl;
////    testLogLikelihood->save(name.c_str());
//  }
  
  double a = .0;
  info("s = %.5f\n", s);
  fprintf(ff, "%d\t%d\t%.9f\t%d\n", _iter, duration(), s / k, k);
  fflush(ff);
  
  // Average log likelihood
  a = s / k;
//cout << k << endl;
  if (validationLikelihood) {
    cout << "Validation: " ;
    cout << setprecision(10) << s << "\t" << k << "\t" << a << endl;
  } else {
    cout << "Log-likelihood: (total,number,average)" << endl;
    cout << "Test: " ;
    cout << setprecision(10) << s << "\t" << k << "\t" << a << endl;
  }
  
  bool stop = false;
  
  // Check validation criterion
  if (validationLikelihood) {
    int why = -1;
    
    // Check stopping criteria every iteration after 30
    if (_iter > 30) {
      
      cout << " Likelihood change: " << fabs((a - _prev_h) / _prev_h) << endl;
      cout << a << " " << a -  _prev_h << endl;
      // If log likelihood increased, is not zero, and it increased less than 0.000001 of the previous value, set why to zero
      if (a >= _prev_h && _prev_h != 0 && fabs((a - _prev_h) / _prev_h) < 0.000001) {
        
        stop = true;
        why = 0;
      }
      // Count the number of times in a row that the likelihood decreased
      else if (a < _prev_h)
        _nh++;
      else if (a > _prev_h)
        _nh = 0;
      
      if (_nh > 5) { // be robust to small fluctuations in predictive likelihood
        why = 1;
        stop = true;
      }
    }
    // Store average log likelihood in _prev_h (previous likelihood)
    _prev_h = a;
    string name = _env.outfname+"/"+_env.prefix +"/max.txt";
    FILE *f = fopen(name.c_str(), "w");
    fprintf(f, "%d\t%d\t%.5f\t%d\n",
            _iter, duration(), a, why);
    fclose(f);
  }
//  delete testLogLikelihood;
  return stop;
}

double
HGAPRec::rating_likelihood(uint32_t p, uint32_t q, yval_t y) const
{
  assert (!_env.bias || (_env.bias && !_env.mle_user && !_env.mle_item));

  const double **etheta = _theta.expected_v().const_data();
  const double **ebeta = _beta.expected_v().const_data();

  const double **td = _theta_mle.data();
  const double **bd = _beta_mle.data();
  
  double s = .0;
  for (uint32_t k = 0; k < _k; ++k) {
    if (_env.mle_user) 
      s += td[p][k] * ebeta[q][k];
    else if (_env.mle_item || _env.canny)
      s += etheta[p][k] * bd[q][k];
    else
      s += etheta[p][k] * ebeta[q][k];
  }
  
  if (_env.bias) {
    const double **ethetabias = _thetabias.expected_v().const_data();
    const double **ebetabias = _betabias.expected_v().const_data();
    s += ethetabias[p][0] + ebetabias[q][0];
  } 
  
  if (s < 1e-30)
    s = 1e-30;
  
  if (_env.binary_data)
    return y == 0 ? -s : log(1 - exp(-s));    
  return y * log(s) - s - log_factorial(y);
}

// Computes the log likelihood for the rating of user p for item q
double
HGAPRec::rating_likelihood_hier(uint32_t u, uint32_t i, yval_t y) const
{
  // Gets the current expected values of theta, beta, sigma, rho, and the expected value of the inverse of xi and eta
  const double **etheta = _htheta.expected_v().const_data();
  const double **ebeta = _hbeta.expected_v().const_data();
  const double **esigma = _hsigma.expected_v().const_data();
  const double **erho = _hrho.expected_v().const_data();
  const double *einvxi = _thetarate.expected_inv().const_data();
  const double *einveta = _betarate.expected_inv().const_data();
  double **itemChar = _ratings._itemObs.data();
  double **userChar = _ratings._userObs.data();
  
  // Finds the sum of dot products \theta_u\beta_i+\sigma_u\x_i+\w_u\rho_i

  
  double s = .0;
  for (uint32_t k = 0; k < _k; ++k) {
    s += etheta[u][k] * ebeta[i][k];
  }
  for (uint32_t l = 0; l < _ic; ++l) {
    s += esigma[u][l] * itemChar[i][l] * einveta[i];
  }
  for (uint32_t m = 0; m < _uc; ++m) {
    s += userChar[u][m] * einvxi[u] * erho[i][m];
  }
  
//  if ( _iter == 130 && i == 263  ) {
//    
//    cout << _betarate.rate_curr().get(i) << " " << _betarate.shape_curr().get(i) << endl;
//    cout << einveta[i] << " " << einvxi[u] << endl;
//    
//    for (uint32_t k = 0; k < _k; ++k) {
//      cout << etheta[u][k] * ebeta[i][k] << " ";
//    }
//        cout << endl;
//    for (uint32_t l = 0; l < _ic; ++l) {
//      cout << esigma[u][l] * itemChar[i][l] * einveta[i] << " ";
//      cout << esigma[u][l] << " " << itemChar[i][l]<< " " << einveta[i] << " " << _betarate.rate_curr().get(i)<< " " << _betarate.shape_curr().get(i) <<": ";
//    }
//        cout << endl;
//    for (uint32_t m = 0; m < _uc; ++m) {
//      cout << userChar[u][m] * einvxi[u] * erho[i][m] << " ";
//      cout << erho[i][m] << " " << userChar[u][m] << " " << einvxi[u] << ": ";
//    }
//    
//    cout << endl;
//    
//    cout << _ratings.seq2user().at(u) << " " << _ratings.seq2movie().at(i) << endl;
//    cout << "y: " << y << " s: " << s << " ll: " << (y * log(s) - s - log_factorial(y)) << endl;
////    cout << "-------------" << i << endl;
//  }
  
  if (_env.bias) {
    const double **ethetabias = _thetabias.expected_v().const_data();
    const double **ebetabias = _betabias.expected_v().const_data();
    s += ethetabias[u][0] + ebetabias[i][0];
  } 
  
  if (s < 1e-30) {
    s = 1e-30;
  }
  
  // Returns log of Poisson pmf
  if (_env.binary_data)
    return y == 0 ? -s : log(1 - exp(-s));
  return y * log(s) - s - log_factorial(y);
}

// Computes the log likelihood for the rating of user p for item q, stores the rate and likelihood in parameters
//double
//HGAPRec::rating_likelihood_hier_return(uint32_t p, uint32_t q, yval_t y, double & rate, double & likelihood) const
//{
//  // Gets the current expected values of theta, beta, sigma, and rho
//  const double **etheta = _htheta.expected_v().const_data();
//  const double **ebeta = _hbeta.expected_v().const_data();
//  const double **esigma = _hsigma.expected_v().const_data();
//  const double **erho = _hrho.expected_v().const_data();
//  double **itemChar = _ratings._itemObs.data();
//  double **userChar = _ratings._userObs.data();
//  
//  // Finds the sum of dot products \theta_u\beta_i+\sigma_u\x_i+\w_u\rho_i
//  double s = .0;
//  for (uint32_t k = 0; k < _k; ++k)
//    s += etheta[p][k] * ebeta[q][k];
//  for (uint32_t l = 0; l < _ic; ++l)
//    s += esigma[p][l] * itemChar[q][l];
//  for (uint32_t m = 0; m < _uc; ++m)
//    s += userChar[p][m] * erho[q][m];
//  
//  
//  if (_env.bias) {
//    const double **ethetabias = _thetabias.expected_v().const_data();
//    const double **ebetabias = _betabias.expected_v().const_data();
//    s += ethetabias[p][0] + ebetabias[q][0];
//  }
//  
//  if (s < 1e-30)
//    s = 1e-30;
//  
//  int y2 = y;
////  cout << "y: " << y2 << " s: " << s << " ll: " << (y * log(s) - s - log_factorial(y)) << endl;
//  
//  rate = s;
//  likelihood = y * log(s) - s - log_factorial(y);
//  
//  // Returns log of Poisson pmf
//  if (_env.binary_data)
//    return y == 0 ? -s : log(1 - exp(-s));
//  return y * log(s) - s - log_factorial(y);
//}

double
HGAPRec::log_factorial(uint32_t n)  const
{
  double v = log(1);
  for (uint32_t i = 2; i <= n; ++i)
    v += log(i);
  return v;
} 

void
HGAPRec::do_on_stop()
{
  save_model();
  gen_ranking_for_users(false);
}

double
HGAPRec::compute_rmse()
{
  string name = _env.outfname+"/"+_env.prefix +"/test_scores.tsv";
  FILE *outf = fopen(name.c_str(), "w");
  if (!_rf)  {
    printf("cannot open test_scores file:%s\n",  strerror(errno));
    exit(-1);
  }  
  //load_beta_and_theta();
  double s = .0;
  for (CountMap::const_iterator i = _ratings._test_map.begin();
       i != _ratings._test_map.end(); ++i) {
    const Rating &r = i->first;
    uint32_t n = r.first;
    uint32_t m = r.second;
    yval_t v = i->second;
    double u = _env.graphchi? prediction_score_chi(n,m) : \
      (_env.hier ? prediction_score_hier(n, m) : prediction_score(n, m));
    s += (u - v) * (u - v);
    fprintf(outf, "%d\t%.5f\n", v, u);
  }
  fprintf(_rf, "%.5f\n", sqrt(s / _ratings._test_map.size()));
  fflush(_rf);
  fclose(outf);
  return sqrt(s / _ratings._test_map.size());
}

//
double
HGAPRec::compute_itemrank(bool final)
{
  // Only runs every 100 iterations, otherwise leaves the method
  if (_iter % 100 == 0 && _iter > 0)
    final = true;

  if (!final)
    return .0;

  uint32_t total_users = 0;
  FILE *f = 0;
  string name = _env.outfname+"/"+_env.prefix +"/itemrank.tsv";
  f = fopen(name.c_str(), "w");
  string name2 = _env.outfname+"/"+_env.prefix +"/meanrank.txt";
  FILE *itemf = fopen(name2.c_str(), "w");
  if (!itemf)  {
    printf("cannot open meanrank file:%s\n",  strerror(errno));
    exit(-1);
  }  

  // Array of predicted ratings
  KVArray mlist(_m);
  // Array of observed ratings
  KVIArray ndcglist(_m);
  
  double sum_rank = .0;
  double sum_reciprocal_rank = .0;
  
  // Iterates over users in sample. Uses the sample of users generated in this::compute_precision()
  for (UserMap::const_iterator itr = _sampled_users.begin();
       itr != _sampled_users.end(); ++itr) {
    
    // Saves number of user in n
    uint32_t n = itr->first;
    
    // Loop over items
    for (uint32_t m = 0; m < _m; ++m) {
      Rating r(n,m);
      
      // Saves zero predicted rating if the observed rating is nonzero or if it is in the validation set, then moves on to the next movie
      if (_ratings.r(n,m) > 0 || is_validation(r)) { // skip training and validation
        mlist[m].first = m;
        mlist[m].second = .0;
        continue;
      }
      
      // Saves the predicted score in u
      double u = .0;
      if (_env.nmf)
        u = prediction_score_nmf(n, m);
      else if (_env.lda || _env.vwlda)
        u = prediction_score_lda(n, m);
      else if (_env.graphchi) {
        u = prediction_score_chi(n, m);
      } else
        u = _env.hier ? prediction_score_hier(n, m) : prediction_score(n, m);
      
      // Saves the predicted rating in mlist. If the rating is in the test set, saves the observed rating in ndcglist
      mlist[m].first = m;
      mlist[m].second = u;
    } // End of loop over items
    
    // Sorts predicted ratings
    mlist.sort_by_value();
    
    double rank_ui = .0;
    double reciprocal_rank_ui = .0;
    uint32_t ntestitems = 0;
    uint32_t nranked = 0;
    
    // Loops again over items
    for (uint32_t j = 0; j < mlist.size(); ++j) {
      // Saves the user's predicted rating for the current item as (m,pred)
      KV &kv = mlist[j];
      uint32_t m = kv.first;
      double pred = kv.second;
      Rating r(n, m);
      
      if (_ratings.r(n,m) == 0) // not in validation or training
        nranked++;
      
      // If the movie is in the test set
      CountMap::const_iterator itr = _ratings._test_map.find(r);
      if (itr != _ratings._test_map.end()) {
        // Compute whether it counts as a hit
        int v_ = itr->second;
        int v;
        if (_ratings.test_hit(v_))
          v = 1;
        else
          v = 0;
        
        // If it is a hit, add ranking to rank_ui, add reciprocal to reciprocal_rank_ui
        if (v > 0) {
          ntestitems++;
          const vector<uint32_t> *t = _ratings.get_users(m);
          fprintf(f, "%d\t%d\t%.5f\t%d\t%d\n", n, m, pred, j, t?t->size():0);
          rank_ui += (j+1);
          reciprocal_rank_ui += 1 / (j+1);
        }
      }
    } // End of loop over items
    
    // Computes average rank and average reciprocal of rank for hits for each user
    if (ntestitems > 0 && nranked > 0) {
      sum_rank += (rank_ui / nranked) / ntestitems;
      sum_reciprocal_rank += reciprocal_rank_ui / ntestitems;
      total_users++;
    }
  }
  fclose(f);
  
  // Saves average rank and average reciprocal rank
  fprintf(itemf, "%d\t%.5f\t%.5f\n",
          total_users,
          (double)sum_rank/total_users,
          (double)sum_reciprocal_rank/total_users);
  fclose(itemf);
  return(sum_rank / total_users);
}

// Computes the number of relevant recommendation among best ranked 10 and 100
void
HGAPRec::compute_precision(bool save_ranking_file)
{
  // Sets boolean to save ranking file as true every 100 iterations
  if (_iter % 100 == 0 && _iter > 0)
    save_ranking_file = true;
  
  // Initializes the variables that will store the number of relevant recommendations
  double mhits10 = 0, mhits100 = 0;
  double cumndcg10 = 0, cumndcg100 = 0; // Unused
  uint32_t total_users = 0;
  FILE *f = 0;
  if (save_ranking_file) {
    string name = _env.outfname+"/"+_env.prefix+"/ranking.tsv";
    f = fopen(name.c_str(), "w");
  }
  
  // Picks randomly a sample with either 1000 or one half of the users. Stores their indices in _sampled_users (a map of booleans)
  if (!save_ranking_file) {
    _sampled_users.clear();
    do {
      uint32_t n = gsl_rng_uniform_int(_r, _n);
      _sampled_users[n] = true;
    } while (_sampled_users.size() < 1000 && _sampled_users.size() < _n / 2);
  }
  
  // Matrix of predicted ratings
  KVArray mlist(_m);
  // Matrix of observed ratings
  KVIArray ndcglist(_m);
  
  // Iterates over users in sample
  for (UserMap::const_iterator itr = _sampled_users.begin();
       itr != _sampled_users.end(); ++itr) {
    
    // Saves number of user in n
    uint32_t n = itr->first;
    
    // Loop over items
    for (uint32_t m = 0; m < _m; ++m) {
      Rating r(n,m);
      
      // Saves zero predicted and observed rating if the observed rating is in the training or validation set, then moves on to the next movie
      if (_ratings.r(n,m) > 0 || is_validation(r)) { // skip training and validation
        mlist[m].first = m;
        mlist[m].second = .0;
        ndcglist[m].first = m;
        ndcglist[m].second = 0;
//                cout << "nonzero" << endl;
        continue;
      }
      
      // Saves the predicted score in u
      double u = .0;
      if (_env.nmf) {
//        cout << "nmf" << endl;
        u = prediction_score_nmf(n, m);
      } else if (_env.lda || _env.vwlda) {
//                cout << "lda vwlda" << endl;
        u = prediction_score_lda(n, m);
      } else if (_env.graphchi) {
//                cout << "graphchi" << endl;
        u = prediction_score_chi(n, m);
      } else if (_env.ctr) {
//                cout << "ctr" << endl;
        u = prediction_score_ctr(n, m);
      } else {
//                cout << "else" << endl;
        u = _env.hier ? prediction_score_hier(n, m) : prediction_score(n, m);
      }
      
      // Saves the predicted rating in mlist. If the rating is in the test set, saves the observed rating in ndcglist
      mlist[m].first = m;
      mlist[m].second = u;
      ndcglist[m].first = m;
      CountMap::const_iterator itr = _ratings._test_map.find(r);
      if (itr != _ratings._test_map.end()) {
        ndcglist[m].second = itr->second;
      } else {
        ndcglist[m].second = 0;
      }
    } // End of loop over items
    
    // Variables to store the number of hits
    uint32_t hits10 = 0, hits100 = 0;
    double   dcg10 = .0, dcg100 = .0;
    
    // Sorts predicted ratings
    mlist.sort_by_value();
    
    // Loops again over items
    for (uint32_t j = 0; j < mlist.size() && j < _topN_by_user; ++j) {
      // Saves the user's predicted rating for the current item as (m,pred)
      KV &kv = mlist[j];
      uint32_t m = kv.first;
      double pred = kv.second;
      Rating r(n, m);
      
      // If save_ranking_file, saves the codes of the user and the item as m2 and n2
      uint32_t m2 = 0, n2 = 0;
      if (save_ranking_file) {
        IDMap::const_iterator it = _ratings.seq2user().find(n);
        assert (it != _ratings.seq2user().end());
        
        IDMap::const_iterator mt = _ratings.seq2movie().find(m);
        if (mt == _ratings.seq2movie().end())
          continue;
        
        m2 = mt->second;
        n2 = it->second;
      }
      
      // If the movie is in the test set
      CountMap::const_iterator itr = _ratings._test_map.find(r);
      if (itr != _ratings._test_map.end()) {
        // Saves the rating from the test set as v_
        int v_ = itr->second;
        
        // v is whether the rating is above the threshold for a hit
        int v;
        if (_ratings.test_hit(v_))
          v = 1;
        else
          v = 0;
        
        // Counts the number of hits in the best ranked 10 and 100 items
        if (j < 10) {
          if (v > 0) { //hit
            hits10++;
            hits100++;
          }
        } else if (j < 100) {
          if (v > 0)
            hits100++;
        }
        
        // Saves the codes, the predicted value, and whether it is a hit
        if (save_ranking_file) {
          if (_ratings.r(n, m) == .0)  {
            //double hol = _env.hier ? rating_likelihood_hier(n,m,v) : rating_likelihood(n,m,v);
            //fprintf(f, "%d\t%d\t%.5f\t%d\t%.5f\n", n2, m2, pred, v,
            //(pow(2.,v_) - 1)/log(j+2));
            fprintf(f, "%d\t%d\t%.5f\t%d\n", n2, m2, pred, v);
          }
        }
      } else { // If not in the test set, just save the codes, the predicted value, and whether it is a hit
        if (save_ranking_file) {
          if (_ratings.r(n, m) == .0) {
            //double hol = _env.hier ? rating_likelihood_hier(n,m,0) : rating_likelihood(n,m,0);
            //fprintf(f, "%d\t%d\t%.5f\t%d\t%.5f\n", n2, m2, pred, 0, .0);
            fprintf(f, "%d\t%d\t%.5f\t%d\n", n2, m2, pred, 0);
          }
        }
      }
    }
    
    // Finds fraction of hits among best 10 and 100 ranked items
    mhits10 += (double)hits10 / 10;
    mhits100 += (double)hits100 / 100;
    total_users++;
    // DCG normalizer
    double dcg10_gt = 0, dcg100_gt = 0;
    bool user_has_test_ratings = true;
    ndcglist.sort_by_value();
    for (uint32_t j = 0; j < ndcglist.size() && j < _topN_by_user; ++j) {
      int v = ndcglist[j].second;
      if(v==0) { //all subsequent docs are irrelevant
        if(j==0)
          user_has_test_ratings = false;
        break;
      }
    }
  }
  
  // Saves average number of hits among best ranked 10 and 100 items
  if (save_ranking_file)
    fclose(f);
  fprintf(_pf, "%d\t%.5f\t%.5f\n", 
          total_users,
          (double)mhits10 / total_users, 
          (double)mhits100 / total_users);
  fflush(_pf);
  
  //fprintf(_df, "%.5f\t%.5f\n", 
  //cumndcg10 / total_users, 
  //cumndcg100 / total_users);
  fflush(_df);
}

double
HGAPRec::prediction_score(uint32_t user, uint32_t movie) const
{
  const double **etheta = _theta.expected_v().const_data();
  const double **ebeta = _beta.expected_v().const_data();
  const double **td = _theta_mle.data();
  const double **bd = _beta_mle.data();

  double s = .0;
  for (uint32_t k = 0; k < _k; ++k)
    if (_env.mle_user)
      s += td[user][k] * ebeta[movie][k];
    else if (_env.mle_item || _env.canny)
      s += etheta[user][k] * bd[movie][k];
    else
      s += etheta[user][k] * ebeta[movie][k];

  if (_env.bias) {
    const double **ethetabias = _thetabias.expected_v().const_data();
    const double **ebetabias = _betabias.expected_v().const_data();
    s += ethetabias[user][0] + ebetabias[movie][0];
  }
  
  if (_use_rate_as_score)
    return s;
  
  if (s < 1e-30)
    s = 1e-30;
  double prob_zero = exp(-s);
  return 1 - prob_zero;
}

double
HGAPRec::prediction_score_nmf(uint32_t user, uint32_t movie) const
{
  const double **etheta = _nmf_theta->const_data();
  const double **ebeta = _nmf_beta->const_data();
  double s = .0;
  for (uint32_t k = 0; k < _k; ++k)
    s += etheta[user][k] * ebeta[movie][k];
  return s;
}

double
HGAPRec::prediction_score_ctr(uint32_t user, uint32_t movie) const
{
  const double **etheta = _ctr_theta->const_data();
  const double **ebeta = _ctr_beta->const_data();

  IDMap::const_iterator it = _ratings.seq2user().find(user);
  assert (it != _ratings.seq2user().end());
	
  IDMap::const_iterator mt = _ratings.seq2movie().find(movie);
  assert (mt != _ratings.seq2movie().end());

  uint32_t n = it->second;      
  uint32_t m = mt->second;

  IDMap::const_iterator it2 = _ctr_user_to_idx.find(n);
  assert (it2 != _ctr_user_to_idx.end());
	
  IDMap::const_iterator mt2 = _ctr_item_to_idx.find(m);
  assert (mt2 != _ctr_item_to_idx.end());

  uint32_t n2 = it2->second;
  uint32_t m2 = mt2->second;

  if (n2 >= _ctr_theta->m() || m2 >= _ctr_beta->m()) {
    lerr("asking score for nonexistent user or item: (%d,%d)", n2, m2);
    return .0;
  }

  double s = .0;
  for (uint32_t k = 0; k < _k; ++k)
    s += etheta[n2][k] * ebeta[m2][k];
  return s;
}

double
HGAPRec::prediction_score_chi(uint32_t user, uint32_t movie) const
{
  const double **etheta = _chi_gamma->const_data();
  const double **ebeta = _chi_beta->const_data();
  if (user == 0)
    debug("user = %d, movie = %d, rating = %d", user, movie, _ratings.r(user, movie));

  double s = .0;
  for (uint32_t k = 0; k < _k; ++k) {
    s += etheta[user][k] * ebeta[movie][k];
    if (user == 0)
      debug("%f, %f", etheta[user][k], ebeta[movie][k]);
  }
  if (user == 0)
    debug("s = %f", s);

  if (_env.bias) {
    const double **ubias = _chi_ubias->const_data();
    const double **vbias = _chi_vbias->const_data();
    const double **global = _chi_global->const_data();
    s += ubias[user][0] + vbias[movie][0];
    debug("%f\t%f\n", ubias[user][0], vbias[movie][0]);
    s += global[0][0];
    debug("%f\n", global[0][0]);
  }
  return s;
}

double
HGAPRec::prediction_score_lda(uint32_t user, uint32_t movie) const
{
  const double **etheta = _lda_gamma->const_data();
  const double **ebeta = _lda_beta->const_data();
  double s = .0;
  for (uint32_t k = 0; k < _k; ++k)
    s += etheta[user][k] * ebeta[k][movie];
  return s;
}


// Returns the predicted score. If _use_rate_as_score, returns the rate of the Poisson random variable, otherwise returns the probability of a nonzero rating
double
HGAPRec::prediction_score_hier(uint32_t user, uint32_t movie) const
{
  // Saves the expected values of theta and beta
  const double **etheta = _htheta.expected_v().const_data();
  const double **ebeta = _hbeta.expected_v().const_data();
  
  // Computes the dot product of the expected values of theta and beta (the rate of the Poisson random variable)
  double s = .0;
  for (uint32_t k = 0; k < _k; ++k)
    s += etheta[user][k] * ebeta[movie][k];

  if (_env.bias) {
    const double **ethetabias = _thetabias.expected_v().const_data();
    const double **ebetabias = _betabias.expected_v().const_data();
    s += ethetabias[user][0] + ebetabias[movie][0];
  }
  
  // Returns the expected value of the Poisson rv (its rate)
  if (_use_rate_as_score) {
//    cout << "Use rate as score" << endl;
    return s;
  }
  
  // Returns the probability of a nonzero rating
  if (s < 1e-30)
    s = 1e-30;
  double prob_zero = exp(-s);
  return 1 - prob_zero;
}

void
HGAPRec::gen_msr_csv()
{
  load_beta_and_theta();

  FILE *f = fopen(Env::file_str("/pred.csv").c_str(), "w");
  if (!f) {
    lerr("cannot open pred.csv");
    return;
  }
  
  fprintf(f, "User\tHeldOutItem\tHeldOutItemIndex\tUserNegatives\tUserCount\tItemCount\n");
  for (uint32_t n = 0; n < _n; ++n) {
    // User
    IDMap::const_iterator it = _ratings.seq2user().find(n);
    assert (it != _ratings.seq2user().end());
    fprintf(f, "%d\t", it->second);
    debug("user: %d", it->second);

    // id of heldout item
    IDMap::const_iterator ct = _ratings.leave_one_out().find(n);
    debug("heldout item for user %d (%d)",   it->second, n);
    assert (ct != _ratings.leave_one_out().end());

    uint32_t test_item_seq = ct->second;
    IDMap::const_iterator pt = _ratings.seq2movie().find(test_item_seq);
    assert (pt != _ratings.seq2movie().end());
    fprintf(f, "%d\t", pt->second);
    debug("test item: %d", pt->second);

    // rank of heldout test item
    uint32_t training = 0, negatives = 0;
    KVArray mlist(_m);
    for (uint32_t m = 0; m < _m-1; ++m) {
      Rating r(n, m);
      if (_ratings.r(n,m) > 0 || is_validation(r)) { // skip training non-zero rating
	mlist[m].first = m;
	mlist[m].second = .0;
	training++;
	continue;
      }
      double u = .0;
      if (_env.nmf)
	u = prediction_score_nmf(n, m);
      else
	u = _env.hier ? prediction_score_hier(n, m) : prediction_score(n, m);
      mlist[m].first = m;
      mlist[m].second = u;
      negatives++;
    }
    mlist.sort_by_value();

    uint32_t c = 0, rank = 0;
    for (uint32_t j = 0; j < mlist.size(); ++j) {
      KV &kv = mlist[j];
      uint32_t m = kv.first;
      double pred = kv.second;
      if (m == test_item_seq) {
	rank = c;
	break;
      }
      c++;
    }
    fprintf(f, "%d\t", rank);
    debug("rank: %d", rank);
    
    // UserNegatives
    fprintf(f, "%d\t", negatives);
    debug("user negatives: %d", negatives);

    // UserCount
    fprintf(f, "%d\t", training);
    debug("user count: %d", training);
    
    const vector<uint32_t> *q = _ratings.get_users(test_item_seq);
    uint32_t ntraining_users = 0;
    if (q)
      ntraining_users = q->size();
    
    uint32_t nvalid_users = 0;
    FreqMap::const_iterator itr = _ratings.validation_users_of_movie().find(test_item_seq);
    if (itr != _ratings.validation_users_of_movie().end())
      nvalid_users = itr->second;
    
    // ItemCount
    fprintf(f, "%d\n", nvalid_users + ntraining_users);
    debug("item count: %d", nvalid_users + ntraining_users);

    printf("\r user %d", n);
    fflush(stdout);
  }
  
}

void
HGAPRec::gen_ranking_for_users(bool load)
{
  if (_env.ctr) {
    lerr("loading CTR files");
    load_ctr_beta_and_theta();
  } else if (load)
    load_beta_and_theta();

  char buf[4096];
  sprintf(buf, "%s/test_users.tsv", _env.datfname.c_str());
  lerr("loading test users from %s", buf);
  FILE *f = fopen(buf, "r");
  if (!f) { 
    lerr("cannot open %s", buf);
    return;
  }
  //assert(f);
  _sampled_users.clear();
  _ratings.read_test_users(f, &_sampled_users);
  fclose(f);
  
  compute_precision(true);
  compute_itemrank(true);
  lerr("DONE writing ranking.tsv in output directory\n");
}

void
HGAPRec::load_beta_and_theta()
{
  if (!_env.hier) {
    _beta.load();
    _theta.load();
  } else {
    _thetarate.load();
    _betarate.load();

    _hbeta.load();
    _htheta.load();
  }
  if (_env.bias) {
    _betabias.load();
    _thetabias.load();
  }
  if (_env.canny || _env.mle_user || _env.mle_item) {
    _theta_mle.load("theta_mle.tsv");
    _beta_mle.load("beta_mle.tsv");
  }
}

// Saves the matrices in output files
void
HGAPRec::save_model()
{
  if (_env.hier) {
    _hbeta.save_state(_ratings.seq2movie(),_env.outfname);
    _betarate.save_state(_ratings.seq2movie(),_env.outfname);
    _htheta.save_state(_ratings.seq2user(),_env.outfname);
    _thetarate.save_state(_ratings.seq2user(),_env.outfname);
    _hsigma.save_state(_ratings.seq2user(),_env.outfname);
    _hrho.save_state(_ratings.seq2movie(),_env.outfname);
  } else {
    _beta.save_state(_ratings.seq2movie(),_env.outfname);
    _theta.save_state(_ratings.seq2user(),_env.outfname);
  }

  if (_env.bias) {
    _betabias.save_state(_ratings.seq2movie(),_env.outfname);
    _thetabias.save_state(_ratings.seq2user(),_env.outfname);
  }
  if (_env.canny || _env.mle_user || _env.mle_item) {
    _theta_mle.save(Env::file_str("/theta_mle.tsv"), _ratings.seq2user());
    _beta_mle.save(Env::file_str("/beta_mle.tsv"), _ratings.seq2movie());
  }
  
//  save_phi();
}

//void
//HGAPRec::save_phi() {
//  string name = string("/phi.tsv");
//  phi.save(_env.outfname+"/"+Env::outfile_str(name));
//}

void
HGAPRec::logl()
{
  uint32_t x;
  if (_env.bias)
    x = _k+2;
  else
    x = _k;

  Array phi(x);
  double s = .0;

  const double  **etheta = NULL;
  const double  **ebeta = NULL; 
  const double  **elogtheta = NULL; 
  const double  **elogbeta = NULL; 
  const double  **eu = NULL;
  const double  **ei = NULL; 
  const double  **elogu = NULL; 
  const double  **elogi = NULL; 

  if (!_env.hier) {
    etheta = _theta.expected_v().const_data();
    ebeta = _beta.expected_v().const_data();
    elogtheta = _theta.expected_logv().const_data();
    elogbeta = _beta.expected_logv().const_data();
  } else {
    etheta = _htheta.expected_v().const_data();
    ebeta = _hbeta.expected_v().const_data();
    elogtheta = _htheta.expected_logv().const_data();
    elogbeta = _hbeta.expected_logv().const_data();
  }

  if (_env.bias) {
    eu = _thetabias.expected_v().const_data();
    ei = _betabias.expected_v().const_data();
    elogu = _thetabias.expected_logv().const_data();
    elogi = _betabias.expected_logv().const_data();
  }

  for (uint32_t n = 0; n < _n; ++n) {
    const vector<uint32_t> *movies = _ratings.get_movies(n);
    for (uint32_t j = 0; j < movies->size(); ++j) {
      uint32_t m = (*movies)[j];
      yval_t y = _ratings.r(n,m);
      
      if (_env.hier) {
	if (!_env.bias) 
	  get_phi(_htheta, n, _hbeta, m, phi);
	else
	  get_phi(_htheta, n, _hbeta, m, elogu[n][0], elogi[m][0], phi);
      } else {
	if (!_env.bias)
	  get_phi(_theta, n, _beta, m, phi);
	else
	  get_phi(_theta, n, _beta, m, elogu[n][0], elogi[m][0], phi);
      }
      if (y > 1)
	phi.scale(y);

      double v = .0;
      for (uint32_t k = 0; k < _k; ++k)
	v += y * phi[k] * (elogtheta[n][k] + elogbeta[m][k] - log(phi[k]));
      s += v;
      if (_env.bias) {
	s += y * phi[_k] * (elogu[n][0] - log(phi[_k]));
	s += y * phi[_k+1] * (elogi[m][0] - log(phi[_k+1]));
      }
      
      for (uint32_t k = 0; k < _k; ++k)
	s -= etheta[n][k] * ebeta[m][k];
      if (_env.bias) {
	s -= eu[n][0];
	s -= ei[m][0];
      }
    }
  }

  if (!_env.hier) {
    s += _theta.compute_elbo_term();
    s += _beta.compute_elbo_term();
  } else {
    s += _htheta.compute_elbo_term();
    s += _hbeta.compute_elbo_term();
    s += _thetarate.compute_elbo_term();
    s += _betarate.compute_elbo_term();
  }

  if (_env.bias) {
    s += _thetabias.compute_elbo_term();
    s += _betabias.compute_elbo_term();
  }
    
  fprintf(_af, "%.5f\n", s);
  fflush(_af);
}

void
HGAPRec::test()
{
  _ratings.read_netflix_metadata("");
  //_ratings.read_movielens_metadata("");

  vector<uint32_t> items;
  items.push_back(118);
  items.push_back(12263);

  for (uint32_t i = 0; i < items.size(); ++i) {
    const IDMap &m = _ratings.movie2seq();
    IDMap::const_iterator itr = m.find(items[i]);
    assert(itr != m.end());
    items[i] = itr->second;
    lerr("%d\n", items[i]); 
    printf("%s\n", _ratings.movie_name(items[i]).c_str());
  }

  load_beta_and_theta();
  printf("loading model state complete\n");
  _theta.set_to_prior();
  _theta.set_to_prior_curr();

  //items.push_back(2403);
  //items.push_back(2404);

  //items.push_back(2398);
  //items.push_back(2384);

  //items.push_back(1091);
  //items.push_back(1080);
  //items.push_back(1107);
  //  items.push_back(1125);

  // items.push_back(2687);
  // items.push_back(2709);
  // items.push_back(3396);
  // items.push_back(3398);
  // items.push_back(3399);

  Array betarowsum(_k);
  _beta.sum_rows(betarowsum);

  printf("predictions\n");
  uint32_t user = 0;
  for (uint32_t itr = 0; itr < 10; itr++) {
    const double  **eloga = _theta.expected_logv().const_data();
    const double  **elogb = _beta.expected_logv().const_data();
    
    for (uint32_t i = 0; i < items.size(); ++i) {
      Array phi(_k);
      phi.zero();
      for (uint32_t k = 0; k < _k; ++k)
	phi[k] = eloga[user][k] + elogb[items[i]][k];
      phi.lognormalize();
      debug("%s", phi.s().c_str());
      _theta.update_shape_next(user, phi);
    }
    _theta.update_rate_next(betarowsum);
    _theta.swap();
    _theta.compute_expectations();
    printf("iteration %d\n", itr);
  }

  // predictions for user
  KVArray mlist(_m);
  for (uint32_t m = 0; m < _m; ++m) {
    Rating r(user,m);
    double u = prediction_score(user, m);
    for (uint32_t i = 0; i < items.size(); ++i) {
      if (items[i] == m)
	continue;
    }
    mlist[m].first = m;
    mlist[m].second = u;
  }

  mlist.sort_by_value();
  uint32_t t =0;
  uint32_t c = 0, rank = 0;
  for (uint32_t j = 0; j < mlist.size(); ++j,t++) {
    if (t > 20)
      break;
    KV &kv = mlist[j];
    uint32_t m = kv.first;
    printf("%s, %s\n", _ratings.movie_name(m).c_str(), _ratings.movie_type(m).c_str());
    lerr("%s, %s\n", _ratings.movie_name(m).c_str(), _ratings.movie_type(m).c_str());
  }
}


// Methods that are not used in the paper
/*
void
HGAPRec::vb()
{
  lerr("running vb()");
  //cout << "running vb()" << endl;
  
  // Initial values of the paremters
  initialize();
  // approx_log_likelihood();
  
  // Initializes the array of multinomial parameters
  // Array is a one dimensional array of doubles of dimension k
  Array phi(_k);
  
  // Optimization loop
  while (1) {
    for (uint32_t n = 0; n < _n; ++n) {
      
      // Pointer with the
      const vector<uint32_t> *movies = _ratings.get_movies(n);
      for (uint32_t j = 0; j < movies->size(); ++j) {
        uint32_t m = (*movies)[j];
        yval_t y = _ratings.r(n,m);
        
        get_phi(_theta, n, _beta, m, phi);
        if (y > 1)
          phi.scale(y);
        
        _theta.update_shape_next(n, phi);
        _beta.update_shape_next(m, phi);
      }
    }
    
    Array betasum(_k);
    _beta.sum_rows(betasum);
    _theta.update_rate_next(betasum);
    
    _theta.swap();
    _theta.compute_expectations();
    
    Array thetasum(_k);
    _theta.sum_rows(thetasum);
    _beta.update_rate_next(thetasum);
    
    _beta.swap();
    _beta.compute_expectations();
    
    printf("\r iteration %d", _iter);
    fflush(stdout);
    if (_iter % _env.reportfreq == 0) {
      // approx_log_likelihood();
      compute_likelihood(true);
      compute_likelihood(false);
      //compute_rmse();
      save_model();
      compute_precision(false);
      compute_itemrank(false);
      //gen_ranking_for_users(false);
      if (_env.logl)
        logl();
    }
    
    if (_env.save_state_now) {
      lerr("Saving state at iteration %d duration %d secs", _iter, duration());
      do_on_stop();
    }
    
    _iter++;
  }
}

void
HGAPRec::vb_mle_user()
{
  initialize();
  
  double **td = _theta_mle.data();
  double **old_td = _old_theta_mle.data();
  
  for (uint32_t n = 0; n < _n; ++n)
    for (uint32_t k = 0; k < _k; ++k)
      old_td[n][k] = (double)1.0/(double)_k;
  
  Array phi(_k);
  while (1) {
    for (uint32_t n = 0; n < _n; ++n) {
      for (uint32_t k = 0; k < _k; ++k)
        td[n][k] = 0;
      
      const vector<uint32_t> *movies = _ratings.get_movies(n);
      for (uint32_t j = 0; j < movies->size(); ++j) {
        uint32_t m = (*movies)[j];
        yval_t y = _ratings.r(n,m);
        
        get_phi(_old_theta_mle, n, _beta, m, phi);
        if (y > 1)
          phi.scale(y);
        
        for (uint32_t k = 0; k < _k; ++k)
          td[n][k] += phi[k];
        
        _beta.update_shape_next(m, phi);
      }
    }
    
    Array betasum(_k);
    _beta.sum_rows(betasum);
    
    double sum = .0;
    for (uint32_t n = 0; n < _n; ++n)
      for (uint32_t k = 0; k < _k; ++k)  {
        td[n][k] /= betasum[k];
        old_td[n][k] = td[n][k];
      }
    
    Array thetasum(_k);
    for (uint32_t n = 0; n < _n; ++n)
      for (uint32_t k = 0; k < _k; ++k)
        thetasum[k] += td[n][k];
    
    _beta.update_rate_next(thetasum);
    _beta.swap();
    _beta.compute_expectations();
    
    printf("\r iteration %d", _iter);
    fflush(stdout);
    if (_iter % _env.reportfreq == 0) {
      // approx_log_likelihood();
      compute_likelihood(true);
      compute_likelihood(false);
      //compute_rmse();
      //save_model();
      compute_precision(false);
      //gen_ranking_for_users(false);
      if (_env.logl)
        logl();
    }
    
    if (_env.save_state_now) {
      lerr("Saving state at iteration %d duration %d secs", _iter, duration());
      do_on_stop();
    }
    
    _iter++;
  }
}


void
HGAPRec::vb_mle_item()
{
  initialize();
  
  double **bd = _beta_mle.data();
  double **old_bd = _old_beta_mle.data();
  
  for (uint32_t m = 0; m < _m; ++m)
    for (uint32_t k = 0; k < _k; ++k)
      old_bd[m][k] = (double)1.0/(double)_k;
  
  Array phi(_k);
  while (1) {
    for (uint32_t m = 0; m < _m; ++m)
      for (uint32_t k = 0; k < _k; ++k)
        bd[m][k] = 0;
    
    for (uint32_t n = 0; n < _n; ++n) {
      const vector<uint32_t> *movies = _ratings.get_movies(n);
      for (uint32_t j = 0; j < movies->size(); ++j) {
        uint32_t m = (*movies)[j];
        yval_t y = _ratings.r(n,m);
        
        get_phi(_theta, n, _old_beta_mle, m, phi);
        if (y > 1)
          phi.scale(y);
        
        for (uint32_t k = 0; k < _k; ++k)
          bd[m][k] += phi[k];
        
        _theta.update_shape_next(n, phi);
      }
    }
    
    Array thetasum(_k);
    _theta.sum_rows(thetasum);
    
    double sum = .0;
    for (uint32_t m = 0; m < _m; ++m)
      for (uint32_t k = 0; k < _k; ++k)  {
        bd[m][k] /= thetasum[k];
        old_bd[m][k] = bd[m][k];
      }
    
    Array betasum(_k);
    for (uint32_t m = 0; m < _m; ++m)
      for (uint32_t k = 0; k < _k; ++k)
        betasum[k] += bd[m][k];
    
    _theta.update_rate_next(betasum);
    _theta.swap();
    _theta.compute_expectations();
    
    printf("\r iteration %d", _iter);
    fflush(stdout);
    if (_iter % _env.reportfreq == 0) {
      compute_likelihood(true);
      compute_likelihood(false);
      compute_precision(false);
      if (_env.logl)
        logl();
    }
    
    if (_env.save_state_now) {
      lerr("Saving state at iteration %d duration %d secs", _iter, duration());
      do_on_stop();
    }
    
    _iter++;
  }
}

void
HGAPRec::vb_canny()
{
  initialize();
  
  double **bd = _beta_mle.data();
  double **old_bd = _old_beta_mle.data();
  
  for (uint32_t m = 0; m < _m; ++m)
    for (uint32_t k = 0; k < _k; ++k)
      old_bd[m][k] = (double)1.0/(double)_m;
  
  Array phi(_k);
  while (1) {
    for (uint32_t m = 0; m < _m; ++m)
      for (uint32_t k = 0; k < _k; ++k)
        bd[m][k] = 0;
    
    for (uint32_t n = 0; n < _n; ++n) {
      const vector<uint32_t> *movies = _ratings.get_movies(n);
      for (uint32_t j = 0; j < movies->size(); ++j) {
        uint32_t m = (*movies)[j];
        yval_t y = _ratings.r(n,m);
        
        get_phi(_theta, n, _old_beta_mle, m, phi);
        if (y > 1)
          phi.scale(y);
        
        for (uint32_t k = 0; k < _k; ++k)
          bd[m][k] += phi[k];
        
        _theta.update_shape_next(n, phi);
      }
    }
    
    Array thetasum(_k);
    _theta.sum_rows(thetasum);
    
    double sum = .0;
    for (uint32_t m = 0; m < _m; ++m)
      for (uint32_t k = 0; k < _k; ++k)  {
        if (bd[m][k] < 1e-10)
          bd[m][k] = 1e-10;
        bd[m][k] = bd[m][k] / thetasum[k];
        old_bd[m][k] = bd[m][k];
      }
    
    for (uint32_t k = 0; k < _k; ++k) {
      double s = .0;
      for (uint32_t m = 0; m < _m; ++m)
        s += old_bd[m][k];
      for (uint32_t m = 0; m < _m; ++m) {
        old_bd[m][k] /= s;
        bd[m][k] = old_bd[m][k];
      }
    }
    
    Array betasum(_k);
    for (uint32_t m = 0; m < _m; ++m)
      for (uint32_t k = 0; k < _k; ++k)
        betasum[k] += bd[m][k];
    
    _theta.update_rate_next(betasum);
    _theta.swap();
    _theta.compute_expectations();
    
    printf("\r iteration %d", _iter);
    fflush(stdout);
    if (_iter % _env.reportfreq == 0) {
      compute_likelihood(true);
      compute_likelihood(false);
      compute_precision(false);
      compute_itemrank(false);
      if (_env.logl)
        logl();
    }
    
    if (_env.save_state_now) {
      lerr("Saving state at iteration %d duration %d secs", _iter, duration());
      do_on_stop();
    }
    
    _iter++;
  }
}


void
HGAPRec::vb_bias()
{
  lerr("running vb_bias()");
  initialize();
  
  Array phi(_k+2);
  while (1) {
    for (uint32_t n = 0; n < _n; ++n) {
      
      const vector<uint32_t> *movies = _ratings.get_movies(n);
      for (uint32_t j = 0; j < movies->size(); ++j) {
        uint32_t m = (*movies)[j];
        yval_t y = _ratings.r(n,m);
        
        const double **tbias = _thetabias.expected_logv().const_data();
        const double **bbias = _betabias.expected_logv().const_data();
        
        get_phi(_theta, n, _beta, m, tbias[n][0], bbias[m][0], phi);
        
        if (y > 1)
          phi.scale(y);
        
        _theta.update_shape_next(n, phi);
        _beta.update_shape_next(m, phi);
        
        _thetabias.update_shape_next3(n, 0, phi[_k]);
        _betabias.update_shape_next3(m, 0, phi[_k+1]);
      }
    }
    
    if (_env.vb) {
      Array betasum(_k);
      _beta.sum_rows(betasum);
      _theta.update_rate_next(betasum);
      
      _theta.swap();
      _theta.compute_expectations();
      
      Array thetasum(_k);
      _theta.sum_rows(thetasum);
      _beta.update_rate_next(thetasum);
      
      _beta.swap();
      _beta.compute_expectations();
      
      _thetabias.update_rate_next_all(0, _m);
      _thetabias.swap();
      _thetabias.compute_expectations();
      
      _betabias.update_rate_next_all(0, _n);
      _betabias.swap();
      _betabias.compute_expectations();
      
      debug("thetabias: %s", _thetabias.expected_v().s().c_str());
      debug("betabias: %s", _betabias.expected_v().s().c_str());
      
    } else {
      
      Array betasum(_k);
      _beta.sum_rows(betasum);
      _theta.update_rate_next(betasum);
      Array thetasum(_k);
      _theta.sum_rows(thetasum);
      _beta.update_rate_next(thetasum);
      
      _thetabias.update_rate_next_all(0, _m);
      _betabias.update_rate_next_all(0, _n);
      
      _theta.swap();
      _beta.swap();
      _thetabias.swap();
      _betabias.swap();
      
      _theta.compute_expectations();
      _beta.compute_expectations();
      _thetabias.compute_expectations();
      _betabias.compute_expectations();
    }
    
    printf("\r iteration %d", _iter);
    fflush(stdout);    
    if (_iter % _env.reportfreq == 0) {
      compute_likelihood(true);
      compute_likelihood(false);
      //compute_rmse();
      save_model();
      compute_precision(false);
      //gen_ranking_for_users(false);
      if (_env.logl)
        logl();
    }
    
    if (_env.save_state_now) {
      lerr("Saving state at iteration %d duration %d secs", _iter, duration());
      do_on_stop();
    }
    
    _iter++;
  }
}
*/
