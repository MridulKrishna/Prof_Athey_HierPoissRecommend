#include "ratings.hh"
#include "log.hh"
#include <wchar.h>
#include <iostream>
#include <string>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <set>

using namespace std;
using namespace boost;

// Reads dataset
int
Ratings::read(string s)
{
  fprintf(stdout, "+ reading ratings dataset from %s\n", s.c_str());
  fflush(stdout);

    if (_env.mode == Env::CREATE_TRAIN_TEST_SETS) {
        if (_env.dataset == Env::NETFLIX) {
            for (uint32_t i = 0; i < _env.m; ++i) {
                if (read_netflix_movie(s,i+1) < 0) {
                    lerr("error adding movie %d\n", i);
                    return -1;
                }
            }
        } else if (_env.dataset == Env::MOVIELENS)
            read_movielens(s);
        else  if (_env.dataset == Env::MENDELEY)
            read_mendeley(s);
        else if (_env.dataset == Env::ECHONEST)
            read_echonest(s);
        else if (_env.dataset == Env::NYT)
            read_nyt(s);
    } else {
        // Runs this line to read the file
        read_generic_train(s);
        write_marginal_distributions();
    }
    
  char st[1024];
  sprintf(st, "read %d users, %d movies, %d ratings", 
	  _curr_user_seq, _curr_movie_seq, _nratings);
  _env.nTrain = _curr_user_seq;
  _env.mTrain = _curr_movie_seq;
  Env::plog("statistics", string(st));

  return 0;
}

// Reads datasets with validation and test sets
void
Ratings::readValidationAndTest(string dir)
{
  
  char buf[4096];
  sprintf(buf, "%s/validation.tsv", _env.datfname.c_str());
  FILE *validf = fopen(buf, "r");
  assert(validf);
  if (_env.dataset == Env::NYT)
    read_nyt_train(validf, &_validation_map);
  else
    // Saves the ratings from the validation file in validation_map
    read_generic(validf, &_validation_map);
  fclose(validf);
  
  // Loop with iterator on the elements saved in _validation_map
  // Builds a histogram with the number of users per movie in _validation_users_of_movie
  for (CountMap::const_iterator i = _validation_map.begin();
       i != _validation_map.end(); ++i) {
    const Rating &r = i->first;
    const Rating r2 = r;
    _validation_users_of_movie[r.second]++;
  }
  
  sprintf(buf, "%s/test.tsv", _env.datfname.c_str());
  FILE *testf = fopen(buf, "r");
  assert(testf);
  if (_env.dataset == Env::NYT)
    read_nyt_train(testf, &_test_map);
  else
    // Saves the ratings from the validation file in test_map
    read_generic(testf, &_test_map);
  fclose(testf);
  
  // XXX: keeps one heldout test item for each user
  // assumes leave-one-out
  // JCC: Loop with iterator on the elements saved in _validation_map
  // Builds a histogram with the number of users per movie in _validation_users_of_movie
  for (CountMap::const_iterator i = _test_map.begin();
       i != _test_map.end(); ++i) {
    const Rating &r = i->first;
    _leave_one_out[r.first] = r.second;
    debug("adding %d -> %d to leave one out", r.first, r.second);
  }
  
  printf("+ loaded validation and test sets from %s\n", _env.datfname.c_str());
  fflush(stdout);
  Env::plog("test ratings", _test_map.size());
  Env::plog("validation ratings", _validation_map.size());
}

// Reads dataset with observed user and item characteristics
void
Ratings::readObserved(string dir)
{
  // Message if there are actually some observed characteristics
  if (_env.uc > 0 || _env.ic > 0) {
    fprintf(stdout, "+ reading observed characteristics from %s\n", dir.c_str());
    fflush(stdout);
  }
  
  // Read user characteristics if the number is nonzero
  if (_env.uc > 0) {
    char buf[1024];
    sprintf(buf, "%s/obsUser.tsv", dir.c_str());

    ifstream infile(buf);
    string line;
    int nLines = 0;
    set<uint64_t> users;
    
    // Loops over the lines of the file and saves the variables
    while ( getline(infile,line)) {
      nLines++;
      assert(nLines <= _env.n);
      vector<string> strs;
      split(strs,line,is_any_of("\t"));
      assert(size(strs)==_env.uc+1);
      
      // Gets the code of the user in that line
      uint64_t uCode = stoi(strs.at(0));
      // Checks that the user is in the list of users
//      cout << uCode << endl;
      assert(_user2seq.find(uCode) != _user2seq.end());
      // Checks that this is the first line for that user
      assert(users.insert(uCode).second);
      uint64_t pos = _user2seq[uCode];
      
//      cout << uCode << " " << pos << endl;
      
      Array row(_env.uc,true);

      // Loops over the variables in the line and saves each one in the matrix of observed characteristics
      for ( int i = 0; i<_env.uc; i++) {
//        _userObs.get(pos,i) = stod(strs.at(i+1));
        _userObs.set(pos,i,stod(strs.at(i+1)));
      }
    }
    assert(nLines==_env.n);
//    _userObs.print();
    
    // Compute vector of scale of observed user characteristics
    if (_env.scale == Env::MEAN) {
      _userObs.colmeans(_userObsScale);
      _userObsScale.scale(_env.scaleFactor);
    } else if (_env.scale == Env::ONES) {
      _userObsScale.set_elements(1);
      _userObsScale.scale(_env.scaleFactor);
    } else if (_env.scale == Env::STD) {
      _userObs.colstds(_userObsScale);
      _userObsScale.scale(_env.scaleFactor);
    } else {
      cout << "Invalid scale type" << endl;
      exit(0);
    }
//    _userObsScale.print();
  }
  
  // Read item characteristics if the number is nonzero
  if (_env.ic > 0) {
    char buf[1024];
    sprintf(buf, "%s/obsItem.tsv", dir.c_str());
    
    ifstream infile(buf);
    string line;
    int nLines = 0;
    set<uint64_t> items;
    
    // Loops over the lines of the file and saves the variables
    while ( getline(infile,line)) {
      nLines++;
      assert(nLines <= _env.m);
      vector<string> strs;
      split(strs,line,is_any_of("\t"));
      assert(size(strs)==_env.ic+1);
      
      // Gets the code of the item in that line
      uint64_t iCode = stol(strs.at(0));
//      cout << iCode << endl;
      // Checks that the user is in the list of users
      assert(_movie2seq.find(iCode) != _movie2seq.end());
      // Checks that this is the first line for that user
      assert(items.insert(iCode).second);
      uint64_t pos = _movie2seq[iCode];
      
//      cout << iCode << " " << pos << endl;
      
      Array row(_env.ic,true);
      
      // Loops over the variables in the line and saves each one in the matrix of observed characteristics
      for ( int i = 0; i<_env.ic; i++) {
        _itemObs.get(pos,i) = stod(strs.at(i+1));
      }
    }
    assert(nLines==_env.m);
//    _itemObs.print();
    
    // Compute vector of scale of observed item characteristics
    if (_env.scale == Env::MEAN) {
      _itemObs.colmeans(_itemObsScale);
      _itemObsScale.scale(_env.scaleFactor);
    } else if (_env.scale == Env::ONES) {
      _itemObsScale.set_elements(1);
      _itemObsScale.scale(_env.scaleFactor);
    } else if (_env.scale == Env::STD) {
      _itemObs.print();
      _itemObs.colstds(_itemObsScale);
      _itemObsScale.scale(_env.scaleFactor);
    } else {
      cout << "Invalid scale type" << endl;
      exit(0);
    }
//    _itemObsScale.print();
  }
//  _userObsScale.print();
//  _itemObsScale.print();
//  _userObsScale.print();
}

// Reads a generic train data file
void
Ratings::read_generic_train(string dir)
{
  char buf[1024];
  sprintf(buf, "%s/train.tsv", dir.c_str());
  
  //cout << "buf" << endl;
  //cout << buf << endl;
  
  FILE *f = fopen(buf, "r");
  if (!f) {
//    cout << "aqui" << endl;
        cout << buf << endl;
    fprintf(stderr, "error: cannot open file %s: %s", buf, strerror(errno));
    fclose(f);
    exit(-1);
  }

  if (_env.dataset == Env::NYT) {
    read_nyt_titles(dir);
    read_nyt_train(f, NULL);
  } else
    // Sends the FILE object pointer to read_generic
    read_generic(f, NULL);
  fclose(f);
  Env::plog("training ratings", _nratings);
}

// Reads a file of users, items, and ratings
// If cmap is NULL, save in the attributes of this. If it is not null, save in cmap
int
Ratings::read_generic(FILE *f, CountMap *cmap)
{
  assert(f);
  //char b[128];
  
  // Integers that will store the values read from each line for the user, item, and rating (mid because of movie?)
  uint64_t mid = 0, uid = 0;
  uint32_t rating;
  
  // Loop that reads each line in the file
  while (!feof(f)) {
    
    // Checks that the line has 3 numbers
    if (fscanf(f, "%llu\t%llu\t%u\n", &uid, &mid, &rating) < 0) {
      printf("error: unexpected lines in file\n");
      fclose(f);
      exit(-1);
    }
   
//    if ( uid == 2169292)
//      cout << uid << endl;
    
    IDMap::iterator it = _user2seq.find(uid);
    IDMap::iterator mt = _movie2seq.find(mid);
    
//    cout << "_curr_user_seq" << endl;
//    cout << _curr_user_seq << endl;
//    cout << "_env.n" << endl;
//    cout << _env.n << endl;
    
    // If all users and movies have been added, leave the loop
    // _curr_user_seq is the number of users added so far. Same for movies.
    if ((it == _user2seq.end() && _curr_user_seq >= _env.n) ||
        (mt == _movie2seq.end() && _curr_movie_seq >= _env.m))
      continue;
    
    // If the rating is a zero, do nothing
    if (input_rating_class(rating) == 0)
      continue;
    
    if (it == _user2seq.end()) {
//      cout << "Adding user" << endl;
      assert(add_user(uid));
    }
    
    if (mt == _movie2seq.end()) {
//            cout << "Adding movie" << endl;
      assert(add_movie(mid));
    }
    
    // Finds the indices for the user and the item
    uint64_t m = _movie2seq[mid];
    uint64_t n = _user2seq[uid];
    
    // If the value is not zero (otherwise do nothing)
    if (input_rating_class(rating) > 0) {
      // If cmap is NULL, i.e., if the ratings should be saved in _users2rating
      if (!cmap) {
        // Increases the counter of ratings
        _nratings++;
        
        // Pointer to the user's ratings
        RatingMap *rm = _users2rating[n];
        
        // Adds the new item to the user's rating
        if (_env.binary_data)
          (*rm)[m] = 1;
        else {
          assert (rating > 0);
          (*rm)[m] = rating;
        }
        
        // Adds the item to the list of items rated by the user, and the user to the list of users rating an item
        _users[n]->push_back(m);
        _movies[m]->push_back(n);
      } else {
        // If cmap is not NULL, i.e., if the ratings should be saved to cmap
        debug("adding test or validation entry for user %d, item %d", n, m);
        // Creates the index for the rating of class Rating (a pair of integers)
        Rating r(n,m);
        assert(cmap);
        
        // Saves the rating in the map
        if (_env.binary_data)
          (*cmap)[r] = 1;
        else
          (*cmap)[r] = rating;
      }
    }
  }
  return 0;
}

int
Ratings::read_nyt_titles(string dir)
{
  char buf[1024];
  sprintf(buf, "%s/nyt-titles.tsv", dir.c_str());
  FILE *f = fopen(buf, "r");
  if (!f) {
    fprintf(stderr, "error: cannot open file %s:%s", buf, strerror(errno));
    fclose(f);
    exit(-1);
  }
  char title[512];
  uint32_t id;
  uint32_t c = 0;
  char *line = (char *)malloc(4096);
  while (!feof(f)) {
    if (fgets(line, 4096, f) == NULL)
      break;
    char *p = line;
    const char r[3]="|";
    char *q = NULL;
    char *d=strtok_r(p, r, &q);
    id = atoi(d);
    IDMap::iterator mt = _movie2seq.find(id);    
    if (mt == _movie2seq.end()) {
      add_movie(id);
      c++;
    }
  }
  fclose(f);
  Env::plog("read titles", c);
  return 0;
}

int
Ratings::read_nyt_train(FILE *f, CountMap *cmap)
{
  assert(f);
  char b[128];
  uint32_t mid = 0, uid = 0, rating = 0;
  while (!feof(f)) {
    if (fscanf(f, "%u\t%u\t%u\n", &uid, &mid, &rating) < 0) {
      printf("error: unexpected lines in file\n");
      fclose(f);
      exit(-1);
    }

    IDMap::iterator it = _user2seq.find(uid);
    IDMap::iterator mt = _movie2seq.find(mid);

    if ((it == _user2seq.end() && _curr_user_seq >= _env.n) ||
	(mt == _movie2seq.end()))
      continue;

    if (input_rating_class(rating) == 0)
      continue;
    
    if (it == _user2seq.end())
      assert(add_user(uid));
    
    assert (mt != _movie2seq.end());

    uint64_t m = _movie2seq[mid];
    uint64_t n = _user2seq[uid];
    
    if (input_rating_class(rating) > 0) {
      if (!cmap) {
	_nratings++;
	RatingMap *rm = _users2rating[n];
	if (_env.binary_data)
	  (*rm)[m] = 1;
	else {
	  assert (rating > 0);
	  (*rm)[m] = rating;
	}
	_users[n]->push_back(m);
	_movies[m]->push_back(n);
      } else {
	debug("adding test or validation entry for user %d, item %d", n, m);
	Rating r(n,m);
	assert(cmap);
	if (_env.binary_data)
	  (*cmap)[r] = 1;
	else
	  (*cmap)[r] = rating;
      }
    }
    if (_nratings % 1000 == 0) {
      printf("\r+ read %d users, %d movies, %d ratings", 
	     _curr_user_seq, _curr_movie_seq, _nratings);
      fflush(stdout);
    }
  }
  return 0;
}


int
Ratings::write_marginal_distributions()
{
  FILE *f = fopen(Env::file_str("/byusers.tsv").c_str(), "w");
  uint32_t x = 0;
  uint32_t nusers = 0;
  for (uint32_t n = 0; n < _env.n; ++n) {
    const vector<uint32_t> *movies = get_movies(n);
    IDMap::const_iterator it = seq2user().find(n);
    if (!movies || movies->size() == 0) {
      debug("0 movies for user %d (%d)", n, it->second);
      x++;
      continue;
    }
    uint32_t t = 0;
    for (uint32_t m = 0; m < movies->size(); m++) {
      uint32_t mov = (*movies)[m];
      yval_t y = r(n,mov);
      t += y;
    }
    x = 0;
    fprintf(f, "%d\t%d\t%d\t%d\n", n, it->second, movies->size(), t);
    nusers++;
  }
  fclose(f);
  //_env.n = nusers;
  lerr("longest sequence of users with no movies: %d", x);

  f = fopen(Env::file_str("/byitems.tsv").c_str(), "w");
  x = 0;
  uint32_t nitems = 0;
  for (uint32_t n = 0; n < _env.m; ++n) {
    const vector<uint32_t> *users = get_users(n);
    IDMap::const_iterator it = seq2movie().find(n);
    if (!users || users->size() == 0) {
      lerr("0 users for movie %d (%d)", n, it->second);
      x++;
      continue;
    }
    uint32_t t = 0;
    for (uint32_t m = 0; m < users->size(); m++) {
      uint32_t u = (*users)[m];
      yval_t y = r(u,n);
      t += y;
    }
    x = 0;
    fprintf(f, "%d\t%d\t%d\t%d\n", n, it->second, users->size(), t);
    nitems++;
  }
  fclose(f);
  //_env.m = nitems;
  lerr("longest sequence of items with no users: %d", x);
  Env::plog("post pruning nusers:", _env.n);
  Env::plog("post pruning nitems:", _env.m);
  return 0;
}

int
Ratings::read_test_users(FILE *f, UserMap *bmap)
{
  assert (bmap);
  uint32_t uid = 0;
  while (!feof(f)) {
    if (fscanf(f, "%u\n", &uid) < 0) {
      printf("error: unexpected lines in file\n");
      exit(-1);
    }

    IDMap::iterator it = _user2seq.find(uid);
    if (it == _user2seq.end())
      continue;
    uint32_t n = _user2seq[uid];
    (*bmap)[n] = true;
  }
  Env::plog("read %d test users", bmap->size());
  return 0;
}

int
Ratings::read_echonest(string dir)
{
  printf("reading echo nest dataset...\n");
  fflush(stdout);
  uint32_t mcurr = 1, scurr = 1;
  char buf[1024];
  sprintf(buf, "%s/train_triplets.txt", dir.c_str());

  FILE *f = fopen(buf, "r");
  if (!f) {
    fprintf(stderr, "error: cannot open file %s:%s", buf, strerror(errno));
    fclose(f);
    exit(-1);
  }
  uint32_t mid = 0, uid = 0, rating = 0;
  char mids[512], uids[512];
  char b[128];
  while (!feof(f)) {
    if (fscanf(f, "%s\t%s\t%u\n", uids, mids, &rating) < 0) {
      printf("error: unexpected lines in file\n");
      fclose(f);
      exit(-1);
    }

    StrMap::iterator uiditr = _str2id.find(uids);
    if (uiditr == _str2id.end()) {
      _str2id[uids] = scurr;
      scurr++;
    }
    uid = _str2id[uids];
    
    StrMap::iterator miditr = _str2id.find(mids);
    if (miditr == _str2id.end()) {
      _str2id[mids] = mcurr;
      mcurr++;
    }
    mid = _str2id[mids];

    IDMap::iterator it = _user2seq.find(uid);
    if (it == _user2seq.end() && !add_user(uid)) {
      printf("error: exceeded user limit %d, %d, %d\n",
	     uid, mid, rating);
      fflush(stdout);
      continue;
    }
    
    IDMap::iterator mt = _movie2seq.find(mid);
    if (mt == _movie2seq.end() && !add_movie(mid)) {
      printf("error: exceeded movie limit %d, %d, %d\n",
	     uid, mid, rating);
      fflush(stdout);
      continue;
    }
    
    uint32_t m = _movie2seq[mid];
    uint32_t n = _user2seq[uid];

    _user2str[n] = uids;
    _movie2str[m] = mids;

    if (rating > 0) {
      _nratings++;
      RatingMap *rm = _users2rating[n];
      (*rm)[m] = rating;
      _users[n]->push_back(m);
      _movies[m]->push_back(n);
      _ratings.push_back(Rating(n,m));
    }
    if (_nratings % 1000 == 0) {
      printf("\r+ read %d users, %d movies, %d ratings", 
	     _curr_user_seq, _curr_movie_seq, _nratings);
      fflush(stdout);
    }
  }
  fclose(f);
  return 0;
}

int
Ratings::read_nyt(string dir)
{
  printf("reading nyt dataset...\n");
  fflush(stdout);
  uint32_t mcurr = 1, scurr = 1;
  char buf[1024];
  sprintf(buf, "%s/nyt-clicks.tsv", dir.c_str());

  FILE *f = fopen(buf, "r");
  if (!f) {
    fprintf(stderr, "error: cannot open file %s:%s", buf, strerror(errno));
    fclose(f);
    exit(-1);
  }
  uint32_t mid = 0, uid = 0, rating = 0;
  char mids[512], uids[512];
  char b[128];
  while (!feof(f)) {
    if (fscanf(f, "%s\t%s\t%u\n", uids, mids, &rating) < 0) {
      printf("error: unexpected lines in file\n");
      fclose(f);
      exit(-1);
    }

    StrMap::iterator uiditr = _str2id.find(uids);
    if (uiditr == _str2id.end()) {
      _str2id[uids] = scurr;
      scurr++;
    }
    uid = _str2id[uids];
    
    StrMap::iterator miditr = _str2id.find(mids);
    if (miditr == _str2id.end()) {
      _str2id[mids] = mcurr;
      mcurr++;
    }
    mid = _str2id[mids];

    IDMap::iterator it = _user2seq.find(uid);
    if (it == _user2seq.end() && !add_user(uid)) {
      printf("error: exceeded user limit %d, %d, %d\n",
	     uid, mid, rating);
      fflush(stdout);
      continue;
    }
    
    IDMap::iterator mt = _movie2seq.find(mid);
    if (mt == _movie2seq.end() && !add_movie(mid)) {
      printf("error: exceeded movie limit %d, %d, %d\n",
	     uid, mid, rating);
      fflush(stdout);
      continue;
    }
    
    uint32_t m = _movie2seq[mid];
    uint32_t n = _user2seq[uid];

    _user2str[n] = uids;
    _movie2str[m] = mids;

    if (rating > 0) {
      _nratings++;
      RatingMap *rm = _users2rating[n];
      (*rm)[m] = rating;
      _users[n]->push_back(m);
      _movies[m]->push_back(n);
      _ratings.push_back(Rating(n,m));
    }
    if (_nratings % 1000 == 0) {
      printf("\r+ read %d users, %d movies, %d ratings", 
	     _curr_user_seq, _curr_movie_seq, _nratings);
      fflush(stdout);
    }
  }
  fclose(f);

  sprintf(buf, "%s/str2id.tsv", dir.c_str());

  uint32_t q = 0;
  f = fopen(buf, "w");
  for (StrMap::const_iterator i = _str2id.begin(); i != _str2id.end(); ++i) {
    if (strlen(i->first.c_str()) >= strlen("10219231518"))  {
      fprintf(f, "%s\t%d\n", i->first.c_str(), i->second);
      q++;
    }
  }
  fclose(f);
  Env::plog("wrote %d str2id entries", q);
  return 0;
}

int
Ratings::read_mendeley(string dir)
{
  char buf[1024];
  sprintf(buf, "%s/users.dat", dir.c_str());
  
  info("reading from %s\n", buf);

  FILE *f = fopen(buf, "r");
  if (!f) {
    fprintf(stderr, "error: cannot open file %s:%s", buf, strerror(errno));
    fclose(f);
    exit(-1);
  }
  
  uint32_t uid = 1, rating = 0;
  char b[128];
  while (!feof(f)) {
    vector<uint32_t> mids;
    uint32_t len = 0;
    if (fscanf(f, "%u\t", &len) < 0) {
      printf("error: unexpected lines in file\n");
      fclose(f);
      exit(-1);
    }

    uint32_t mid = 0;
    for (uint32_t i = 0; i < len; ++i) {
      if (i == len - 1) {
	if (fscanf(f, "%u\t", &mid) < 0) {
	  printf("error: unexpected lines in file\n");
	  fclose(f);
	  exit(-1);
	}
	mids.push_back(mid);
      } else {
	if (fscanf(f, "%u", &mid) < 0) {
	  printf("error: unexpected lines in file\n");
	  fclose(f);
	  exit(-1);
	}
	mids.push_back(mid);
      }
    }
    
    IDMap::iterator it = _user2seq.find(uid);
    if (it == _user2seq.end() && !add_user(uid)) {
      printf("error: exceeded user limit %d, %d, %d\n",
	     uid, mid, rating);
      fflush(stdout);
      continue;
    }
    
    for (uint32_t idx = 0; idx < mids.size(); ++idx) {
      uint32_t mid = mids[idx];
      IDMap::iterator mt = _movie2seq.find(mid);
      if (mt == _movie2seq.end() && !add_movie(mid)) {
	printf("error: exceeded movie limit %d, %d, %d\n",
	       uid, mid, rating);
	fflush(stdout);
	continue;
      }
      uint32_t m = _movie2seq[mid];
      uint32_t n = _user2seq[uid];

      yval_t rating = 1.0;
      _nratings++;
      RatingMap *rm = _users2rating[n];
      (*rm)[m] = rating;
      _users[n]->push_back(m);
      _movies[m]->push_back(n);
      _ratings.push_back(Rating(n,m));
    }
    uid++;
    if (_nratings % 1000 == 0) {
      printf("\r+ read %d users, %d movies, %d ratings", 
	     _curr_user_seq, _curr_movie_seq, _nratings);
      fflush(stdout);
    }
  }
  fclose(f);
  return 0;
}

int
Ratings::read_netflix_movie(string dir, uint32_t movie)
{
  char buf[1024];
  sprintf(buf, "%s/mv_%.7d.txt", dir.c_str(), movie);

  info("reading from %s\n", buf);

  FILE *f = fopen(buf, "r");
  if (!f) {
    fprintf(stderr, "error: cannot open file %s:%s", buf, strerror(errno));
    fclose(f);
    exit(-1);
  }

  uint32_t mid = 0;
  if (!fscanf(f, "%d:\n", &mid)) {
    fclose(f);
    return -1;
  }
  assert (mid == movie);
  
  IDMap::iterator mt = _movie2seq.find(mid);
  if (mt == _movie2seq.end() && !add_movie(mid)) {
    fclose(f);
    return 0;
  }
  
  uint32_t m = _movie2seq[mid];
  uint32_t uid = 0, rating = 0;
  char b[128];
  while (!feof(f)) {
    if (fscanf(f, "%u,%u,%*s\n", &uid, &rating, b) < 0) {
	printf("error: unexpected lines in file\n");
	fclose(f);
	exit(-1);
    }

    IDMap::iterator it = _user2seq.find(uid);
    if (it == _user2seq.end() && !add_user(uid))
      continue;

    uint32_t n = _user2seq[uid];

    if (rating > 0) {
      _nratings++;
      RatingMap *rm = _users2rating[n];
      (*rm)[m] = rating;
      _users[n]->push_back(m);
      _movies[m]->push_back(n);
      _ratings.push_back(Rating(n,m));
    }
  }
  fclose(f);
  printf("\r+ read %d users, %d movies, %d ratings", 
	 _curr_user_seq, _curr_movie_seq, _nratings);
  fflush(stdout);
  return 0;
}

// Method that reads the data from the movielens database
int
Ratings::read_movielens(string dir)
{
  // Store name of train data file in buf
  
  char buf[1024];
  sprintf(buf, "%s/ml-1m_train.tsv", dir.c_str());

  info("reading from %s\n", buf);

  FILE *f = fopen(buf, "r");
  if (!f) {
    fprintf(stderr, "error: cannot open file %s:%s", buf, strerror(errno));
    fclose(f);
    exit(-1);
  }
  
  uint32_t mid = 0, uid = 0, rating = 0;
  char b[128];
  while (!feof(f)) {
    if (fscanf(f, "%u\t%u\t%u\n", &uid, &mid, &rating) < 0) {
      printf("error: unexpected lines in file\n");
      fclose(f);
      exit(-1);
    }

    IDMap::iterator it = _user2seq.find(uid);
    if (it == _user2seq.end() && !add_user(uid)) {
      printf("error: exceeded user limit %d, %d, %d\n",
	     uid, mid, rating);
      fflush(stdout);
      continue;
    }
    
    IDMap::iterator mt = _movie2seq.find(mid);
    if (mt == _movie2seq.end() && !add_movie(mid)) {
      printf("error: exceeded movie limit %d, %d, %d\n",
	     uid, mid, rating);
      fflush(stdout);
      continue;
    }
    
    uint32_t m = _movie2seq[mid];
    uint32_t n = _user2seq[uid];

    if (rating > 0) {
      _nratings++;
      RatingMap *rm = _users2rating[n];
      (*rm)[m] = rating;
      _users[n]->push_back(m);
      _movies[m]->push_back(n);
      _ratings.push_back(Rating(n,m));
    }
  }
  fclose(f);
  return 0;
}


void
Ratings::load_movies_metadata(string s)
{
  if (_env.dataset == Env::MOVIELENS)
    read_movielens_metadata(s);
  else if (_env.dataset == Env::NETFLIX)
    read_netflix_metadata(s);
  else if (_env.dataset == Env::MENDELEY)
    read_mendeley_metadata(s);
}

int
Ratings::read_movielens_metadata(string dir)
{
  uint32_t n = 0;
  char buf[1024];
  sprintf(buf, "movies.tsv", dir.c_str());
  FILE *f = fopen(buf, "r");
  assert(f);
  uint32_t id;
  char name[4096];
  char type[4096];
  char *line = (char *)malloc(4096);
  while (!feof(f)) {
    if (fgets(line, 4096, f) == NULL)
      break;
    uint32_t k = 0;
    char *p = line;
    const char r[3] = "#";
    do {
      char *q = NULL;
      char *d = strtok_r(p, r, &q);
      if (q == p)
	break;
      if (k == 0) {
	id = atoi(d);
	id = _movie2seq[id];
      } else if (k == 1) {
	strcpy(name, d);
	_movie_names[id] = name;
	debug("%d -> %s", id, name);
      } else if (k == 2) {
	strcpy(type, d);
	_movie_types[id] = type;
	debug("%d -> %s", id, type);
      }
      p = q;
      k++;
    } while (p != NULL);
    n++;
    debug("read %d lines\n", n);
    memset(line, 0, 4096);
  }
  free(line);
  return 0;
}

int
Ratings::read_netflix_metadata(string dir)
{
  uint32_t n = 0;
  char buf[1024];
  sprintf(buf, "movie_titles.txt", dir.c_str());
  FILE *f = fopen(buf, "r");
  assert(f);
  uint32_t id, year;
  char name[4096];
  char *line = (char *)malloc(4096);
  while (!feof(f)) {
    if (fgets(line, 4096, f) == NULL)
      break;
    uint32_t k = 0;
    char *p = line;
    const char r[3] = ",";
    do {
      char *q = NULL;
      char *d = strtok_r(p, r, &q);
      if (q == p)
	break;
      if (k == 0) {
	id = atoi(d);
	lerr("%d: ", id);
	id = _movie2seq[id];
	lerr("%d -> ", id);
      } else if (k == 1) {
	year  = atoi(d); // skip
	lerr("%d", year);
      } else if (k == 2) {
	strcpy(name, d);
	_movie_names[id] = name;
	_movie_types[id] = "";
	lerr("%d -> %s", id, name);
      }
      p = q;
      k++;
    } while (p != NULL);
    n++;
    debug("read %d lines\n", n);
    memset(line, 0, 4096);
  }
  free(line);
  return 0;
}

int
Ratings::read_mendeley_metadata(string dir)
{
  uint32_t n = 0;
  char buf[1024];
  sprintf(buf, "%s/titles.dat", dir.c_str());
  FILE *f = fopen(buf, "r");
  assert(f);
  char name[4096];
  char *line = (char *)malloc(4096);
  uint32_t id = 0;
  while (!feof(f)) {
    if (fgets(line, 4096, f) == NULL)
      break;
    strcpy(name, line);
    uint32_t seq = _movie2seq[id];
    _movie_names[seq] = name;
    id++;
  }
  lerr("read %d lines\n", n);
  free(line);
    return 0;
}


string
Ratings::movies_by_user_s() const
{
  ostringstream sa;
  sa << "\n[\n";
  for (uint32_t i = 0; i < _users.size(); ++i) {
    IDMap::const_iterator it = _seq2user.find(i);
    sa << it->second << ":";
    vector<uint32_t> *v = _users[i];
    if (v)  {
      for (uint32_t j = 0; j < v->size(); ++j) {
	uint32_t m = v->at(j);
	IDMap::const_iterator mt = _seq2movie.find(m);
	sa << mt->second;
	if (j < v->size() - 1)
	  sa << ", ";
      }
      sa << "\n";
    }
  }
  sa << "]";
  return sa.str();
}

FreqMap
Ratings::validation_users_of_movie() {
  return _validation_users_of_movie;
}

IDMap
Ratings::leave_one_out() {
  return _leave_one_out;
}