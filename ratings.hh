#ifndef RATINGS_HH
#define RATINGS_HH

#include <string>
#include <vector>
#include <queue>
#include <map>
#include <stdint.h>
#include "matrix.hh"
#include "env.hh"
//#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf.h>

using namespace std;

typedef std::map<Rating, D1Array<uint64_t>> AvailabilityMap;

class Ratings {
public:
  Ratings(Env &env, uint64_t* (*fptr) (uint64_t, uint64_t, uint32_t &)):
    _users2rating(env.n),
    _users(env.n),
    _movies(env.m),
    _userObs(env.n,env.uc,true),
    _itemObs(env.m,env.ic,true),
    _userObsScale(env.uc),
    _itemObsScale(env.ic),
    _env(env),
    _curr_user_seq(0), 
    _curr_movie_seq(0),
    _nratings(0),
    _likes(0),
    _offset(env.offset){
	getAvailableItems = fptr;
    }
  ~Ratings() { }

  int read(string s);
  void readObserved(string s);
  void readValidationAndTest(string s);
  uint32_t input_rating_class(uint32_t v) const;
  bool test_hit(uint32_t v) const;
  int write_marginal_distributions();
  
  const SparseMatrix &users() const { return _users; }
  SparseMatrix &users() { return _users; }

  const SparseMatrix &movies() const { return _movies; }
  SparseMatrix &movies() { return _movies; }
  
  uint32_t n() const;
  uint32_t m() const;
  uint32_t r(uint32_t i, uint32_t j) const;
  uint32_t nratings() const { return _nratings; }
  uint32_t likes() const { return _likes; }
  const vector<Rating> &allratings() const { return _ratings; }
 
  const vector<uint32_t> *get_users(uint32_t a);
  const vector<uint32_t> *get_movies(uint32_t a);

  const IDMap &user2seq() const { return _user2seq; }
  const IDMap &seq2user() const { return _seq2user; }

  const IDMap &movie2seq() const { return _movie2seq; }
  const IDMap &seq2movie() const { return _seq2movie; }
  int read_generic(FILE *f, CountMap *m);
  
  int read_nyt_train(FILE *f, CountMap *m);
  int read_nyt_titles(string dir);
  
  void load_movies_metadata(string dir);
  int read_test_users(FILE *f, UserMap *);

  string movie_type(uint32_t movie_seq) const;
  string movie_name(uint32_t movie_seq) const;
  int read_netflix_movie(string dir, uint32_t movie);
  int read_netflix_metadata(string dir);
  int read_movielens_metadata(string dir);
  double getAvailability(uint32_t user, uint32_t item){ //Returns value stored in avblty[(user,item)]
	uint64_t uid = _seq2user[user];
	uint64_t itemid = _seq2movie[item];
	Rating elem(uid, itemid);
	return avblty[elem];
  }

  FreqMap validation_users_of_movie();
  IDMap leave_one_out();
  
  Matrix _userObs;    // User characteristics
  Matrix _itemObs;    // Item characteristics
  
  Array _userObsScale;
  Array _itemObsScale;
  
  CountMap _validation_map;
  CountMap _test_map;
  
  uint32_t totRating;
  
private:
  uint64_t* (*getAvailableItems) (uint64_t uid, uint64_t sid, uint32_t &numItems);
  void read_generic_train(string dir);
  int read_movielens(string dir);
  int read_mendeley(string dir);
  int read_echonest(string dir);
  int read_nyt(string dir);
  int read_mendeley_metadata(string dir);
  string movies_by_user_s() const;
  bool add_movie(uint64_t id);
  bool add_user(uint64_t id);
  
  int _offset;

  SparseMatrixR _users2rating;
  SparseMatrix _users;
  SparseMatrix _movies;
  vector<Rating> _ratings;

  Env &_env;
  IDMap _user2seq;
  IDMap _movie2seq;
  IDMap _seq2user;
  IDMap _seq2movie;
  StrMap _str2id;
  StrMapInv _user2str;
  StrMapInv _movie2str;
  uint32_t _curr_user_seq;
  uint32_t _curr_movie_seq;
  uint32_t _nratings;
  uint32_t _likes;
  StrMapInv _movie_names;
  StrMapInv _movie_types;
  
  FreqMap _validation_users_of_movie;
  IDMap _leave_one_out;
  AvailabilityMap userSess2avblItems;
  ValueMap avblty; // This measure maps Rating class (std::pair of user and item) to total availability over all sessions
};

inline uint32_t
Ratings::n() const
{
  return _users.size();
}

inline uint32_t
Ratings::m() const
{
  return _movies.size();
}

// Adds a new user to the lists
inline bool
Ratings::add_user(uint64_t id)
{
  if (_curr_user_seq >= _env.n) {
    debug("max users %d reached", _env.n);
    cout << "max users reached" << id << endl; 
    return false;
  }
    
  // Adds new values to the lists to convert from user to index
  _user2seq[id] = _curr_user_seq;
  _seq2user[_curr_user_seq] = id;

  assert (!_users[_curr_user_seq]);
  
  // Pointer to the vector of users
  std::vector<uint32_t> **v = _users.data();
  // Adds the new user to the vector
  v[_curr_user_seq] = new vector<uint32_t>;
  
  // Pointer to the map of ratings
  RatingMap **rm = _users2rating.data();
  // Adds a new map of ratings for the new user
  rm[_curr_user_seq] = new RatingMap;
  
  //Icreasees the number of users
  _curr_user_seq++;
  return true;
}

// Adds a new user to the list
inline bool
Ratings::add_movie(uint64_t id)
{
  if (_curr_movie_seq >= _env.m) {
    debug("max movies %d reached", _env.m);
    return false;
  }
  
  // Adds new values to the lists to convert from item to index
  _movie2seq[id] = _curr_movie_seq;
  _seq2movie[_curr_movie_seq] = id;

  assert (!_movies[_curr_movie_seq]);
  
  // Pointer to the vector of users
  std::vector<uint32_t> **v = _movies.data();
  // Adds the new item to the vector
  v[_curr_movie_seq] = new vector<uint32_t>;
  _curr_movie_seq++;
  return true;
}

// Finds the rating for a user and an item given their indices, a for user and b for item
inline uint32_t
Ratings::r(uint32_t a, uint32_t b) const
{
  // Checks that the indices don't exceed the total number of users or items
  assert (a < _env.n && b < _env.m);
  // Pointer to user a's rating
  const RatingMap *rm = _users2rating[a];
  assert(rm);
  
  // Saves *rm as rmc by reference
  const RatingMap &rmc = *rm;
  RatingMap::const_iterator itr = rmc.find(b);
  if (itr == rmc.end())
    // Returns zero if b isn't found in a's RagingMap
    return 0;
  else
    // Returns a's rating for b
    return itr->second;
}

inline const vector<uint32_t> *
Ratings::get_users(uint32_t a)
{
  assert (a < _movies.size());
  const vector<uint32_t> *v = _movies[a];
  return v;
}

inline const vector<uint32_t> *
Ratings::get_movies(uint32_t a)
{
  assert (a < _users.size());
  const vector<uint32_t> *v = _users[a];
  return v;
}

inline bool
Ratings::test_hit(uint32_t v) const
{
  if (_env.binary_data)
    return v >= 1;
  return v >= _env.rating_threshold;
}

// Converts original rating into binary rating
inline uint32_t
Ratings::input_rating_class(uint32_t v) const
{
  if (!_env.binary_data)
    return v;
  return v >= _env.rating_threshold  ? 1 : 0;
}

inline string
Ratings::movie_name(uint32_t movie_seq) const
{
  assert (movie_seq < _env.m);
  StrMapInv::const_iterator i = _movie_names.find(movie_seq);
  if (i != _movie_names.end())
    return i->second;
  return "";
}

inline string
Ratings::movie_type(uint32_t movie_seq) const
{
  assert (movie_seq < _env.m);
  StrMapInv::const_iterator i = _movie_types.find(movie_seq);
  if (i != _movie_types.end())
    return i->second;
  return "";
}

#endif
