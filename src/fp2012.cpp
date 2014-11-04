/* author : Shane Neph
   date:    2008
*/

#include <algorithm>
#include <cstddef>
#include <deque>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>


namespace {

  struct Bed; // forward decl

  // typedefs
  typedef double T;
  typedef std::deque<T> AllInRange;
  typedef std::deque<Bed> BedList;
  typedef std::pair< double, std::pair<BedList::const_iterator, double> > PType;


  // Local globals; defaults overridden by user
  // flankmin and flankmax must each be greater than 1
  std::size_t flankmin = 1;
  std::size_t flankmax = 2;
  std::size_t cmin = 1;
  std::size_t cmax = 2;
  const int maxrecurses = 20; // applies to each of left/right in reanalyze
  const double BADSCORE = 10000;


  //=====
  // Bed
  //=====
  struct Bed {
    Bed() : score_(BADSCORE), leftmean_(-1000), centermean_(1000), rightmean_(-1000),
            leftmost_(-1), rightmost_(-1), innerright_(-10000),
            pos_(std::numeric_limits<unsigned long>::max()) { /* */ }

    friend std::ostream& operator<<(std::ostream& os, const Bed& b) {
      // leftmost_ & rightmost_ & innerright_ are absolute counts of bases
      //  doesn't make sense to think of them as 0-based
      // must deal with that here
      os << (1 + b.pos_ - b.leftmost_)
         << "\t"
         << (b.pos_ + 1)
         << "\t"
         << (b.pos_ + b.innerright_ )
         << "\t"
         << (b.pos_ + b.rightmost_)
         << "\t"
         << b.score_
         << "\t"
         << b.leftmean_
         << "\t"
         << b.centermean_
         << "\t"
         << b.rightmean_;
      return(os);
    }

    double score_;
    double leftmean_, centermean_, rightmean_;
    long leftmost_, rightmost_, innerright_;
    unsigned long pos_;
  };

  //=========
  // usage()
  //=========
  std::string usage() {
    std::string msg =  "Usage: \n";
                msg += "\t[--help]\n";
                msg += "\t[--flankmin <bases>  = 6]\n";
                msg += "\t[--flankmax <bases>  = 12]\n";
                msg += "\t[--centermin <bases> = 6]\n";
                msg += "\t[--centermax <bases> = 100]\n";
                msg += "\t[--maxthold <value>  = 10]\n";
                msg += "\tFile-Full-O-Integers\n";
    return(msg);
  }

  struct Help {};

  //=======
  // Input
  //=======
  struct Input {
    Input(int argc, char** argv) {
      if ( argc < 2 )
        throw(false);

      vals_.insert(std::make_pair("--flankmin", 6));
      vals_.insert(std::make_pair("--flankmax", 12));
      vals_.insert(std::make_pair("--centermin", 6));
      vals_.insert(std::make_pair("--centermax", 100));
      vals_.insert(std::make_pair("--maxthold", 10));

      for ( int i = 1; i < argc - 1; i += 2 ) {
        if ( argv[i] == std::string("--help") )
          throw(Help());
        else if ( argv[i+1] == std::string("--help") )
          throw(Help());
        else if ( vals_.find(argv[i]) == vals_.end() ) {
          std::cerr << "Unknown Argument: " << argv[i] << std::endl;
          throw (false);
        }
        std::stringstream con(argv[i + 1]);
        int tmp;
        con >> tmp;
        vals_[argv[i]] = tmp;
      } // for
      file_ = argv[argc-1];
      if ( file_ == "--help" )
        throw(Help());
      else if ( !check() ) {
        std::cerr << "Bad numeric argument value(s) detected" << std::endl;
        throw(false);
      }
    }

    int CenterMax() const
      { return(static_cast<int>(vals_.find("--centermax")->second)); }
    int CenterMin() const
      { return(static_cast<int>(vals_.find("--centermin")->second)); }
    int FlankMax() const
      { return(static_cast<int>(vals_.find("--flankmax")->second)); }
    int FlankMin() const
      { return(static_cast<int>(vals_.find("--flankmin")->second)); }
    double MaxThold() const
      { return(static_cast<int>(vals_.find("--maxthold")->second)); }
    std::string File() const { return(file_); }

  private:
    bool check() {
      if ( CenterMin() <= 1 ) // make sure at least 2
        return(false);
      else if ( FlankMax() < FlankMin() )
        return(false);
      else if ( FlankMin() <= 1 ) // make sure at least 2
        return(false);
      else if ( CenterMax() < CenterMin() )
        return(false);
      else if ( MaxThold() <= 0 )
        return(false);
      return(true);
    }

  private:
    std::map<std::string, double> vals_;
    std::string file_;
  };

  //===================
  // overlap_centers()
  //===================
  inline bool overlap_centers(const Bed& db1, const Bed& db2) {
    if ( db1.pos_ < db2.pos_ ) // +1 applies to both -> safe to ignore
      return(db1.pos_ + db1.innerright_ > db2.pos_);
    return(db2.pos_ + db2.innerright_ > db1.pos_);
  }
/*
  //=================
  // operator< : Bed
  //=================
  inline bool operator<(const Bed& db1, const Bed& db2) {
    // center overlaps are treated as equality; though this isn't the same as operator== below
    if ( overlap_centers(db1, db2) )
      return(false);
    return(db1.pos_ < db2.pos_);
  }
*/
  //=================
  // operator== : Bed
  //=================
  inline bool operator==(const Bed& db1, const Bed& db2) {
    return(
      db1.score_      == db2.score_      &&
      db1.leftmean_   == db2.leftmean_   &&
      db1.centermean_ == db2.centermean_ &&
      db1.rightmean_  == db2.rightmean_  &&
      db1.leftmost_   == db2.leftmost_   &&
      db1.rightmost_  == db2.rightmost_  &&
      db1.innerright_ == db2.innerright_ &&
      db1.pos_        == db2.pos_
          );
  }

  //============
  // max_mean()
  //============
  template <typename RAIter>
  std::pair<double, RAIter> max_mean(RAIter mv, RAIter end, int minmoves) {
    // surprisingly (to me) std::numeric_limits<double>::min() is a positive value!!  What the...?
    double mxmn = static_cast<double>(std::numeric_limits<long>::min());
    RAIter endmark = mv;
    double sum = 0, cntr = 0;
    while ( mv != end ) {
      sum += *mv;
      if ( ++cntr >= minmoves ) {
        double tmp = sum / cntr;
        if ( tmp >= mxmn ) {
          mxmn = tmp;
          endmark = mv + 1;
        }
      }
      ++mv;
    } // while
    return(std::make_pair(mxmn, endmark));
  }

  //=============
  // score_first
  //=============
  PType score_first(const BedList& bl, AllInRange::const_iterator beg) {
    // only scoring first element of 'bl', which is our L representative
    BedList::const_iterator rightiter = bl.begin() + 1;
    double mnstat = BADSCORE, stat;
    double L = bl.front().leftmean_, C = BADSCORE, R;
    double bestC = C;
    if ( L > 0 && bl.size() > cmin ) {
      for ( BedList::const_iterator right = bl.begin() + 1 + cmin; right != bl.end(); ++right ) {
        if ( bl.size() >= 1 + cmax && right > bl.begin() + 1 + cmax )
          break;
        else if ( (R = right->rightmean_) > 0 ) {
          std::size_t length = right - bl.begin() - 1;
          C = 1.0 + std::accumulate(beg + 1, beg + 1 + length, 0.0) / length;
          if ( L > C && R > C ) {
            stat = C/L + C/R; // 1 was already added to C
            if ( stat <= mnstat ) {
              mnstat = stat;
              rightiter = right;
              bestC = C - 1.0;
            }
          }
        }
      } // for
    }
    return(std::make_pair(mnstat, std::make_pair(rightiter, bestC)));
  }

  //=============
  // reanalyze() : recursively look for footprints that have been trumped, but are still real
  //=============
  BedList reanalyze(BedList& cache, AllInRange& cacheRef, const Bed& currCacheFootprint, const Bed& currFootprint, double thold, int maxrec, int rec) {
    if ( rec >= maxrec )
      return(BedList());

    // For every item of cache, optimal L has already been determined
    // Want to add constraints of currFootprint being removed and calculate items on each side of partition
    // Assume that currCacheFootprint is an element of cache (which it better be)
    BedList toRtn;
    BedList left, right;
    std::size_t rightRefCntr = 0;
    for ( BedList::iterator i = cache.begin(); i != cache.end(); ++i ) {
      if ( *i == currCacheFootprint ) {
        std::copy(cache.begin(), i+1, std::back_inserter(left)); // currCacheFootprint can serve as R element and keep disjoint/non-adjacent property
        // blank out central footprint region from 'right'
        int counts = 0;
        int nummoves = currFootprint.innerright_;
        while ( i != cache.end() && counts < nummoves )
          ++i, ++counts, ++rightRefCntr;
        std::copy(i, cache.end(), std::back_inserter(right));
        break;
      }
      ++rightRefCntr;
    } // for


    // Left side
    { // local scope for memory
      AllInRange leftRefCache;
      Bed bestLeftFull, bestLeft;
      bestLeftFull.score_ = thold + 1;
      bestLeft.score_ = thold + 1;
      std::size_t leftRefCntr = 0;
      BedList tmpLeft(left.begin(), left.end());
      while ( tmpLeft.size() ) {
        BedList::iterator i = tmpLeft.begin();
        PType pt = score_first(tmpLeft, cacheRef.begin() + leftRefCntr);
        if ( pt.first < thold && pt.first < bestLeftFull.score_ ) {
          if ( leftRefCache.empty() )
            std::copy(cacheRef.begin() + leftRefCntr, cacheRef.end(), std::back_inserter(leftRefCache));
          Bed candidate = *i;
          candidate.score_ = pt.first;
          candidate.rightmost_ = (pt.second.first - i) + pt.second.first->rightmost_;
          candidate.rightmean_ = pt.second.first->rightmean_;
          candidate.innerright_ = (pt.second.first - i);
          candidate.centermean_ = pt.second.second;
          bestLeftFull = candidate;
          bestLeft = *i;
        }
        tmpLeft.pop_front();
        ++leftRefCntr;
      } // for

      // Curse again
      if ( bestLeftFull.score_ < thold ) {
        BedList tmp = reanalyze(left, leftRefCache, bestLeft, bestLeftFull, thold, maxrec, rec+1);
        std::copy(tmp.begin(), tmp.end(), std::back_inserter(toRtn));
      }
      left.clear();
    } // local scope

    // Add currFootprint now to keep genomic ordering proper
    toRtn.push_back(currFootprint);


    // Right side
    { // local scope for memory
      AllInRange rightRefCache;
      Bed bestRightFull, bestRight;
      bestRightFull.score_ = thold + 1;
      bestRight.score_ = thold + 1;
      BedList tmpRight(right.begin(), right.end());
      while ( tmpRight.size() ) {
        BedList::iterator i = tmpRight.begin();
        PType pt = score_first(tmpRight, cacheRef.begin() + rightRefCntr);
        if ( pt.first < thold && pt.first < bestRightFull.score_ ) {
          if ( rightRefCache.empty() )
            std::copy(cacheRef.begin() + rightRefCntr, cacheRef.end(), std::back_inserter(rightRefCache));
          Bed candidate = *i;
          candidate.score_ = pt.first;
          candidate.rightmost_ = (pt.second.first - i) + pt.second.first->rightmost_;
          candidate.rightmean_ = pt.second.first->rightmean_;
          candidate.innerright_ = (pt.second.first - i);
          candidate.centermean_ = pt.second.second;
          bestRightFull = candidate;
          bestRight = *i;
        }
        tmpRight.pop_front();
        ++rightRefCntr;
      } // for

      // Curse again
      if ( bestRightFull.score_ < thold ) {
        BedList tmp = reanalyze(right, rightRefCache, bestRight, bestRightFull, thold, maxrec, rec+1);
        std::copy(tmp.begin(), tmp.end(), std::back_inserter(toRtn));
      }
      right.clear();
    } // local scope


    // Clean up and return
    cache.clear();
    cacheRef.clear();
    return(toRtn);
  }

  //========
  // read()
  //========
  void read(std::istream_iterator<T> start, std::istream_iterator<T> end, double thold) {
    // flankmin and flankmax must each be greater than one
    AllInRange ref, cacheRef;
    BedList current, cache;

    const std::size_t min2score = flankmin + cmax + flankmax; // element flankmin-1 (0-based)
    std::size_t totalcount = 0;
    while ( start != end ) {
      ref.push_back(*start++);
      if ( ++totalcount == flankmin + flankmax - 1 ) // |L|min + |R|max
        break;
    } // while

    if ( start == end )
      return;

    AllInRange::const_reverse_iterator lstart = ref.rbegin() + flankmax - 1;
    AllInRange::const_reverse_iterator lstop = lstart + std::min(flankmax, static_cast<std::size_t>(ref.rend() - lstart));
    std::pair<double, AllInRange::const_reverse_iterator> lrtn = max_mean(lstart, lstop, flankmin);

    AllInRange::const_iterator rstart = ref.begin() + flankmin - 1; // points to same actual base as lstart
    AllInRange::const_iterator rstop = rstart + std::min(flankmax, static_cast<std::size_t>(ref.end() - rstart));
    std::pair<double, AllInRange::const_iterator> rrtn = max_mean(rstart, rstop, flankmin);

    Bed tmp;
    tmp.leftmean_ = lrtn.first;
    tmp.rightmean_ = rrtn.first;
    tmp.leftmost_ = lrtn.second - lstart;
    tmp.rightmost_ = rrtn.second - rstart;
    tmp.pos_ = flankmin - 1;
    current.push_back(tmp);

    Bed currFootprint, currCacheFootprint; // construction creates uninteresting obj
    bool first = true;

    std::size_t poscount = flankmin; // one more than current 'pos_' of 'tmp'
    while ( start != end ) {
      ref.push_back(*start++);

      //------------
      // Left Side
      //  All of my speedup attempts have failed in similar ways as the example shown
      //   below for the right side.
      lstart = ref.rbegin() + flankmax - 1;
      lstop = lstart + std::min(flankmax, static_cast<std::size_t>(ref.rend() - lstart));
      lrtn = max_mean(lstart, lstop, flankmin);
      tmp.leftmean_ = lrtn.first;
      tmp.leftmost_ = lrtn.second - lstart;

      //------------
      // Right Side
      //  Cannot seem to speed up robustly, consider:
      //    initial = { 5, 6, 7, 8, 9, 10, 10, 9, 8, 1 } ; mean = 8, where 1 is beyond end.
      //    move = { 6, 7, 8, 9, 10, 10, 9, 8, 1 } ; mean = 8.428, where 8 is beyond end.
      //  Notice that after dropping the 5, the 8 no longer belongs to the optimal mean.
      //  This is not a simple consequence of allowing ties to extend the optimal R-value.
      rstart = ref.begin() + (ref.size() - 1 - (lstart - ref.rbegin())); // point to same base as lstart
      rstop = rstart + std::min(flankmax, static_cast<std::size_t>(ref.end() - rstart));
      rrtn = max_mean(rstart, rstop, flankmin);
      tmp.rightmean_ = rrtn.first;
      tmp.rightmost_ = rrtn.second - rstart;
      tmp.pos_ = poscount++;

      current.push_back(tmp);

      if ( ++totalcount >= min2score ) { // we have enough elements to derive a score
        // give the element to score:
        //  o a statistic score
        //  o a new rightmost_ postion == full footprint's rightmost base
        //  o a innerright_ position == where the right flanking bases start
        //  o a centermean_ value == lower part of the footprint
        //  o reassign rightmean_ value
        // there are more elements in 'ref' than in current -> adjust
        PType pt = score_first(current, ref.begin() + flankmin - 1);
        if ( !cache.empty() ) {
          cache.push_back(current.front());
          cacheRef.push_back(ref.back());
        }

        if ( pt.first < thold ) {
          if ( cache.empty() ) {
            cache.push_back(current.front());
            std::copy(ref.begin() + flankmin - 1, ref.end(), std::back_inserter(cacheRef));
          }
          Bed candidate = current.front();
          candidate.score_ = pt.first;
          candidate.rightmost_ = (pt.second.first - current.begin()) + pt.second.first->rightmost_;
          candidate.rightmean_ = pt.second.first->rightmean_;
          candidate.innerright_ = (pt.second.first - current.begin());
          candidate.centermean_ = pt.second.second;

          if ( first || overlap_centers(currFootprint, candidate) ) {
            if ( currFootprint.score_ > candidate.score_ ) // new score is better; must not have >= here or miss correct central fp calls on LHS
              currFootprint = candidate, currCacheFootprint = current.front();
          } else {
            if ( !first ) { // spit results
              BedList entries = reanalyze(cache, cacheRef, currCacheFootprint, currFootprint, thold, maxrecurses, 1);
              for ( std::size_t i = 0; i < entries.size(); ++i )
                std::cout << entries[i] << std::endl;
            }
            currFootprint = candidate;
          }
          first = false;
        }

        current.pop_front();
        ref.pop_front();
      }
    } // while

    if ( !first ) { // still more to spit
      BedList entries = reanalyze(cache, cacheRef, currCacheFootprint, currFootprint, thold, maxrecurses, 1);
      for ( std::size_t i = 0; i < entries.size(); ++i )
        std::cout << entries[i] << std::endl;
    }
  }

} // unnamed



//========
// main()
//========
int main(int argc, char** argv) {
  try {
    // check inputs
    Input input(argc, argv);

    // get user arguments
    std::ifstream infile(input.File().c_str());
    if ( !infile && input.File() != "-" ) {
      std::cerr << "Unable to find input file: " << input.File() << std::endl;
      return(EXIT_FAILURE);
    }
    std::istream_iterator<T> ins(input.File() == "-" ? std::cin : infile), eos;

    flankmin = input.FlankMin();
    flankmax = input.FlankMax();
    cmin = input.CenterMin();
    cmax = input.CenterMax();
    double maxthold = input.MaxThold();

    // do work
    read(ins, eos, maxthold);

    return(EXIT_SUCCESS);
  } catch(Help h) {
    std::cout << usage() << std::endl;
    return(EXIT_SUCCESS);
  } catch(bool b) {
    std::cerr << usage() << std::endl;
  } catch(std::exception& s) {
    std::cerr << "Standard Exception: " << s.what() << std::endl;
  } catch(...) {
    std::cerr << "Unknown Exception" << std::endl;
  }
  return(EXIT_FAILURE);
}
