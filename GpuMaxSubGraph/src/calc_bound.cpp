#include <vector>
#include <algorithm> // for std::min
#include "main.hpp"


/*
  # maximum possible mapping size = current mapping size +
  # sum of(min(num of nodes with specified label in g, num of nodes with specified label in h)) for each label
  # present in both graphs
*/
int calc_bound( std::vector<LabelClass> label_classes) {
  int bound = 0;
  for ( const LabelClass& label_class : label_classes) {
      if( label_class.g.size() > label_class.h.size() ) bound = bound + label_class.h.size();
      else bound = bound + label_class.g.size();
  }
  return bound;
}
