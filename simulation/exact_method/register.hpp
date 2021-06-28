#ifndef REGISTER
#define REGISTER

class Register{
public:
  int size_least_loaded_class;
  int mutnum_least_loaded_class;
  int mutnum_least_single;
  int mutnum_least_duplicated;
  double mean_mutnum_single;
  double mean_mutnum_duplicated;
  double mean_fitness;
  double mean_heterosite;
  double var_heterosite;
};

#endif
