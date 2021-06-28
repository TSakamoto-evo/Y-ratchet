#ifndef PARAMETERS
#define PARAMETERS

#include <vector>
#include <random>

class Parameters{
public:
  int pop_size;
  int num_loci_single;
  int num_loci_duplicated;
  double mu;
  double dom;

  double mean_sele;
  double conv;

  std::vector<double> s;
  std::vector<double> t1;
  std::vector<double> t2;

  void set_sele_coef();
  std::vector<double> return_s() const{ return(s); };
  std::vector<double> return_t1() const{ return(t1); };
  std::vector<double> return_t2() const{ return(t2); };
  void set_sele_coef_externally(std::vector<double> input_s);
};

#endif
