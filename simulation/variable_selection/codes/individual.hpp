#ifndef INDIVIDUAL
#define INDIVIDUAL

#include <vector>

class Individual{
private:
  std::vector<bool> single;
  std::vector<bool> dup1;
  std::vector<bool> dup2;

public:
  Individual(const int num_single, const int num_dup);
  double return_log_fitness(const int num_single, const int num_dup,
                            const std::vector<double>& log_single0,
                            const std::vector<double>& log_dup2,
                            const std::vector<double>& log_dup0) const;
  void mutation(const int locus_class, const int locus_index);
  void conversion(const int locus_class, const int locus_index);

  bool return_single_state(const int locus_index){
    return(single.at(locus_index));
  };
  bool return_dup_state(const int locus_index){
    return(dup1.at(locus_index) && dup2.at(locus_index));
  };
};

#endif
