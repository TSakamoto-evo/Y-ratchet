#include "individual.hpp"

Individual::Individual(const int num_single, const int num_dup){
  std::vector<bool> tmp_single(num_single, 0);
  std::vector<bool> tmp_dup(num_dup, 0);

  single = tmp_single;
  dup1 = tmp_dup;
  dup2 = tmp_dup;
}

double Individual::return_log_fitness(const int num_single, const int num_dup,
                                      const std::vector<double>& log_single0,
                                      const std::vector<double>& log_dup2,
                                      const std::vector<double>& log_dup0)
                                      const{
  double log_fitness = 0;
  for(int i = 0; i < num_single; i++){
    if(single.at(i)){
      log_fitness += log_single0.at(i);
    }
  }

  for(int i = 0; i < num_dup; i++){
    if(dup1.at(i) && dup2.at(i)){
      log_fitness += log_dup0.at(i);
    }else if(!dup1.at(i) && !dup2.at(i)){
      log_fitness += log_dup2.at(i);
    }
  }

  return(log_fitness);
}

void Individual::mutation(const int locus_class, const int locus_index){
  if(locus_class == 0){
    single.at(locus_index) = 1;
  }else if(locus_class == 1){
    dup1.at(locus_index) = 1;
  }else{
    dup2.at(locus_index) = 1;
  }
}

void Individual::conversion(const int locus_class, const int locus_index){
  if(locus_class == 0){
    dup1.at(locus_index) = dup2.at(locus_index);
  }else{
    dup2.at(locus_index) = dup1.at(locus_index);
  }
}
