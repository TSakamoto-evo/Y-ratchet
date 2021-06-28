#include "parameters.hpp"

void Parameters::set_sele_coef(){
  s.clear();
  t1.clear();
  t2.clear();

  std::random_device seed;
  std::mt19937 mt(seed());

  std::exponential_distribution<> dist(-1.0 / mean_sele);
  std::uniform_real_distribution<> dist2(0.0, 1.0);

  for(int i = 0; i < num_loci_single; i++){
    double add = -dist(mt);
    while(add < -1.0){
      add = -dist(mt);
    }

    s.push_back(add);
  }

  for(int i = 0; i < num_loci_duplicated; i++){
    double add = -dist(mt);
    while(add < -1.0){
      add = -dist(mt);
    }

    t2.push_back(add);
  }

  for(int i = 0; i < num_loci_duplicated; i++){
    double scale = dist(mt);
    t1.push_back(t2.at(i) * scale);
  }
}

void Parameters::set_sele_coef_externally(std::vector<double> input_s){
  s.clear();
  t1.clear();
  t2.clear();

  for(int i = 0; i < num_loci_single; i++){
    s.push_back(input_s.at(i));
  }
  for(int i = num_loci_single; i < num_loci_single + num_loci_duplicated; i++){
    t1.push_back(input_s.at(i) * dom);
    t2.push_back(input_s.at(i));
  }
}
