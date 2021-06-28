#ifndef POPULATION
#define POPULATION

#include <vector>
#include <random>
#include <cmath>
#include <iostream>
#include "parameters.hpp"
#include "register.hpp"
#include "individual.hpp"

class Population{
private:
  Parameters paras;
  std::vector<Individual> pop;
  std::mt19937 mt;

  std::vector<double> log_single0;
  std::vector<double> log_dup2;
  std::vector<double> log_dup0;

public:
  Population(const Parameters input_paras);
  void selection(Register& regi);
  void mutation();
  void conversion();
  void one_generation(Register& regi);
  void erase_single();
  void erase_duplicated();

  std::vector<int> return_irreversible_single();
  std::vector<int> return_irreversible_duplicated();
};

#endif
