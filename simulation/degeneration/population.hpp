#ifndef POPULATION
#define POPULATION

#include <vector>
#include <random>
#include <cmath>
#include <iostream>
#include "parameters.hpp"
#include "register.hpp"

class Population{
private:
  Parameters paras;
  std::vector<int> single;
  std::vector<int> one_intact;
  std::vector<int> zero_intact;
  std::mt19937 mt;

public:
  Population(const Parameters input_paras);
  void selection(Register& regi);
  void mutation();
  void conversion();
  void one_generation(Register& regi);
  void erase_single();
  void erase_duplicated();
};

#endif
