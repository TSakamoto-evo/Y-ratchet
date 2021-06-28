#include "population.hpp"

Population::Population(const Parameters input_paras){
  paras = input_paras;

  std::vector<int> tmp(paras.pop_size, 0);
  single = tmp;
  one_intact = tmp;
  zero_intact = tmp;

  std::random_device seed;
  std::mt19937 mt_tmp(seed());
  mt = mt_tmp;
}

void Population::selection(Register& regi){
  regi.size_least_loaded_class = 0;
  regi.mutnum_least_loaded_class = paras.num_loci_single + paras.num_loci_duplicated;
  regi.mutnum_least_single = paras.num_loci_single;
  regi.mutnum_least_duplicated = paras.num_loci_duplicated;
  regi.mean_mutnum_single = 0.0;
  regi.mean_mutnum_duplicated = 0.0;
  regi.mean_fitness = 0.0;
  regi.mean_heterosite = 0.0;
  regi.var_heterosite = 0.0;

  std::vector<double> fitness;
  for(int i = 0; i < paras.pop_size; i++){
    double tmp_fit = std::pow(1.0 - paras.sele_s, single.at(i)) *
                      std::pow(1.0 + paras.sele_t1, (paras.num_loci_duplicated - one_intact.at(i) - zero_intact.at(i))) *
                      std::pow(1.0 - paras.sele_t2, zero_intact.at(i));

    fitness.push_back(tmp_fit);

    if(single.at(i) + zero_intact.at(i) == regi.mutnum_least_loaded_class){
      regi.size_least_loaded_class++;
    }else if(single.at(i) + zero_intact.at(i) < regi.mutnum_least_loaded_class){
      regi.mutnum_least_loaded_class = single.at(i) + zero_intact.at(i);
      regi.size_least_loaded_class = 1;
    }

    if(single.at(i) < regi.mutnum_least_single){
      regi.mutnum_least_single = single.at(i);
    }

    if(zero_intact.at(i) < regi.mutnum_least_duplicated){
      regi.mutnum_least_duplicated = zero_intact.at(i);
    }

    regi.mean_mutnum_single += single.at(i);
    regi.mean_mutnum_duplicated += zero_intact.at(i);
    regi.mean_fitness += tmp_fit;
    regi.mean_heterosite += one_intact.at(i);
    regi.var_heterosite += one_intact.at(i) * one_intact.at(i);
  }
  regi.mean_mutnum_single /= paras.pop_size;
  regi.mean_mutnum_duplicated /= paras.pop_size;
  regi.mean_fitness /= paras.pop_size;
  regi.mean_heterosite /= paras.pop_size;
  regi.var_heterosite = regi.var_heterosite / paras.pop_size -
                          regi.mean_heterosite * regi.mean_heterosite;

  std::discrete_distribution<> dist(fitness.begin(), fitness.end());
  std::vector<int> next_single;
  std::vector<int> next_one_intact;
  std::vector<int> next_zero_intact;

  for(int i = 0; i < paras.pop_size; i++){
    int parent = dist(mt);
    next_single.push_back(single.at(parent));
    next_one_intact.push_back(one_intact.at(parent));
    next_zero_intact.push_back(zero_intact.at(parent));
  }

  single = next_single;
  one_intact = next_one_intact;
  zero_intact = next_zero_intact;
}

void Population::mutation(){
  std::poisson_distribution<> mut_gen((paras.pop_size * paras.num_loci_single +
                    2 * paras.pop_size * paras.num_loci_duplicated) * paras.mu);
  std::uniform_int_distribution<> choose_ind(0, paras.pop_size - 1);
  std::uniform_int_distribution<> choose_mut_locus(1, paras.num_loci_single +
                    2 * paras.num_loci_duplicated);

  int mut_num = mut_gen(mt);

  for(int i = 0; i < mut_num; i++){
    int ind = choose_ind(mt);
    int locus = choose_mut_locus(mt);

    if(locus <= paras.num_loci_single){
      if(locus > single.at(ind)){
        single.at(ind)++;
      }
    }else if(locus <= paras.num_loci_single + paras.num_loci_duplicated){
      locus -= paras.num_loci_single;
      if(locus > zero_intact.at(ind) + one_intact.at(ind)){
        one_intact.at(ind)++;
      }
    }else{
      locus -= (paras.num_loci_single + paras.num_loci_duplicated);
      if(locus > zero_intact.at(ind) && locus <= zero_intact.at(ind) + one_intact.at(ind)){
        one_intact.at(ind)--;
        zero_intact.at(ind)++;
      }else if(locus > zero_intact.at(ind) + one_intact.at(ind)){
        one_intact.at(ind)++;
      }
    }
  }
}

void Population::conversion(){
  std::poisson_distribution<> conv_gen(2 * paras.pop_size *
                                        paras.num_loci_duplicated * paras.conv);
  std::uniform_int_distribution<> choose_ind(0, paras.pop_size - 1);
  std::uniform_int_distribution<> choose_conv_locus(1, 2 * paras.num_loci_duplicated);

  int conv_num = conv_gen(mt);

  for(int i = 0; i < conv_num; i++){
    int ind = choose_ind(mt);
    int locus = choose_conv_locus(mt);

    if(locus <= paras.num_loci_duplicated){
      if(locus > zero_intact.at(ind) && locus <= zero_intact.at(ind) + one_intact.at(ind)){
        one_intact.at(ind) --;
        zero_intact.at(ind) ++;
      }
    }else{
      locus -= paras.num_loci_duplicated;
      if(locus > zero_intact.at(ind) && locus <= zero_intact.at(ind) + one_intact.at(ind)){
        one_intact.at(ind) --;
      }
    }
  }
}

void Population::one_generation(Register& regi){
  selection(regi);
  mutation();
  conversion();
}

void Population::erase_single(){
  for(auto& i: single){
    i--;
  }
}

void Population::erase_duplicated(){
  for(auto& i: zero_intact){
    i--;
  }
}
