#include "population.hpp"

Population::Population(const Parameters input_paras){
  paras = input_paras;

  std::vector<Individual> tmp;
  for(int i = 0; i < paras.pop_size; i++){
    tmp.emplace_back(paras.num_loci_single, paras.num_loci_duplicated);
  }
  pop = tmp;

  std::random_device seed;
  std::mt19937 mt_tmp(seed());
  mt = mt_tmp;

  log_single0.clear();
  log_dup2.clear();
  log_dup0.clear();

  for(int i = 0; i < paras.num_loci_single; i++){
    log_single0.push_back(std::log(1.0 - paras.s.at(i)));
  }
  for(int i = 0; i < paras.num_loci_duplicated; i++){
    log_dup2.push_back(std::log(1.0 + paras.t1.at(i)));
    log_dup0.push_back(std::log(1.0 - paras.t2.at(i)));
  }
}

void Population::selection(Register& regi){

  std::vector<double> fitness;
  double log_fit_base = regi.log_max_fitness;
  regi.mean_fitness = 0.0;

  for(int i = 0; i < paras.pop_size; i++){
    double tmp_log_fit = pop.at(i).return_log_fitness(paras.num_loci_single,
      paras.num_loci_duplicated, log_single0, log_dup2, log_dup0);

    if(i == 0){
      regi.log_max_fitness = tmp_log_fit;
    }else if(tmp_log_fit > regi.log_max_fitness){
      regi.log_max_fitness = tmp_log_fit;
    }

    double tmp_fit = std::exp(tmp_log_fit - log_fit_base);
    fitness.push_back( tmp_fit );
    regi.mean_fitness += tmp_fit;
  }
  regi.mean_fitness /= paras.pop_size;
  regi.mean_fitness *= std::exp(log_fit_base);

  std::discrete_distribution<> dist(fitness.begin(), fitness.end());

  std::vector<Individual> next_gen;
  for(int i = 0; i < paras.pop_size; i++){
    int parent = dist(mt);
    next_gen.push_back(pop.at(parent));
  }

  pop = next_gen;
}

void Population::mutation(){
  std::poisson_distribution<> mut_gen((paras.pop_size * paras.num_loci_single +
                    2 * paras.pop_size * paras.num_loci_duplicated) * paras.mu);
  std::uniform_int_distribution<> choose_ind(0, paras.pop_size - 1);
  std::uniform_int_distribution<> choose_mut_locus(0, paras.num_loci_single +
                    2 * paras.num_loci_duplicated - 1);

  int mut_num = mut_gen(mt);

  for(int i = 0; i < mut_num; i++){
    int chosen_ind = choose_ind(mt);
    int chosen_locus = choose_mut_locus(mt);

    int locus_class = 0;
    int locus_index = 0;

    if(chosen_locus < paras.num_loci_single){
      locus_class = 0;
      locus_index = chosen_locus;
    }else if(chosen_locus < paras.num_loci_single + paras.num_loci_duplicated){
      locus_class = 1;
      locus_index = chosen_locus - paras.num_loci_single;
    }else{
      locus_class = 2;
      locus_index = chosen_locus - paras.num_loci_single - paras.num_loci_duplicated;
    }

    pop.at(chosen_ind).mutation(locus_class, locus_index);
  }
}

void Population::conversion(){
  std::poisson_distribution<> conv_gen(2 * paras.pop_size *
                                        paras.num_loci_duplicated * paras.conv);
  std::uniform_int_distribution<> choose_ind(0, paras.pop_size - 1);
  std::uniform_int_distribution<> choose_conv_locus(0, 2 * paras.num_loci_duplicated - 1);

  int conv_num = conv_gen(mt);

  for(int i = 0; i < conv_num; i++){
    int chosen_ind = choose_ind(mt);
    int chosen_locus = choose_conv_locus(mt);

    int locus_class = 0;
    int locus_index = 0;

    if(chosen_locus < paras.num_loci_duplicated){
      locus_class = 0;
      locus_index = chosen_locus;
    }else{
      locus_class = 1;
      locus_index = chosen_locus - paras.num_loci_duplicated;
    }

    pop.at(chosen_ind).conversion(locus_class, locus_index);
  }
}

void Population::one_generation(Register& regi){
  selection(regi);
  mutation();
  conversion();
}

std::vector<int> Population::return_irreversible_single(){
  std::vector<int> ret;

  for(int i = 0; i < paras.num_loci_single; i++){
    for(int j = 0; j < paras.pop_size; j++){
      if(pop.at(j).return_single_state(i) != 1){
        break;
      }
      if(j == paras.pop_size - 1){
        ret.push_back(i);
      }
    }
  }

  return(ret);
}

std::vector<int> Population::return_irreversible_duplicated(){
  std::vector<int> ret;

  for(int i = 0; i < paras.num_loci_duplicated; i++){
    for(int j = 0; j < paras.pop_size; j++){
      if(pop.at(j).return_dup_state(i) != 1){
        break;
      }
      if(j == paras.pop_size - 1){
        ret.push_back(i);
      }
    }
  }

  return(ret);
}
