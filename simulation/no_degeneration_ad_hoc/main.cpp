#include <iostream>
#include <fstream>
#include <string>
#include "population.hpp"
#include "parameters.hpp"

int main(int argc, char *argv[]){
  double tmp_conv;

  if(argc == 3){
    sscanf(argv[2], "%lf", &tmp_conv);
  }else{
    std::cout << "Error! The number of parameters is different!" << std::endl;
    return(1);
  }

  Parameters paras;

  paras.pop_size = 10000;
  paras.num_loci_single = 1000;
  paras.num_loci_duplicated = 1000;
  paras.mu = 0.00001;

  paras.sele_s = 0.01;

  paras.sele_t2 = 0.01;
  paras.sele_t1 = 0.0;
  paras.conv = tmp_conv;

  int max_rep = 1;

  std::string file_index = argv[1];
  std::string txt = ".txt";

  std::string string0 = "parameters";
  std::string string1 = "fit_change";
  std::string string4 = "click_time";
  std::string string5 = "click_time_single";
  std::string string6 = "click_time_duplicated";
  std::string string7 = "progress";

  std::ofstream ofs0(string0 + file_index + txt);
  std::ofstream ofs1(string1 + file_index + txt);
  std::ofstream ofs4(string4 + file_index + txt);
  std::ofstream ofs5(string5 + file_index + txt);
  std::ofstream ofs6(string6 + file_index + txt);
  std::ofstream ofs7(string7 + file_index + txt);

  {
    ofs0 << "Population size: " << paras.pop_size << std::endl;
    ofs0 << "Loci single: " << paras.num_loci_single << std::endl;
    ofs0 << "Loci duplicated: " << paras.num_loci_duplicated << std::endl;
    ofs0 << "mu: " << paras.mu << std::endl;
    ofs0 << "Selection s: " << paras.sele_s << std::endl;
    ofs0 << "Selection t1: " << paras.sele_t1 << std::endl;
    ofs0 << "Selection t2: " << paras.sele_t2 << std::endl;
    ofs0 << "Conversion: " << paras.conv << std::endl;
  }


  for(int i = 0; i < max_rep; i++){
    Population pop(paras);
    int click = 0;
    Register regi;
    regi.mutnum_least_loaded_class = 0;
    regi.mutnum_least_single = 0;
    regi.mutnum_least_duplicated = 0;

    for(int j = -100 * paras.pop_size; j < 10000 * paras.pop_size; j++){
      int ini_min_class_num = regi.mutnum_least_loaded_class;
      int ini_min_single_class_num = regi.mutnum_least_single;
      int ini_min_duplicated_class_num = regi.mutnum_least_duplicated;
      pop.one_generation(regi);

      if(ini_min_class_num != regi.mutnum_least_loaded_class){
        ofs4 << i << "\t" << j << "\t" << regi.mutnum_least_loaded_class << std::endl;
        if(j >= 0){
          click += (regi.mutnum_least_loaded_class - ini_min_class_num);
        }
      }

      if(ini_min_single_class_num != regi.mutnum_least_single){
        ofs5 << i << "\t" << j << "\t" << regi.mutnum_least_single << std::endl;
        int diff_single = regi.mutnum_least_single - ini_min_single_class_num;

        for(int k = 0; k < diff_single; k++){
          pop.erase_single();
          regi.mutnum_least_loaded_class --;
          regi.mutnum_least_single --;
        }
      }

      if(ini_min_duplicated_class_num != regi.mutnum_least_duplicated){
        ofs6 << i << "\t" << j << "\t" << regi.mutnum_least_duplicated << std::endl;
        int diff_dup = regi.mutnum_least_duplicated - ini_min_duplicated_class_num;

        for(int k = 0; k < diff_dup; k++){
          pop.erase_duplicated();
          regi.mutnum_least_loaded_class--;
          regi.mutnum_least_duplicated--;
        }
      }

      if(j % 10000 == 0){
        ofs7 << j << std::endl;
      }

      if(j % 10000 == 0 && j > 0){
        ofs1 << regi.size_least_loaded_class << "\t" <<
        regi.mutnum_least_loaded_class << "\t" <<
        regi.mean_mutnum_single << "\t" <<
        regi.mean_mutnum_duplicated << "\t" <<
        regi.mean_fitness << "\t" <<
        regi.mean_heterosite << "\t" <<
        regi.var_heterosite << std::endl;
      }

      if(click == 10000){
        break;
      }
    }
  }

  return(0);
}
