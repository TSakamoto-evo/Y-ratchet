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

  int target_single = 1000;
  int target_duplicated = 1000;

  paras.pop_size = 500;
  paras.num_loci_single = target_single;
  paras.num_loci_duplicated = target_duplicated;
  paras.mu = 0.0002;

  paras.sele_s = 0.2;

  paras.sele_t2 = 0.2;
  paras.sele_t1 = 0.0;
  paras.conv = tmp_conv;

  int prerun_round = 10;
  int rep_prerun = 10;
  int max_rep = 100000;
  int success_count1 = 0;
  int success_count2 = 0;

  std::string file_index = argv[1];
  std::string txt = ".txt";

  std::string string0 = "parameters";
  std::string string5 = "click_time_single";
  std::string string6 = "click_time_duplicated";
  std::string string9 = "progress";

  std::ofstream ofs0(string0 + file_index + txt);
  std::ofstream ofs5(string5 + file_index + txt);
  std::ofstream ofs6(string6 + file_index + txt);
  std::ofstream ofs9(string9 + file_index + txt);

  {
    ofs0 << "Population size: " << paras.pop_size << std::endl;
    ofs0 << "Target loci single: " << paras.num_loci_single << std::endl;
    ofs0 << "Target loci duplicated: " << paras.num_loci_duplicated << std::endl;
    ofs0 << "mu: " << paras.mu << std::endl;
    ofs0 << "Selection s: " << paras.sele_s << std::endl;
    ofs0 << "Selection t1: " << paras.sele_t1 << std::endl;
    ofs0 << "Selection t2: " << paras.sele_t2 << std::endl;
    ofs0 << "Conversion: " << paras.conv << std::endl;
  }

  //prerun
  {
    for(int pre_round = 0; pre_round < prerun_round; pre_round++){
      int total_single = 0;
      int total_duplicated = 0;

      for(int i = 0; i < rep_prerun; i++){
        Population pop(paras);
        Register regi;

        for(int j = 0; j < 100 * paras.pop_size; j++){
          pop.one_generation(regi);
        }

        total_single += regi.mutnum_least_single;
        total_duplicated += regi.mutnum_least_duplicated;
      }

      int loss_single = total_single / rep_prerun;
      int loss_duplicated = total_duplicated / rep_prerun;

      double margin_single = total_single * 1.0 / rep_prerun - loss_single;
      double margin_duplicated = total_duplicated * 1.0 / rep_prerun - loss_duplicated;

      if(margin_single < margin_duplicated){
        loss_duplicated++;
      }else{
        loss_single++;
      }

      ofs9 << pre_round << "\t" << loss_single << "\t" << paras.num_loci_single <<
      "\t" << loss_duplicated << "\t" << paras.num_loci_duplicated << std::endl;

      paras.num_loci_single = (target_single + loss_single);
      paras.num_loci_duplicated = (target_duplicated + loss_duplicated);

      if(pre_round == prerun_round - 1){
        ofs0 << "add single: " << loss_single << std::endl;
        ofs0 << "add_duplicated: " << loss_duplicated << std::endl;
      }
    }
  }
  ofs9 << std::endl;


  for(int i = 0; i < max_rep; i++){
    Population pop(paras);
    Register regi;
    regi.mutnum_least_loaded_class = 0;
    regi.mutnum_least_single = 0;
    regi.mutnum_least_duplicated = 0;

    int j = 0;
    int regi_time = 0;

    while(1){
      int ini_min_single_class_num = regi.mutnum_least_single;
      int ini_min_duplicated_class_num = regi.mutnum_least_duplicated;
      pop.one_generation(regi);

      if(ini_min_single_class_num != regi.mutnum_least_single){
        if(target_single == paras.num_loci_single - ini_min_single_class_num &&
           target_duplicated == paras.num_loci_duplicated - ini_min_duplicated_class_num){
          ofs5 << i << "\t" << j - regi_time << std::endl;
          success_count1++;
        }
      }

      if(ini_min_duplicated_class_num != regi.mutnum_least_duplicated){
        if(target_single == paras.num_loci_single - ini_min_single_class_num &&
           target_duplicated == paras.num_loci_duplicated - ini_min_duplicated_class_num){
          ofs6 << i << "\t" << j - regi_time << std::endl;
          success_count2++;
        }
      }

      if(regi.mutnum_least_single == paras.num_loci_single - target_single &&
         regi.mutnum_least_duplicated == paras.num_loci_duplicated - target_duplicated &&
         regi_time == 0){
        regi_time = j;
      }

      if(target_single > paras.num_loci_single - regi.mutnum_least_single ||
         target_duplicated > paras.num_loci_duplicated - regi.mutnum_least_duplicated){

        ofs9 << i << "\t" << j << "\t" << paras.num_loci_single - regi.mutnum_least_single << "\t" << paras.num_loci_duplicated - regi.mutnum_least_duplicated << std::endl;
        break;
      }

      /*
      if(j % 10000 == 0){
        ofs9 << i << "\t" << j << "\t" << paras.num_loci_single - regi.mutnum_least_single << "\t" << paras.num_loci_duplicated - regi.mutnum_least_duplicated << std::endl;
      }
      */

      j++;
    }

    if(success_count1 + success_count2 >= 2000){
      break;
    }
  }

  return(0);
}
