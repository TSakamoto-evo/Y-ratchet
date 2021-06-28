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
  paras.conv = tmp_conv;
  paras.dom = 1.0;
  //paras.mean_sele = -0.015;

  std::ifstream ifs_s("../dist_s.txt");

  std::vector<double> input_s;

  std::string str;
  if (ifs_s.fail()){
    std::cerr << "Failed to open file." << std::endl;
    return(-1);
  }
  while(getline(ifs_s, str)){
    double tmp = std::stod(str);
    input_s.push_back(tmp);
  }

  paras.set_sele_coef_externally(input_s);

  std::string file_index = argv[1];
  std::string txt = ".txt";

  std::string string0 = "parameters";
  std::string string1 = "fit_change";
  std::string string2 = "regi_s_";
  std::string string3 = "regi_t1_";
  std::string string4 = "regi_t2_";
  std::string string5 = "fixed_single";
  std::string string6 = "fixed_dup";

  std::string string7 = "progress";

  std::ofstream ofs0(string0 + file_index + txt);
  std::ofstream ofs1(string1 + file_index + txt);
  std::ofstream ofs2(string2 + file_index + txt);
  std::ofstream ofs3(string3 + file_index + txt);
  std::ofstream ofs4(string4 + file_index + txt);
  std::ofstream ofs5(string5 + file_index + txt);
  std::ofstream ofs6(string6 + file_index + txt);
  std::ofstream ofs7(string7 + file_index + txt);

  {
    ofs0 << "Population size: " << paras.pop_size << std::endl;
    ofs0 << "Loci single: " << paras.num_loci_single << std::endl;
    ofs0 << "Loci duplicated: " << paras.num_loci_duplicated << std::endl;
    ofs0 << "mu: " << paras.mu << std::endl;
    //ofs0 << "mean_sele: " << paras.mean_sele << std::endl;
    ofs0 << "Conversion: " << paras.conv << std::endl;
  }

  {
    ofs2 << "selection" << std::endl;
    std::vector<double> tmp_s = paras.return_s();
    for(const auto& i: tmp_s){
      ofs2 << i << std::endl;
    }

    ofs3 << "selection" << std::endl;
    std::vector<double> tmp_t1 = paras.return_t1();
    for(const auto& i: tmp_t1){
      ofs3 << i << std::endl;
    }

    ofs4 << "selection" << std::endl;
    std::vector<double> tmp_t2 = paras.return_t2();
    for(const auto& i: tmp_t2){
      ofs4 << i << std::endl;
    }
  }

  Population pop(paras);
  Register regi;
  regi.log_max_fitness = 0.0;
  regi.mean_fitness = 0.0;

  for(int i = 0; i < 1000 * paras.pop_size; i++){
    pop.one_generation(regi);

    if(i % 100 == 0){
      ofs7 << i << std::endl;
    }

    if(i % 10000 == 0){
      ofs1 << i << "\t" << regi.mean_fitness << "\t" << std::exp(regi.log_max_fitness) << std::endl;

      std::vector<int> tmp_fix_single = pop.return_irreversible_single();
      std::vector<int> tmp_fix_dup = pop.return_irreversible_duplicated();

      ofs5 << i << std::flush;
      for(const auto& j: tmp_fix_single){
        ofs5 << "\t" << j << std::flush;
      }
      ofs5 << std::endl;

      ofs6 << i << std::flush;
      for(const auto& j: tmp_fix_dup){
        ofs6 << "\t" << j << std::flush;
      }
      ofs6 << std::endl;
    }
  }

  return(0);
}
