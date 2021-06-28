#include <iostream>
#include <fstream>
#include <random>

int main(){
  std::random_device seed;
  std::mt19937 mt(seed());

  double mean_s = 0.0025;
  int num_loci = 2000;

  std::ofstream ofs("dist_s.txt");

  std::exponential_distribution<> dist(1.0 / mean_s);

  for(int i = 0; i < num_loci; i++){
    ofs << dist(mt) << std::endl;
  }

  return(0);
}
