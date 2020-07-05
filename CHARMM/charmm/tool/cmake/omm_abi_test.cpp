#include<iostream>
#include<OpenMM.h>

int main() {
  OpenMM::Platform::getPlatformByName("CUDA");
  std::cout << "Hi there!" << std::endl;
}


