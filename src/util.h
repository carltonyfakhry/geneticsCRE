#ifndef GCRE_UTIL_H
#define GCRE_UTIL_H

#include <unistd.h>
#include <cstdint>
#include <limits>
#include <vector>
#include <chrono>
#include <ctime>

class Timer {

public:

  static void print_header() {
    std::printf("\nTIME:PID IMPL METHOD WIDTH LENGTH PATHS PERMS THREADS CYCLES MS\n\n");
  }

  Timer(const JoinExec& exec, int path_length, uint64_t total_paths) : exec(exec), path_length(path_length), total_paths(total_paths) {}

  ~Timer(){
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
    std::printf("\nTIME:%d %s m%d %d %d %lu %d %d %lu %lu\n\n", getpid(), gs_instr_label.c_str(), exec.method, exec.width_ul * 64,
      path_length, total_paths, exec.iterations, exec.nthreads, std::clock() - clock_start, (unsigned long) time.count());
  }

private:

  const clock_t clock_start = std::clock();
  const std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
  const JoinExec& exec;
  const int path_length;
  const uint64_t total_paths;

};

#endif