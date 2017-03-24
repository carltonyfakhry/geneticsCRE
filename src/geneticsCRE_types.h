#ifndef GCRE_TYPES_H
#define GCRE_TYPES_H

struct paths_type {};

struct paths_base : paths_type {
  int size = 0;
  int width_ul = 0;
  int num_cases = 0;
};

struct join_config {
  int num_cases = 0;
  int num_controls = 0;
  int top_k = 0;
  int path_length = 0;
  int iterations = 0;
  int nthreads = 0;
};

// TODO must be a better way to define overrides
class JoinMethod {
public:
  virtual paths_type* createPathSet() const = 0;
  virtual paths_type* createPathSet(vec2d_i& data, int num_cases, int num_controls) const = 0;
  virtual paths_type* createPathSet(vec2d_u64& pos, vec2d_u64& neg, int num_cases, int num_controls) const = 0;
  virtual joined_res join(join_config& conf, vector<uid_ref>& uids, vector<int>& join_gene_signs, vec2d_d& value_table, vec2d_u16& permute_cases, paths_type* p_paths0, paths_type* p_paths1, paths_type* p_paths_res, uint64_t total_paths) const = 0;
};

class JoinMethod1Native : public JoinMethod {
public:
  virtual paths_type* createPathSet() const;
  virtual paths_type* createPathSet(vec2d_i& data, int num_cases, int num_controls) const;
  virtual paths_type* createPathSet(vec2d_u64& pos, vec2d_u64& neg, int num_cases, int num_controls) const;
  virtual joined_res join(join_config& conf, vector<uid_ref>& uids, vector<int>& join_gene_signs, vec2d_d& value_table, vec2d_u16& permute_cases, paths_type* p_paths0, paths_type* p_paths1, paths_type* p_paths_res, uint64_t total_paths) const;
};

class JoinMethod2Native : public JoinMethod {
public:
  virtual paths_type* createPathSet() const;
  virtual paths_type* createPathSet(vec2d_i& data, int num_cases, int num_controls) const;
  virtual paths_type* createPathSet(vec2d_u64& pos, vec2d_u64& neg, int num_cases, int num_controls) const;
  virtual joined_res join(join_config& conf, vector<uid_ref>& uids, vector<int>& join_gene_signs, vec2d_d& value_table, vec2d_u16& permute_cases, paths_type* p_paths0, paths_type* p_paths1, paths_type* p_paths_res, uint64_t total_paths) const;
};

#endif
