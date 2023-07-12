#include <cstddef>

#include <algorithm>
#include <atomic>
#include <optional>
#include <utility>

#include <parlay/primitives.h>
#include <parlay/sequence.h>

// **************************************************************
// A simple concurrent hash-based map
// Supports concurrent linearizable, insert, find and remove.
// size(), keys() don't linearize with updates.
// Requires the capacity to be  specified on construction.
// No more than the capacity distinct keys can ever be added.
// Once a key is added, removing it will empty the value and mark
// the key as deleted, but only a value with the same key can use the
// same slot (i.e. it still counts towards the capacity).
// It uses locks, but holds them very briefly.
// **************************************************************
constexpr size_t kResizableTableCacheLineSz = 128;
template <typename K,
	  typename V,
	  typename Hash = parlay::hash<K>>
struct hash_table {
 public:
  using KV = std::pair<K,V>;
  using index = unsigned long;
  long m;
  long mask;
  KV empty;
  index first_index(K k) { return hash(k) & mask;}
  index next_index(index h) { return (h + 1)& mask;}
  Hash hash;
  bool equal(KV a, KV b){return 
    parlay::internal::get_key(a)==parlay::internal::get_key(b) &&
    parlay::internal::get_val(a)==parlay::internal::get_val(b);}
  inline bool atomic_compare_and_swap(KV *a, KV oldval, KV newval){
    uint64_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(KV));
    std::memcpy(&r_nval, &newval, sizeof(KV));
    return __sync_bool_compare_and_swap(reinterpret_cast<uint64_t *>(a), r_oval,
                                        r_nval);
  }

  parlay::sequence<KV> H;
  parlay::sequence<size_t> cts;
  size_t ne;
  bool overfull;
  size_t size(){
    for (int i = 0;i<parlay::num_workers(); i++){
        ne+= cts[i*kResizableTableCacheLineSz];
        cts[i*parlay::num_workers()]=0;
    }
    if (ne >= m){overfull=true;}
    return ne;
  }

  hash_table(long size,  KV _empty, Hash&& hash = {}) 
    : m((size_t)1 << parlay::log2_up((size_t)(2 * size))),
      mask(m-1), overfull(false), ne(0), empty(_empty),
      hash(hash), 
      H(parlay::sequence<KV>(m)),
      cts(parlay::sequence<size_t>(parlay::num_workers()*kResizableTableCacheLineSz)) {
        init();
      }
  void init(){
    parlay::parallel_for(0, m, [&](size_t i){H[i]=empty;});
    for (size_t i = 0; i<parlay::num_workers();i++){
        cts[i*kResizableTableCacheLineSz]=0;
    }
  }

  bool insert(const K& k, const V& v) {
    KV kv = std::make_pair(k,v);
    index i = first_index(k);
    long count = 0;
    for (count; count < 2000l; count++) {
      if (equal(H[i], empty) && atomic_compare_and_swap(&H[i],empty,kv)){
        size_t wn = parlay::worker_id();
        cts[wn * kResizableTableCacheLineSz]++;
        return true;
      }
      if (equal(H[i], kv)) return false;
	  i = next_index(i);
	}
    std::cout << "Hash table overfull" << std::endl;
    overfull=true;
    return false;
  }

  bool contain(const KV kv) {
    K k = parlay::internal::get_key(kv);
    index i = first_index(k);
    while (true) {
      if (equal(H[i], empty)) return false;
      if (equal(H[i], kv)) return true;
      i = next_index(i);
    }
  }

  template <class F>
  void map(F& f) {
    parlay::parallel_for(0, m, [&](size_t i) {
      if (!equal(H[i], empty)) {
        f(H[i]);
      }
    });
  }

};

template <typename K,
    typename V, typename T>
struct iter_k{
    using KV = std::pair<K,V>;
    using index = unsigned long;
    K k;
    index i;
    size_t num_prob;
    T& table; 
    iter_k(K _k, T& _table):
        k(_k), table(_table) {}
    bool equal(KV a, KV b){return 
      parlay::internal::get_key(a)==parlay::internal::get_key(b) &&
      parlay::internal::get_val(a)==parlay::internal::get_val(b);}
    bool init(){
        i = table.first_index(k);
        while (true){
            if (equal(table.H[i],table.empty)){
                return false;
            }
            if (parlay::internal::get_key(table.H[i]) == k){
                num_prob=0;
                return true;
            }
            i = table.next_index(i);
        }
    }

    bool has_next(){
      // printf("        key: %d, i: %d, num_prob: %d\n", k, i, num_prob);
        while (num_prob<table.m){
            i = table.next_index(i);
            num_prob++;
            if (equal(table.H[i],table.empty)){
                return false;
            }
            if (parlay::internal::get_key(table.H[i]) ==k){
                return true;
            }
        }
        return false;
    }

    V get(){
        return parlay::internal::get_val(table.H[i]);
    }
};