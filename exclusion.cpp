// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>

using namespace Rcpp;
using namespace RcppParallel;

// ---------------- Helpers ----------------

// central exclusion window as bases (1-based in spec, 0-based internally)
inline void compute_exclusion_base_range(int len, int exclusion_half_width_bases, int &E0, int &E1){
    if (exclusion_half_width_bases <= 0){ E0 = 1; E1 = 0; return; } // empty interval
    int center_1b = (len + 1) / 2;
    int S1 = std::max(1,   center_1b - exclusion_half_width_bases);
    int S2 = std::min(len, center_1b + exclusion_half_width_bases);
    E0 = S1 - 1; E1 = S2 - 1; // 0-based inclusive
}

// map an index in the "spliced" string (gap removed) back to original
inline int map_spliced_to_original(int u_spliced, int E0, int gapLen){
    // indices < E0 stay; indices >= E0 shift right by gapLen
    return (u_spliced < E0) ? u_spliced : (u_spliced + gapLen);
}

// char accessor over the spliced view
inline char at_spliced(const std::vector<char>& s, int u_spliced, int E0, int gapLen){
    return s[ map_spliced_to_original(u_spliced, E0, gapLen) ];
}

// Hamming over spliced ranges 
inline int hamming_distance_spliced(const std::vector<char>& A, int a0,
                                    const std::vector<char>& B, int b0,
                                    int k, int E0, int gapLen){
    int d = 0;
    for (int off = 0; off < k; ++off){
        if (at_spliced(A, a0 + off, E0, gapLen) != at_spliced(B, b0 + off, E0, gapLen)) ++d;
    }
    return d;
}

// Finds min distance between k-mers at equivalent positions (same start)
struct MinHammingResult {
    int position; // 1-based index for R compatibility
    int distance;
};

MinHammingResult find_min_hamming_equivalent_cpp_spliced(
    const std::vector<char>& A,
    const std::vector<char>& B,
    int k,
    int len,
    int E0,
    int E1
) {
    const int gapLen = (E1 >= E0) ? (E1 - E0 + 1) : 0;
    const int Ls = len - gapLen; // spliced length

    if (k <= 0 || k > Ls) {
        // match original guard style: return a benign pair
        return {1, 0};
    }

    int min_dist = std::numeric_limits<int>::max();
    int best_pos_spliced = 0; // 0-based spliced

    const int n_kmers = Ls - k + 1;
    for (int t = 0; t < n_kmers; ++t){
        int d = hamming_distance_spliced(A, t, B, t, k, E0, gapLen);
        if (d < min_dist){
            min_dist = d;
            best_pos_spliced = t;
        }
    }

    // map back to original index (0-based), then to 1-based
    int start_original0 = map_spliced_to_original(best_pos_spliced, E0, gapLen);
    return { start_original0 + 1, min_dist };
}

// ---------------- Worker ----------------
// Parallelization

struct DissimilarityWorkerRevised : public Worker {
    // Inputs
    const std::vector<std::vector<char>>& seq_chars_cpp;
    const int k;
    const int window_size;
    const int len;
    const int n;
    const int E0;
    const int E1;

    // Outputs
    RMatrix<double> basic_diss;
    RMatrix<double> window_diss;
    RMatrix<int>    min_positions;

    DissimilarityWorkerRevised(
        const std::vector<std::vector<char>>& seq_chars_cpp,
        int k, int window_size, int len, int n, int E0, int E1,
        NumericMatrix basic_diss_r, NumericMatrix window_diss_r, IntegerMatrix min_positions_r
    )
    : seq_chars_cpp(seq_chars_cpp), k(k), window_size(window_size), len(len), n(n), E0(E0), E1(E1),
      basic_diss(basic_diss_r), window_diss(window_diss_r), min_positions(min_positions_r) {}

    void operator()(std::size_t begin, std::size_t end) {
        const int gapLen = (E1 >= E0) ? (E1 - E0 + 1) : 0;
        const int Ls = len - gapLen; // spliced length

        for (std::size_t i = begin; i < end; ++i) {
            const std::vector<char>& chars_i = seq_chars_cpp[i];

            for (int j = (int)i + 1; j < n; ++j) {
                const std::vector<char>& chars_j = seq_chars_cpp[j];

                // ---------- BASIC (same-start on spliced coordinate) ----------
                MinHammingResult res = find_min_hamming_equivalent_cpp_spliced(chars_i, chars_j, k, len, E0, E1);
                basic_diss(i, j) = basic_diss(j, i) = res.distance;
                int min_pos_1based = res.position;
                min_positions(i, j) = min_positions(j, i) = min_pos_1based;

                // ---------- WINDOW (symmetric with offsets s = 0..window_size-1 in both directions, spliced) ----------
                int min_window_dist = std::numeric_limits<int>::max();

                // Track the positions that achieve the global min (spliced starts)
                int pos_i_best_sp = 0;
                int pos_j_best_sp = 0;
                bool have_best = false;

                if (k <= Ls && window_size > 0) {
                    // i vs j (positive shifts)
                    for (int s = 0; s < window_size; ++s) {
                        int max_t = (Ls - k) - s;
                        if (max_t < 0) break;

                        for (int t = 0; t <= max_t; ++t) {
                            int d = hamming_distance_spliced(chars_i, t, chars_j, t + s, k, E0, gapLen);
                            if (d < min_window_dist) {
                                min_window_dist = d;
                                pos_i_best_sp = t;
                                pos_j_best_sp = t + s;
                                have_best = true;
                            }
                        }
                    }

                    // j vs i (positive shifts)
                    for (int s = 0; s < window_size; ++s) {
                        int max_t = (Ls - k) - s;
                        if (max_t < 0) break;

                        for (int t = 0; t <= max_t; ++t) {
                            int d = hamming_distance_spliced(chars_j, t, chars_i, t + s, k, E0, gapLen);
                            if (d < min_window_dist) {
                                min_window_dist = d;
                                pos_i_best_sp = t + s;
                                pos_j_best_sp = t;
                                have_best = true;
                            }
                        }
                    }

                } else if (k <= Ls) {
                    // fallback: window not usable, reuse basic distance and position
                    min_window_dist = res.distance;
                    pos_i_best_sp = 0;
                    pos_j_best_sp = 0;
                    have_best = true;
                } else {
                    // k > spliced length: degenerate
                    min_window_dist = 0;
                    pos_i_best_sp = 0;
                    pos_j_best_sp = 0;
                    have_best = true;
                }

                if (!have_best || min_window_dist == std::numeric_limits<int>::max()) {
                    // final guard
                    min_window_dist = (k > 0 && k <= Ls) ? res.distance : 0;
                    pos_i_best_sp = 0;
                    pos_j_best_sp = 0;
                }

                // penalty is |pos_i - pos_j| (absolute difference on original coordinates)
                int pos_i_best_orig0 = map_spliced_to_original(pos_i_best_sp, E0, gapLen);
                int pos_j_best_orig0 = map_spliced_to_original(pos_j_best_sp, E0, gapLen);
                int penalty = std::abs(pos_i_best_orig0 - pos_j_best_orig0);

                double adjusted;
                if (window_size >= 2) {
                    adjusted = static_cast<double>(min_window_dist) * static_cast<double>(window_size - 1)
                               + static_cast<double>(penalty);
                } else {
                    // window_size == 0 or 1: no offset allowed, keep raw distance
                    adjusted = static_cast<double>(min_window_dist);
                }

                window_diss(i, j) = window_diss(j, i) = adjusted;
                // ---------- end window ----------
            }
        }
    }
};

// [[Rcpp::export]]
List calculate_window_dissimilarity_revised_cpp_parallel(
    List seq_list_r,
    int k,
    int window_size,
    int exclusion_half_width_bases
) {
    int n = seq_list_r.size();
    if (n == 0) {
        Rcpp::warning("Input sequence list is empty.");
        return List::create();
    }

    // Convert R list to C++ char vectors; check equal lengths
    std::vector<std::vector<char>> seq_chars_cpp(n);
    int len = 0;

    CharacterVector first_seq_r = as<CharacterVector>(seq_list_r[0]);
    if (first_seq_r.size() != 1) Rcpp::stop("Each element in seq_list_r must be a character vector of length 1.");
    std::string first_seq_str = as<std::string>(first_seq_r[0]);
    len = (int)first_seq_str.length();
    if (len == 0) Rcpp::stop("Sequence 1 is empty.");
    seq_chars_cpp[0].assign(first_seq_str.begin(), first_seq_str.end());

    for (int i = 1; i < n; ++i) {
        CharacterVector current_seq_r = as<CharacterVector>(seq_list_r[i]);
        if (current_seq_r.size() != 1) Rcpp::stop("Each element in seq_list_r must be a character vector of length 1.");
        std::string current_seq_str = as<std::string>(current_seq_r[0]);
        if ((int)current_seq_str.length() != len) {
            Rcpp::stop("Sequences must all have the same length. Sequence %d has length %d, expected %d.",
                       i + 1, (int)current_seq_str.length(), len);
        }
        seq_chars_cpp[i].assign(current_seq_str.begin(), current_seq_str.end());
    }

    // Exclusion window (spliced view)
    int E0, E1;
    compute_exclusion_base_range(len, exclusion_half_width_bases, E0, E1);

    // Outputs
    NumericMatrix basic_diss_r(n, n);
    NumericMatrix window_diss_r(n, n);
    IntegerMatrix min_positions_r(n, n);

    for (int i = 0; i < n; ++i){
        basic_diss_r(i,i)  = 0.0;
        window_diss_r(i,i) = 0.0;
        min_positions_r(i,i) = 1; // arbitrary diagonal
    }

    // Run
    DissimilarityWorkerRevised worker(
        seq_chars_cpp, k, window_size, len, n, E0, E1,
        basic_diss_r, window_diss_r, min_positions_r
    );
    parallelFor(0, n, worker);

    return List::create(
        _["basic_dissimilarity_matrix"]  = basic_diss_r,
        _["window_dissimilarity_matrix"] = window_diss_r,
        _["min_positions_matrix"]        = min_positions_r,
        _["excluded_base_start_1b"]      = E0 + 1,
        _["excluded_base_end_1b"]        = E1 + 1
    );
}
