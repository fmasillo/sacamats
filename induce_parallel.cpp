#undef max

#include <iostream>
#include <vector>
#include <string>
#include <array>
#include <limits>
#include <numeric>
#include <iomanip>
#include <future>
#include <ranges>
#include <execution>
#include <omp.h>
#include <thread>
#include <queue>
#include <functional>
//#include <boost/core/noinit_adaptor.hpp>

namespace psais {
#define L_TYPE 0
#define S_TYPE 1
#define NUM_THREADS 16u
//#define INDUCE_NUM_THREADS 16u

static uint32_t INDUCE_NUM_THREADS = 16u;

constexpr auto BLOCK_SIZE = 1u << 16;

constexpr inline auto mask = std::array<uint8_t, 8>{0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01};

template <typename T>
using NoInitVector = std::vector<T, std::allocator<T>>;

template <typename T>
constexpr auto EMPTY = std::numeric_limits<T>::max();

template<typename T>
T&& wrapper(std::remove_reference_t<T>&& t) { return std::move(t); }

static uint64_t counterL = 0;
static uint64_t counterS = 0;

template <typename Func, typename ... Args>
void parallel_do (
	std::integral auto n_jobs,
	std::integral auto n_threads,
	Func&& func,
	Args&& ...args
) {
	std::vector<std::jthread> threads;
	threads.reserve(n_jobs);
	auto counts = n_jobs / n_threads;
	auto remain = n_jobs % n_threads;

	for (decltype(n_jobs) tid = 0, L = 0; tid < n_threads; tid++) {
		auto block_size = counts + (tid < remain);
		if (block_size == 0) break;

		auto R = std::min(n_jobs, L + block_size);
		threads.emplace_back(std::forward<Func>(func), L, R, tid, wrapper<Args>(args)...);
		L = R;
	}
}

// #parallel_init
template <std::ranges::input_range R>
void parallel_init(R&& r, const auto& value) {
	std::for_each(std::execution::par, r.begin(), r.end(),
		[&value](auto &x) { x = value; }
	);
}


class ThreadPool {
  public:
	ThreadPool(int cnt) : thread_cnt(cnt) {
		auto init_thread = [this](std::stop_token tok) {
			auto job = std::packaged_task<void()>{};

			while (!tok.stop_requested()) {
				if (auto lock = std::scoped_lock(q_m); jobs.empty()) {
					continue;
				} else {
					job = std::move(jobs.front());
					jobs.pop();
				}
				job();
			}
		};

		for (int i = 0; i < thread_cnt; i++) {
			threads.emplace_back(init_thread);
		}
	}

	template <typename F, typename ... Args>
	requires std::invocable<F, Args...>
	auto enqueue(F&& f, Args&& ... args) {
		using res_type = typename std::invoke_result_t<F, Args...>;

		auto task = std::packaged_task<res_type()>(
			std::bind(std::forward<F>(f), std::forward<Args>(args)...)
		);
		auto res = task.get_future();

		auto lock = std::unique_lock(q_m);
		jobs.emplace(std::move(task));
		return res;
	}

  private:
	std::queue<std::packaged_task<void()>> jobs;
	std::vector<std::jthread> threads;
	std::mutex q_m;

	int thread_cnt;
};

// ##prepare
template<typename IndexType, typename IndexSmallType, typename CharType, typename ElementType>
void prepare_no_type(
	const size_t L,
	const std::ranges::random_access_range auto &S,
	const NoInitVector<ElementType> &SA,
    const NoInitVector<IndexType> &docBoundaries,
	NoInitVector<std::pair<CharType, uint8_t>> &RB
) {
	if (L >= SA.size()) return;
	decltype(L) R = std::min(SA.size(), L + BLOCK_SIZE);

	#pragma omp parallel for num_threads(INDUCE_NUM_THREADS)
	for (auto i = L; i < R; i++) {
		//IndexType induced_idx = SA[i].first + docBoundaries[SA[i].second-1] - 1 ;
		if (SA[i].first == 0 or 0 > SA[i].second) {
			RB[i - L] = {EMPTY<CharType>, 0};
		} else {
			RB[i - L] = {S[SA[i].first + docBoundaries[SA[i].second-1] - 1], S[SA[i].first + docBoundaries[SA[i].second-1] - 2]};
		}
	}
}

// ##update !!!NOT USED!!!
template<typename IndexType, typename IndexSmallType, typename ElementType>
void update(
	const size_t L,
	const NoInitVector<std::pair<IndexType, ElementType>> &WB,
	NoInitVector<ElementType> &SA
) {
	if (L >= SA.size()) return;
	decltype(L) R = std::min(SA.size(), L + BLOCK_SIZE);

	#pragma omp parallel for num_threads(INDUCE_NUM_THREADS / 2)
	for (auto i = L; i < R; i++) {
		auto& [idx, val] = WB[i - L];
		if (idx != EMPTY<IndexType>) {
			SA[idx] = val;
		}
	}
}

// ##induce_impl
template<auto InduceType, typename IndexType, typename IndexSmallType, typename CharType, typename ElementType>
void induce_impl_no_type (
	const std::ranges::random_access_range auto &S,
	const std::ranges::input_range auto &rng,
	const IndexType L,
	NoInitVector<ElementType> &SA,
    NoInitVector<IndexType> &docBoundaries,
	NoInitVector<std::pair<CharType, uint8_t>> &RB,
	//NoInitVector<std::pair<IndexType, ElementType>> &WB,
	NoInitVector<IndexType> &ptr
) {
	for (IndexType i : rng) {
		if(SA[i].first > 0){
			ElementType j = SA[i];
			if constexpr (InduceType == L_TYPE) SA[i].second = ~SA[i].second;

			if(0 < j.second){
				auto [c0, c1] = RB[i - L];
				//CharType c0 = EMPTY<CharType>;
				//CharType c1 = EMPTY<CharType>;
				CharType chr;
				if constexpr (InduceType == L_TYPE){
					if(c0 == EMPTY<CharType> && S[j.first + docBoundaries[j.second-1] - 1] >= S[j.first + docBoundaries[j.second-1]]){
						c0 = S[j.first + docBoundaries[j.second-1] - 1];
						c1 = S[j.first + docBoundaries[j.second-1] - 2];
					}
				}
				else{
					if(c0 == EMPTY<CharType> && S[j.first + docBoundaries[j.second-1] - 1] <= S[j.first + docBoundaries[j.second-1]]){
						//assert(S[j.first + docBoundaries[j.second-1] - 1] <= S[j.first + docBoundaries[j.second-1]]);
						//if(ptr[S[j.first + docBoundaries[j.second-1] - 1]] < i){
							c0 = S[j.first + docBoundaries[j.second-1] - 1];
							c1 = S[j.first + docBoundaries[j.second-1] - 2];
						//}
					}
				}
				chr = c0;

				ElementType tuple;
				if constexpr (InduceType == L_TYPE) {
					tuple = ((j.first - 1 > 0) && (c1 < c0)) ? ElementType(j.first - 1, ~j.second) : ElementType(j.first - 1, j.second);
				} else {
					tuple = ((j.first - 1 == 0) || (c1 > c0)) ? ElementType(j.first - 1, ~j.second) : ElementType(j.first - 1, j.second);
				} 

				//bool is_adjacent;
				auto pos = ptr[chr];
				if constexpr (InduceType == L_TYPE) {
					ptr[chr] += 1;
					//is_adjacent = pos < L + (BLOCK_SIZE << 1);
					//is_adjacent = pos >= L & pos < L + BLOCK_SIZE;
				} else {
					ptr[chr] -= 1;
					//is_adjacent = pos + BLOCK_SIZE >= L;
					//is_adjacent = pos >= L & pos < L + BLOCK_SIZE;
				}

				// if pos is in adjacent block -> directly write it
				// otherwise, write it to WB
				//if (is_adjacent) {
					SA[pos] = tuple;//induced_idx;
					//WB[i - L].first = EMPTY<IndexType>;
					//if (tuple.second > 0) {
						//if (pos >= L & pos < L + BLOCK_SIZE) {
					//		RB[pos - L] = {S[tuple.first + docBoundaries[tuple.second-1] - 1], S[tuple.first + docBoundaries[tuple.second-1] - 2]};
						//}
					//}
				//} else {
				//	WB[i - L] = {pos, tuple};
				//}
			}
			else if constexpr (InduceType == S_TYPE) SA[i].second = ~SA[i].second;
		}
		else if constexpr (InduceType == S_TYPE){
			if (0 > SA[i].second) SA[i].second = ~SA[i].second;
		} 		
	}
}

// ##induce
template<auto InduceType, typename IndexType, typename IndexSmallType, typename CharType, typename ElementType>
void induce (
	const std::ranges::random_access_range auto &S,
	//const TypeVector &T,
	NoInitVector<ElementType> &SA,
    NoInitVector<IndexType> &docBoundaries,
	NoInitVector<std::pair<CharType, uint8_t>> &RBP,
	NoInitVector<std::pair<CharType, uint8_t>> &RBI,
	//NoInitVector<std::pair<IndexType, ElementType>> &WBU,
	//NoInitVector<std::pair<IndexType, ElementType>> &WBI,
	NoInitVector<IndexType> &ptr
) {
	// views
	constexpr auto iter_view = [] {
		if constexpr (InduceType == L_TYPE) {
			return std::views::all;
		} else {
			return std::views::reverse;
		}
	}();

	IndexType size = SA.size();
	auto blocks = std::views::iota(IndexType(0), size)
		| std::views::filter([](IndexType n) { return n % BLOCK_SIZE == 0; });

	// prepare for first block
	if constexpr (InduceType == L_TYPE) {
		//prepare<IndexType, IndexSmallType, CharType, ElementType, decltype(S)>(0, S, SA, docBoundaries, T, RBP);
		prepare_no_type<IndexType, IndexSmallType, CharType, ElementType, decltype(S)>(0, S, SA, docBoundaries, RBP);
	} else {
		//prepare<IndexType, IndexSmallType, CharType, ElementType, decltype(S)>(size / BLOCK_SIZE * BLOCK_SIZE, S, SA, docBoundaries, T, RBP);
		prepare_no_type<IndexType, IndexSmallType, CharType, ElementType, decltype(S)>(size / BLOCK_SIZE * BLOCK_SIZE, S, SA, docBoundaries, RBP);
	}

	auto pool = ThreadPool(2);
	auto stage = std::array<std::future<void>, 2>{};

	// pipeline
	for (IndexType L : blocks | iter_view) {
		for (auto &s : stage) if (s.valid()) s.wait();
		RBI.swap(RBP);
		//WBI.swap(WBU);

		// prepare && update
		IndexType P_L = L + BLOCK_SIZE;
		IndexType U_L = L - BLOCK_SIZE;
		if constexpr (InduceType == S_TYPE) {
			std::swap(P_L, U_L);
		}

		//stage[0] = pool.enqueue(prepare<IndexType, IndexSmallType, CharType, ElementType, decltype(S)>, P_L,
		//		std::ref(S), std::ref(SA), std::ref(docBoundaries), std::ref(T), std::ref(RBP));
		stage[0] = pool.enqueue(prepare_no_type<IndexType, IndexSmallType, CharType, ElementType, decltype(S)>, P_L,
				std::ref(S), std::ref(SA), std::ref(docBoundaries), std::ref(RBP));
		//stage[1] = pool.enqueue(update<IndexType, IndexSmallType, ElementType>, U_L,
		//		std::ref(WBU), std::ref(SA));

		// induce
		auto rng = std::views::iota(L, std::min(L + BLOCK_SIZE, size)) | iter_view;
		//induce_impl<InduceType, IndexType, IndexSmallType, CharType, ElementType>(S, T, rng, L, SA, docBoundaries, RBI, WBI, ptr);
		induce_impl_no_type<InduceType, IndexType, IndexSmallType, CharType, ElementType>(S, rng, L, SA, docBoundaries, RBI, ptr);
	}
}

// ##induce_sort
// S : input string
// T : type vector
// SA1 : SA of recursion
// LMS : LMS positions
// BA : bucket array
// SA : output SA
template<typename IndexType, typename IndexSmallType, typename ElementType>
void induce_sort(
	//const std::ranges::random_access_range auto &S,
	const std::string &S,
	//const TypeVector &T,
	NoInitVector<IndexType> &BA,
	NoInitVector<ElementType> &SA,
    NoInitVector<IndexType> &docBoundaries,
	uint32_t n_threads
) {
	INDUCE_NUM_THREADS = n_threads;
	using CharType = decltype(S.begin())::value_type;

	// induce LMS
	//put_lms(S, LMS, SA1, BA, SA);

	// declare ptr, RBP, RBI, WBI, WBU. 
    // RBP->Read Buffer Prepare, RBI->Read Buffer Induce, WBU->Write Buffer Update, WBI->Write Buffer Induce
	NoInitVector<IndexType> ptr(BA.size());
	NoInitVector<std::pair<CharType, uint8_t>> RBP(BLOCK_SIZE), RBI(BLOCK_SIZE);
	//NoInitVector<std::pair<IndexType, ElementType>> WBU(BLOCK_SIZE), WBI(BLOCK_SIZE);

	// init buffer
    parallel_init(RBP, std::pair{EMPTY<CharType>, uint8_t(0)});
	parallel_init(RBI, std::pair{EMPTY<CharType>, uint8_t(0)});
	//parallel_init(WBU, std::pair{EMPTY<IndexType>, EMPTY<ElementType>});
	//parallel_init(WBI, std::pair{EMPTY<IndexType>, EMPTY<ElementType>});

	// induce L
	//ptr[0] = 0;
	std::transform(std::execution::par_unseq, BA.begin(), BA.end() - 1, ptr.begin(),
		[](IndexType &b) { return b; }
	);
	//induce<L_TYPE, IndexType, IndexSmallType, CharType, ElementType>(S, T, SA, docBoundaries, RBP, RBI, WBU, WBI, ptr);
	induce<L_TYPE, IndexType, IndexSmallType, CharType, ElementType>(S, SA, docBoundaries, RBP, RBI, ptr);

    std::cerr << "induce L done" << std::endl;
	//std::cerr << "counterL = " << counterL << std::endl;
    //for(IndexType i = 0; i < SA.size(); i++){
    //    std::cout << "SA[" << i << "] = " << SA[i].first << ", " << SA[i].second << std::endl;
    //}
	// init buffer
	parallel_init(RBP, std::pair{EMPTY<CharType>, uint8_t(0)});
	parallel_init(RBI, std::pair{EMPTY<CharType>, uint8_t(0)});
	//parallel_init(WBU, std::pair{EMPTY<IndexType>, EMPTY<ElementType>});
	//parallel_init(WBI, std::pair{EMPTY<IndexType>, EMPTY<ElementType>});

	// clean S_TYPE
	//#pragma omp parallel for
	//for(IndexType i = docBoundaries.size(); i < SA.size(); i++){
	//	if ((IndexType)SA[i].first + docBoundaries[SA[i].second - 1] != EMPTY<IndexType> and T.get((IndexType)SA[i].first + docBoundaries[SA[i].second - 1]) == S_TYPE) {
	//		SA[i].first = EMPTY<IndexSmallType>;
	//	}
	//}
	//std::for_each(std::execution::par_unseq, SA.begin() + docBoundaries.size(), SA.end(),
	//	[&T, &docBoundaries](ElementType &idx) {
	//		if ((IndexType)idx.first + docBoundaries[idx.second - 1] != EMPTY<IndexType> and T.get((IndexType)idx.first + docBoundaries[idx.second - 1]) == S_TYPE) {
	//			idx.first = EMPTY<IndexSmallType>;
	//		}
	//	}
	//);
    //std::cerr << "cleaned s-type" << std::endl;
	//for(IndexType i = 0; i < SA.size(); i++){
    //    std::cout << "SA[" << i << "] = " << SA[i].first << ", " << SA[i].second << std::endl;
    //}
	
	for(size_t i = 0; i < 256; i++){
		ptr[i] = BA[i+1]-1;
	}
	// induce S
	//std::transform(std::execution::par_unseq, BA.begin(), BA.end(), ptr.begin(),
	//	[](auto &bucket) { return bucket - 1; }
	//);
    std::cerr << "reset BA" << std::endl;
	//induce<S_TYPE, IndexType, IndexSmallType, CharType, ElementType>(S, T, SA, docBoundaries, RBP, RBI, WBU, WBI, ptr);
	induce<S_TYPE, IndexType, IndexSmallType, CharType, ElementType>(S, SA, docBoundaries, RBP, RBI, ptr);
    std::cerr << "induce S done" << std::endl;
	std::cerr << "counterS = " << counterS << std::endl;
	std::cerr << "counterTotal = " << counterL + counterS << std::endl;
    //for(IndexType i = 0; i < SA.size(); i++){
    //    std::cout << "SA[" << i << "] = " << SA[i].first << ", " << SA[i].second << std::endl;
    //}
}

}
