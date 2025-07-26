
#include <cstdlib>
#include <string>

// Sets the number of threads for Accelerate via the environment variable VECLIB_MAXIMUM_THREADS.
void accelerate_set_num_threads(int num_threads) {
    // Respect user-defined number of threads for Accelerate and set if unset
    const char* env_p = std::getenv("VECLIB_MAXIMUM_THREADS");
    if (!env_p)
        setenv("VECLIB_MAXIMUM_THREADS", std::to_string(num_threads).c_str(), 1);
}