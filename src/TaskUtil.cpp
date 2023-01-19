#include "TaskUtil.h"

uint64_t intRand(const uint64_t & min, const uint64_t & max) {
    static thread_local std::mt19937 generator;
    std::uniform_int_distribution<uint64_t> distribution(min,max);
    return (uint64_t)distribution(generator);
}