#pragma once


#include "decomposition-types.h"


#if defined(NDEBUG)
#define DEBUG_LOG_NO_BREAK(...)
#define DEBUG_LOG(...)
#define DEBUG_ONLY(...)
#else
#include <iostream>
#define DEBUG_LOG_NO_BREAK(...) std::cout << __VA_ARGS__;
#define DEBUG_LOG(...) std::cout << __VA_ARGS__ << std::endl;
#define DEBUG_ONLY(...) __VA_ARGS__
#endif


#if defined(DECOMPOSITION_TEST)
#define DECOMPOSITION_API static
#else
#define DECOMPOSITION_API
#endif