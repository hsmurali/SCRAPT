#ifndef UTILITY_HPP
#define UTILITY_HPP
#include <vector>
template<typename ElementType>

inline const ElementType *getPointerOfVector(const std::vector<ElementType> &v)
{
    return v.empty() ? 0 : &v[0];
}


template<typename ElementType>
inline ElementType *getPointerOfVector(std::vector<ElementType> &v)
{
    return v.empty() ? 0 : &v[0];
}

#endif // UTILITY_HPP
