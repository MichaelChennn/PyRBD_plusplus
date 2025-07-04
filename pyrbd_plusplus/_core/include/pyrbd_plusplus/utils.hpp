#pragma once
#include <pyrbd_plusplus/common.hpp>
#include <algorithm>

namespace pyrbdpp::utils
{   
    /**
     * @brief Check if set1 is a subset of set2.
     * @param set1 The potential subset.
     * @param set2 The set to check against.
     * @return True if set1 is a subset of set2, false otherwise.
     */
    bool isSubSet(const Set &set1, const Set &set2);
}