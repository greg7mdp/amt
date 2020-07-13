// compile with: CL  -IC:/greg/github/cpp-btree -O2 lower_bound.cpp

#include <cstdio>
#include <ctime>
#include <vector>
#include <random>

#include <btree_map.h>
#include "../amt/amt.h"


using std::vector;

// --------------------------------------------------------------------------
template <class T = uint32_t>
T myrand()
{
    return ((T)rand()  << 13) ^ (T)rand();
}

// --------------------------------------------------------------------------
template <class T> 
void _fill(vector<T> &v)
{
    srand(1);   // for a fair/deterministic comparison 
    for (size_t i = 0, sz = v.size(); i < sz; ++i) 
        v[i] = myrand<T>();
}

// --------------------------------------------------------------------------
template <class T> 
void _shuffle(vector<T> &v)
{
    for (size_t n = v.size(); n >= 2; --n)
        std::swap(v[n - 1], v[static_cast<unsigned>(myrand<T>()) % n]);
}

// --------------------------------------------------------------------------
template <class M>
void pr(const M &m)
{
    printf("[");
    for (const auto& x : m)
        printf("%d, ", x.first);
    printf("]\n");
}

// --------------------------------------------------------------------------
template <class M1, class M2>
void test(const M1 &m1, const M2 &m2)
{
    // m2 is the reference map
    if (m2.empty())
        return;

    auto it = m2.end();
    --it;
    auto max = it->first;

    for (uint32_t i=0; i<max; ++i)
    {
        bool lb = m1.lower_bound(i)->first == m2.lower_bound(i)->first;
        bool ub = m1.upper_bound(i)->first == m2.upper_bound(i)->first;
        if (!(lb && ub))
            pr(m1);
        assert(lb && ub);
    }
}

// --------------------------------------------------------------------------
int main()
{
    using M1 = amt::amt<uint32_t, uint32_t>;
    using M2 = btree::btree_map<uint32_t, uint32_t>;

    M1 m1;
    M2 m2;

    for (uint32_t j=2; j<100; ++j)
    {
        m1.clear();
        m2.clear();

        for (uint32_t i=1; i<j; ++i) 
        {
            uint32_t x = myrand<uint32_t>() % (j * 10);
            m1.insert(typename M1::value_type(x, x));
            m2.insert(typename M2::value_type(x, x));
        }
        test(m1, m2);
    }
};

