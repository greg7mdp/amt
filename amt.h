#if !defined(amt_h_guard_)
#define amt_h_guard_

// ---------------------------------------------------------------------------
// Copyright (c) 2020, Gregory Popovitch - greg7mdp@gmail.com
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// ---------------------------------------------------------------------------

#ifdef _MSC_VER
    #pragma warning(push)  

    #pragma warning(disable : 4127) // conditional expression is constant
    #pragma warning(disable : 4324) // structure was padded due to alignment specifier
    #pragma warning(disable : 4514) // unreferenced inline function has been removed
    #pragma warning(disable : 4623) // default constructor was implicitly defined as deleted
    #pragma warning(disable : 4625) // copy constructor was implicitly defined as deleted
    #pragma warning(disable : 4626) // assignment operator was implicitly defined as deleted
    #pragma warning(disable : 4710) // function not inlined
    #pragma warning(disable : 4711) // selected for automatic inline expansion
    #pragma warning(disable : 4820) // '6' bytes padding added after data member
    #pragma warning(disable : 4868) // compiler may not enforce left-to-right evaluation order in braced initializer list
    #pragma warning(disable : 5027) // move assignment operator was implicitly defined as deleted
    #pragma warning(disable : 5045) // Compiler will insert Spectre mitigation for memory load if /Qspectre switch specified
#endif



#include <algorithm>
#include <cmath>
#include <cstring>
#include <iterator>
#include <limits>
#include <memory>
#include <array>
#include <utility>
#include <type_traits>
#include <cassert>

namespace amt_   // internal namespace
{

// ---------------------------------------------------------------------------
// popcount stuff
// ---------------------------------------------------------------------------
static inline uint32_t s_amt_popcount_default(uint32_t i) noexcept
{
    i = i - ((i >> 1) & 0x55555555);
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
    return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

static inline uint32_t s_amt_popcount_default(uint64_t x) noexcept
{
    const uint64_t m1  = uint64_t(0x5555555555555555); // binary: 0101...
    const uint64_t m2  = uint64_t(0x3333333333333333); // binary: 00110011..
    const uint64_t m4  = uint64_t(0x0f0f0f0f0f0f0f0f); // binary:  4 zeros,  4 ones ...
    const uint64_t h01 = uint64_t(0x0101010101010101); // the sum of 256 to the power of 0,1,2,3...

    x -= (x >> 1) & m1;             // put count of each 2 bits into those 2 bits
    x = (x & m2) + ((x >> 2) & m2); // put count of each 4 bits into those 4 bits 
    x = (x + (x >> 4)) & m4;        // put count of each 8 bits into those 8 bits 
    return (x * h01)>>56;           // returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24)+...
}

#if defined __clang__

    #if defined(i386)
        #include <cpuid.h>
        inline void amt_cpuid(int info[4], int InfoType) {
            __cpuid_count(InfoType, 0, info[0], info[1], info[2], info[3]);
        }
    #endif

    #define AMT_POPCNT   __builtin_popcount
    #define AMT_POPCNT64 __builtin_popcountll

#elif defined __GNUC__

    #if defined(i386)
        #include <cpuid.h>
        inline void amt_cpuid(int info[4], int InfoType) {
            __cpuid_count(InfoType, 0, info[0], info[1], info[2], info[3]);
        }
    #endif

    // __POPCNT__ defined when the compiled with popcount support
    // (-mpopcnt compiler option is given for example)
    #ifdef __POPCNT__
        // slower unless compiled iwith -mpopcnt
        #define AMT_POPCNT   __builtin_popcount
        #define AMT_POPCNT64 __builtin_popcountll
    #endif


#elif defined _MSC_VER

    #include <intrin.h>                     // for __popcnt()

    #define AMT_POPCNT_CHECK  // slower when defined, but we have to check!
    #define amt_cpuid(info, x)    __cpuid(info, x)

    #define AMT_POPCNT __popcnt
    #if (INTPTR_MAX == INT64_MAX)
        #define AMT_POPCNT64 __popcnt64
    #endif

#endif

#if defined(AMT_POPCNT_CHECK)
static inline bool amt_popcount_check()
{
    int cpuInfo[4] = { -1 };
    amt_cpuid(cpuInfo, 1);
    if (cpuInfo[2] & (1 << 23))
        return true;   // means AMT_POPCNT supported
    return false;
}
#endif

#if defined(AMT_POPCNT_CHECK) && defined(AMT_POPCNT)

static inline uint32_t amt_popcount(uint32_t i)
{
    static const bool s_ok = amt_popcount_check();
    return s_ok ? AMT_POPCNT(i) : s_amt_popcount_default(i);
}

#else

static inline uint32_t amt_popcount(uint32_t i)
{
#if defined(AMT_POPCNT)
    return static_cast<uint32_t>(AMT_POPCNT(i));
#else
    return s_amt_popcount_default(i);
#endif
}

#endif

// --------------------------------------------------------------------------
// for leaves, SV and V are the same
// for non-leaves, SV is a sparsegroup *, and V is the actual amt mapped_type
// --------------------------------------------------------------------------
template <class K, class SV, class V, uint32_t N>  
class sparsegroup
{
    using bm_type    = uint16_t;
    using size_type  = uint8_t;
    using group      = sparsegroup<K, SV, V, N>;
    using group_ptr  = group *;
    using leaf_group = sparsegroup<K, V, V, N>;
    using leaf_group_ptr = SV;

public:
    struct locator
    {
        locator(sparsegroup *g, uint32_t idx) :
             _group(g), _idx(idx)
        {}

        sparsegroup   *_group;
        uint32_t       _idx;
    };

    struct insert_locator : public locator
    {
        insert_locator(sparsegroup *g, uint32_t idx, bool created) :
            locator(g, idx), _created(created)
        {}

        bool _created;
    };

    sparsegroup(sparsegroup *parent, uint32_t depth, K partial_key) :
        _bitmap(0),
        _num_val(0),
        _num_alloc(N),
        _depth(depth),
        _partial_key(partial_key),
        _parent(parent)
    {
    }

    K        partial_key() const { return _partial_key; }
    bool     is_leaf() const     { return _depth == max_depth; }
    uint32_t depth()             { return _depth; }
    

    void shrink()
    {
        assert(is_leaf() && _parent);
        if (_num_alloc - _num_val <= 1)
            return;
        resize(_num_val);
    }

    void swap(sparsegroup &o)
    {
        using std::swap;
        for (uint32_t i=0; i<_num_val; ++i)
           swap(_values[i], o._values[i]);

        // swap _num_val
        uint32_t tmp = _num_val;
        _num_val = o._num_val;
        o._num_val = tmp;
    }

    void resize(uint32_t sz)
    {
        if (_depth == 0)
            return; // don't resize root sparsegroup as we can't update the parent
        assert(sz >= _num_val);
        group_ptr grp = allocate_group(sz, _parent, _depth, _partial_key);
        grp->_bitmap = _bitmap;

        if (is_leaf())
            reinterpret_cast<leaf_group_ptr>(this)->swap(*reinterpret_cast<leaf_group_ptr>(grp));
        else
            this->swap(*grp);

        size_type parent_idx = _parent->_nibble(_partial_key);
        assert((group_ptr)_parent->get_value(parent_idx) == this);
        _parent->set_value(parent_idx, (SV)grp);
        
        deallocate_group();
    }

    locator begin() const
    {
        assert(0); // todo
        return { nullptr, 0 };
    }
        
    locator find(K key) const
    {
        size_type n = _nibble(key);

        if (_bmtest(n)) 
        {
            // found it
            size_type offset = _pos_to_offset(n);
            if (_depth == max_depth)
                return { this, offset };
            else
                return reinterpret_cast<group_ptr>(_values[offset])->find(key);
        }
        else
            return { nullptr, 0 };
    }

    insert_locator find_or_prepare_insert(K key)
    {
        size_type n = _nibble(key);

        if (_bmtest(n)) // found it
        {
            size_type offset = _pos_to_offset(n);
            if (_depth == max_depth)
                return insert_locator(this, offset, false);
            else
                return reinterpret_cast<group_ptr>(_values[offset])->find_or_prepare_insert(key);
        }
        else
        {
            // need to insert
            
            
        }
            
        
        
        return { nullptr, 0, false };
    }

    bool match(K key) const 
    {
        assert(_depth == max_depth); // we are matching only leaf groups
        return (key & ~bm_mask) == _partial_key;
    }

    template <class Val>
    void set_value(uint32_t idx, Val &&v)
    {
        _values[idx] = std::forward<Val>(v);
    }

    SV& get_value(uint32_t idx)
    {
        return _values[idx];
    }
    
    static group_ptr allocate_group(uint32_t cnt, group_ptr parent, uint32_t depth, K partial_key)
    {
        bool leaf = (depth == max_depth);
        if (!leaf)
        {
            group_ptr grp = (group_ptr)malloc(group_size_in_bytes<SV>(cnt));
            new (grp) group(parent, depth, partial_key);
            for (uint32_t i=0; i<cnt; ++i)
                grp->_values[i] = nullptr;
            return grp;
        }
        else
        {
            leaf_group_ptr grp = (leaf_group_ptr)malloc(group_size_in_bytes<V>(cnt));
            new (grp) leaf_group((leaf_group_ptr)parent, depth, partial_key);
            grp->construct_values();
            return (group_ptr)grp;
        }
    }

    void construct_values()
    {
        for (uint32_t i=0; i<_num_val; ++i)
            new (&_values[i]) SV();
    }

    // deallocates the group and all its children
    void deallocate_group()
    {
        bool leaf = (_depth == max_depth);
        if (leaf)
            reinterpret_cast<leaf_group_ptr>(this)->deallocate_leaf();
        else
        {
            for (uint32_t i=0; i<_num_val; ++i)
                reinterpret_cast<group_ptr>(_values[i])->deallocate_group();
            free(this);
        }
    }

    void deallocate_leaf()
    {
        for (uint32_t i=0; i<_num_val; ++i)
            _values[i].~V();
        free(this);
    }
    
private:

    bool _bmtest(size_type i) const   { return !!(_bitmap & (static_cast<uint32_t>(1) << i)); }
    void _bmset(size_type i)          { _bitmap |= static_cast<uint32_t>(1) << i; }
    void _bmclear(size_type i)        { _bitmap &= ~(static_cast<uint32_t>(1) << i); }

    template <class Val>
    static size_t group_size_in_bytes(uint32_t cnt)   
    {
        assert(cnt >= 1);
        size_t x = sizeof(group) + (cnt - 1) * sizeof(Val);
        return (x + 0xF) & ~0xF;
    }

    size_type _pos_to_offset(size_type pos)
    {
        return static_cast<size_type>(amt_popcount(_bitmap & ((static_cast<uint32_t>(1) << pos) - 1)));
    }

    size_type _offset_to_pos(size_type offset)
    {
        uint32_t bm = _bitmap;

        for (; offset > 0; offset--)
            bm &= (bm-1);  // remove right-most set bit

        // Clear all bits to the left of the rightmost bit (the &),
        // and then clear the rightmost bit but set all bits to the
        // right of it (the -1).
        // --------------------------------------------------------
        bm = (bm & -bm) - 1;
        return  static_cast<size_type>(amt_popcount(bm));
    }

    size_type _nibble(K key) 
    {
        return (size_type)((key >> ((max_depth - _depth) * bm_shift)) & bm_mask);
    }

    static constexpr const uint32_t bm_shift = 4;
    static constexpr const uint32_t bm_mask  = 0xF;
    static constexpr const uint32_t max_depth = (sizeof(K) * 8) / bm_shift - 1;
    
    uint32_t     _bitmap : 16;
    uint32_t     _num_val : 5;
    uint32_t     _num_alloc : 5;
    uint32_t     _depth : 6;
    K            _partial_key;  // only highest bits 0 -> 4 * depth are set
    group_ptr    _parent;
    SV           _values[N];
};


}; // amt_

namespace amt
{
using namespace amt_;

// -----------------------------------------------------------------------------
// An Array Mapped Tree (amt) is an ordered associative container which mapping 
// unsigned integer keys to values. It is optimized for both speed and memory 
// footprint in most common use cases. 
// see paper: "Fast And Space Efficient Trie Searches" by Phil Bagwell
//     http://lampwww.epfl.ch/papers/triesearches.pdf.gz
//
// Its interface is similar to that of `std::unordered_map<K, V>`.
// -----------------------------------------------------------------------------
template <class K, class V> 
class amt 
{
    using leaf_group     = sparsegroup<K, V, V, 1>;
    using leaf_group_ptr = leaf_group *;

    using group          = sparsegroup<K, leaf_group_ptr, V, 1>;
    using group_ptr      = group *;

    using locator        = typename group::locator;
    using insert_locator = typename group::insert_locator;

public:

    static_assert(std::is_integral<K>::value && std::is_unsigned<K>::value &&
                  (sizeof(K) == 2 || sizeof(K) == 4 || sizeof(K) == 8), "supporting 16, 32 or 64 bit unsigned integer keys.");

    static bool _equal(const K &a, const K &b) { return a == b; }
    
    using key_type        = K;
    using value_type      = std::pair<const K, const V>;
    using mapped_type     = V;
    using size_type       = std::size_t;
    using difference_type = std::ptrdiff_t;
    using key_equal       = std::equal_to<K>;


    // --------------------- iterator ---------------------------------------
    class iterator 
    {
        friend class amt;
        using group_ptr         = typename amt::group_ptr;
        using leaf_group_ptr    = typename amt::leaf_group_ptr;
        using locator           = typename amt::locator;

    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type        = typename amt::value_type;
        using reference         = value_type&;
        using const_reference   = const value_type&;
        using pointer           = value_type *;
        using const_pointer     = const value_type *;
        using difference_type   = typename amt::difference_type;

        iterator(group_ptr group = nullptr, uint32_t idx = 0) : 
            _group(group), _idx(idx) 
        {
            assert(!_group || _group->is_leaf());
        }

        reference operator*() const { return _v; }

        pointer operator->() const { return &operator*(); }

        iterator& operator++() 
        {
            assert(0);
            return *this;
        }

        iterator operator++(int) 
        {
            auto tmp = *this;
            ++*this;
            return tmp;
        }

        friend bool operator==(const iterator& a, const iterator& b) 
        {
            return a._group == b._group && a._idx == b._idx;
        }

        friend bool operator!=(const iterator& a, const iterator& b) 
        {
            return !(a == b);
        }

    private:

        group_ptr _group;
        uint32_t  _idx;

        // anonymous union to avoid uninitialized member warnings
        union {
            value_type _v;
        };
    };

    // --------------- const iterator ---------------------------------------
    class const_iterator 
    {
        friend class amt;

    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type        = typename amt::value_type;
        using reference         = const value_type&;
        using pointer           = const value_type *;
        using difference_type   = typename amt::difference_type;

        const_iterator() {}
        const_iterator(iterator it) : _it(std::move(it)) {}

        reference operator*() const { return *_it; }
        pointer operator->() const  { return _it.operator->(); }

        const_iterator& operator++() 
        {
            ++_it;
            return *this;
        }
        const_iterator operator++(int) { return _it++; }

        friend bool operator==(const const_iterator& a, const const_iterator& b) 
        {
            return a._it == b._it;
        }

        friend bool operator!=(const const_iterator& a, const const_iterator& b) 
        {
            return !(a == b);
        }

    private:
        iterator _it;
    };
    
    // --------------------- constructors -----------------------------------
    amt() noexcept :
        _last_insert(nullptr),
        _last_lookup(nullptr),
        _root(nullptr),
        _size(0)
    {
        _init();
    }
    
    template <class InputIter>
    amt(InputIter first, InputIter last) : amt() 
    {
        insert(first, last);
    }

    amt(std::initializer_list<value_type> init) : amt(init.begin(), init.end())
    {}

    amt(const amt &o) : amt() 
    {
        assert(0); //todo
    }

    amt(amt &&o) noexcept : amt() 
    {
        swap(o);
    }

    amt& operator=(const amt &o)
    {
        amt tmp(o);
        swap(tmp);
        return *this;
    }

    amt& operator=(amt &&o)
    {
        amt tmp(std::move(o));
        swap(tmp);
        return *this;
    }

    // --------------------- desstructor -----------------------------------
    ~amt()
    {
        _cleanup();
    }

    // --------------------- apis ------------------------------------------
    iterator begin() 
    {
        auto loc = _root->begin();
        return iterator(loc._group, loc._idx);
    }
    
    iterator end() 
    {
        return iterator(nullptr, 0); 
    }

    const_iterator begin() const  { return const_cast<amt *>(this)->begin();  }
    const_iterator end() const    { return const_cast<amt *>(this)->end(); }
    const_iterator cbegin() const { return begin(); }
    const_iterator cend() const   { return end(); }

    bool empty() const            { return !size(); }
    size_type size() const        { return _size; }
    size_type capacity() const    { return max_size(); }
    size_type max_size() const    { return (std::numeric_limits<size_type>::max)(); }

    void clear()
    {
        _cleanup();
        _init();
    }

    // ------------------------------------ erase ------------------------------
    iterator erase(const_iterator cit)
    {
        return erase(cit._it_);
    }

    iterator erase(iterator it) 
    { 
        _erase(it++); 
        return it; 
    }

    size_type erase(const key_type& key)
    {
        auto it = find(key);
        if (it == end()) 
            return 0;
        erase(it);
        return 1;
    }

    iterator erase(const_iterator first, const_iterator last)
    {
        while (first != last)
            _erase(first++);
        return last._it;
    }
    
    // ------------------------------------ insert ------------------------------
    template <class VT>
    std::pair<iterator, bool> insert(VT&& value) 
    {
        auto loc = insert_impl(value.first, value.second);
        return { iterator(loc._group, loc._idx), loc._created };
    }

    template <class VT>
    iterator insert(const_iterator hint, VT&& value) 
    { 
        return insert(std::forward<VT>(value)).first; 
    }

    template <class InputIt>
    void insert(InputIt first, InputIt last) 
    {
        for (; first != last; ++first) 
            insert(*first);
    }

    void insert(std::initializer_list<value_type> ilist)
    {
        insert(ilist.begin(), ilist.end());
    }

    template <class... Args>
    std::pair<iterator, bool> emplace(Args&&... args) 
    {
        return insert(value_type(std::forward<Args>(args)...));
    }
    
    template <class... Args>
    iterator emplace_hint(const_iterator hint, Args&&... args) 
    {
        return insert(hint, value_type(std::forward<Args>(args)...));        
    }
    
    template <class Key, class... Args>
    std::pair<iterator, bool> try_emplace(Key&& key, Args&&... args) 
    {
        auto loc =  insert_impl(key, std::piecewise_construct, 
                                std::forward_as_tuple(std::forward<Key>(key)), 
                                std::forward_as_tuple(std::forward<Args>(args)...));
        return { iterator(loc._group, loc._idx), loc._created };
    }
    
    template <class Key, class... Args>
    iterator try_emplace(const_iterator hint, Key&& key, Args&&... args) 
    { 
        return try_emplace(std::forward<Key>(key), std::forward<Args>(args)...).first;
    }
    

    // ------------------------------------ insert_or_assign ------------------------------
    template <class M = mapped_type>
    std::pair<iterator, bool> insert_or_assign(key_type&& k, M&& v) 
    {
        return insert_or_assign_impl(std::forward<K>(k), std::forward<M>(v));
    }

    template <class M = mapped_type>
    std::pair<iterator, bool> insert_or_assign(key_type&& k, const M& v) 
    {
        return insert_or_assign_impl(std::forward<K>(k), v);
    }

    template <class M = mapped_type>
    std::pair<iterator, bool> insert_or_assign(const key_type& k, M&& v) 
    {
        return insert_or_assign_impl(k, std::forward<M>(v));
    }

    template <class M = mapped_type>
    std::pair<iterator, bool> insert_or_assign(const key_type& k, const M& v) 
    {
        return insert_or_assign_impl(k, v);
    }

    template <class M = mapped_type>
    iterator insert_or_assign(const_iterator, key_type&& k, M&& v) 
    {
        return insert_or_assign(std::forward<K>(k), std::forward<M>(v)).first;
    }

    template <class M = mapped_type>
    iterator insert_or_assign(const_iterator, key_type&& k, const M& v) 
    {
        return insert_or_assign(std::forward<K>(k), v).first;
    }

    template <class M = mapped_type>
    iterator insert_or_assign(const_iterator, const key_type& k, M&& v) 
    {
        return insert_or_assign(k, std::forward<M>(v)).first;
    }

    template <class M = mapped_type>
    iterator insert_or_assign(const_iterator, const key_type& k, const M& v) 
    {
        return insert_or_assign(k, v).first;
    }

#if 0
    extract;
    merge;
#endif

    void reserve(size_type count) {}
    void rehash(size_type count) {}

    void swap(amt& o) noexcept
    {
        std::swap(_root, o._root);
        std::swap(_last_insert, o._last_insert);
        std::swap(_last_lookup, o._last_lookup);
        std::swap(_size, o._size);
    }

    V& at(const key_type& key)
    {
        auto loc = this->find_impl(key);
        if (!loc._group) 
            throw std::out_of_range("amt at(): lookup non-existent key");
        return reinterpret_cast<leaf_group_ptr>(loc._group)->get_value(loc._idx);;
    }

    const V& at(const key_type& key) const
    {
        return const_cast<amt *>(this)->at(key);
    }

    bool contains(const key_type& key) const
    {
        return this->find(key) != this->end();
    }

    size_type count(const key_type& key) const
    {
        return this->contains(key) ? 1 : 0;
    }
    
    std::pair<iterator, iterator> equal_range(const key_type& key) 
    {
        auto it = find(key);
        if (it != end()) 
            return {it, std::next(it)};
        return {it, it};
    }

    std::pair<const_iterator, const_iterator> equal_range(const key_type& key) const
    {
        auto it = find(key);
        if (it != end()) 
            return {it, std::next(it)};
        return {it, it};
    }

    iterator find(const key_type& key)
    {
        auto loc = find_impl(key);
        return iterator(loc._group, loc._idx);
    }

    const_iterator find(const key_type& key) const
    {
        return const_cast<amt *>(this)->find(key);
    }

    mapped_type& operator[](const key_type& key)
    {
        auto loc = insert_impl(key);
        return reinterpret_cast<leaf_group_ptr>(loc._group)->get_value(loc._idx);
    }

    size_type bucket_count() const { return _size; } // not meaningful for amt
    float load_factor() const      { return 0.5f; }  // not meaningful for amt
    float max_load_factor() const  { return 1.0f; }  // not meaningful for amt
    void max_load_factor(float)    {}                // not meaningful for amt
    
    key_equal key_eq() const { return std::equal_to<K>(); }

private:
    void _init()
    {
        assert(_root == nullptr);
        _root = group::allocate_group(16, nullptr, 0, 0);
    }

    void _cleanup()
    {
        _root->deallocate_group();
        _root = nullptr;
        _last_insert = nullptr;
        _last_lookup = nullptr;
        _size = 0;
    }

    locator find_impl(const key_type& k)
    {
        bool last_matches = _last_lookup && _last_lookup->match(k);
        if (last_matches)
            return _last_lookup->find(k);
        auto loc = _root->find(k);
        _last_lookup = loc._group;
        return loc;
    }
   
    template <class Key, class Val>
    insert_locator insert_impl(Key&& k, Val&& v, bool force_insert = false) 
    {
        bool last_matches = _last_insert && _last_insert->match(k);

        auto loc = last_matches ? _last_insert->find_or_prepare_insert(k) : _root->find_or_prepare_insert(k);
        if (loc._created || force_insert)
        {
            reinterpret_cast<leaf_group_ptr>(loc._group)->set_value(loc._idx,  std::forward<Val>(v));
            ++_size;
        }
        if (!last_matches)
        {
            if (_last_insert)
                _last_insert->shrink();
            _last_insert = loc._group;
        }
        return loc;
    }

    void _erase(iterator it) 
    {
        assert(it != end());
        
        assert(0); //todo
    }
    
    void _erase(const_iterator cit) 
    {
        _erase(cit._it);
    }
    
    friend class iterator;

    group_ptr      _last_insert;  // always a leaf if set
    group_ptr      _last_lookup;  // always a leaf if set
    group_ptr      _root;  // never a leaf
    size_type      _size;
};



}  // namespace amt

#ifdef _MSC_VER
     #pragma warning(pop)  
#endif


#endif // amt_h_guard_
