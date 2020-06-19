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
#include <utility>
#include <type_traits>
#include <cassert>

namespace amt {


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
public:

    static_assert(std::is_integral<K>::value && std::is_unsigned<K>::value &&
                  (sizeof(K) == 2 || sizeof(K) == 4 || sizeof(K) == 8), "supporting 16, 32 or 64 bit unsigned integer keys.");

    static bool _equal(const K &a, const K &b) { return a == b; }
    
    using key_type    = K;
    using value_type  = std::pair<const K, const V>;
    using mapped_type = V;
    using size_type   = std::size_t;
    using difference_type = std::ptrdiff_t;
    using key_equal   = std::equal_to<K>;

    // --------------------- iterator ---------------------------------------
    class iterator 
    {
        friend class amt;

    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type        = typename amt::value_type;
        using reference         = value_type&;
        using const_reference   = const value_type&;
        using pointer           = value_type *;
        using const_pointer     = const value_type *;
        using difference_type   = typename amt::difference_type;

        iterator() {}

        reference operator*() const { assert(0); }

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
            assert(0);
            return false;
        }

        friend bool operator!=(const iterator& a, const iterator& b) 
        {
            return !(a == b);
        }

    private:
        // anonymous union to avoid uninitialized member warnings
        union {
            value_type v;
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
    amt() noexcept {}
    
    template <class InputIter>
    amt(InputIter first, InputIter last) : amt() 
    {
        insert(first, last);
    }

    amt(std::initializer_list<value_type> init) : amt(init.begin(), init.end())
    {}

    amt(const amt &o)
    {
        assert(0); //todo
    }

    amt(amt &&o) noexcept
    {
        assert(0); //todo
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
        assert(0); //todo
    }

    // --------------------- apis ------------------------------------------
    iterator begin() 
    {
        iterator it;
        assert(0); //todo
        return it;
    }
    
    iterator end() 
    {
        iterator it;
        assert(0); //todo
        return it; 
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
        assert(0); //todo
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
        return insert_impl(value.first, std::forward<VT>(value));
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

    template<class... Args>
    std::pair<iterator, bool> emplace(Args&&... args) 
    {
        return insert(value_type(std::forward<Args>(args)...));
    }
    
    template<class... Args>
    iterator emplace_hint(const_iterator hint, Args&&... args) 
    {
        return insert(hint, value_type(std::forward<Args>(args)...));        
    }
    
    template<class Key, class... Args>
    std::pair<iterator, bool> try_emplace(Key&& key, Args&&... args) 
    {
        return insert_impl(key, std::piecewise_construct, 
                                std::forward_as_tuple(std::forward<Key>(key)), 
                                std::forward_as_tuple(std::forward<Args>(args)...));
    }
    
    template<class Key, class... Args>
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
    emplace;
    emplace_hint;
    try_emplace;
    extract;
    merge;
#endif

    void reserve(size_type count) {}
    void rehash(size_type count) {}

    void swap(amt& other) noexcept
    {
        assert(0); //todo
    }

    V& at(const key_type& key)
    {
        auto it = this->find(key);
        if (it == this->end()) 
            throw std::out_of_range("amt at(): lookup non-existent key");
        return *it;
    }

    const V& at(const key_type& key) const
    {
        auto it = this->find(key);
        if (it == this->end()) 
            throw std::out_of_range("amt at(): lookup non-existent key");
        return *it;
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
        iterator it;
        assert(0); //todo
        return it;
    }

    const_iterator find(const key_type& key) const
    {
        return const_cast<amt *>(this)->find(key);
    }

    mapped_type& operator[](const key_type& key)
    {
        iterator it;
        assert(0); //todo
        return *it;
    }

    size_type bucket_count() const { return _size; } // not meaningful for amt
    float load_factor() const      { return 0.5f; }  // not meaningful for amt
    float max_load_factor() const  { return 1.0f; }  // not meaningful for amt
    void max_load_factor(float)    {}                // not meaningful for amt
    
    key_equal key_eq() const { return std::equal_to<K>(); }

private:
    template <class Key, class Val>
    std::pair<iterator, bool> insert_impl(Key&& k, Val&& v) 
    {
    }

    template <class Key = key_type, class... Args>
    std::pair<iterator, bool> try_emplace_impl(Key&& k, Args&&... args) 
    {
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

    size_type _size;
};



}  // namespace amt

#ifdef _MSC_VER
     #pragma warning(pop)  
#endif


#endif // amt_h_guard_
