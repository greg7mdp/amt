# AMT - Fast integer trie map

This is based on  Phil Bagwell's `Array Mapped Tree` from ["Fast And Space Efficient Trie Searches"](https://idea.popcount.org/2012-07-25-introduction-to-hamt/triesearches.pdf)

It supports unsigned integer keys (32, 64 or 128 bits) and has been specially optimised for the case where keys are somewhat monotonic.

## Overview

- **Header only**: the file `amt.h` is all you need

- **ordered container**

- **very fast** -  performance similar to hash table

- **mostly same interface as `std::map`**. However the value_type has both the `key` and the `value` const.

- **low memory usage**


