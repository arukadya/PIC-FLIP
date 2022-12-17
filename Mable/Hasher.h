//
//  Hasher.h
//  HashMapClasstaring
//
//  Created by 須之内俊樹 on 2022/10/29.
//

#ifndef Hasher_h
#define Hasher_h
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <unordered_map>
#include <type_traits>//std::remove_cvref_t(C++20)
#include <functional>//std::hash

template <int N>
struct ArrayHasher{
    template <typename T>
    std::size_t operator()(T const& key) const{
        //std::unorderd_map<key,value>では、キーにはハッシュ値をデフォで計算できる型しか使えないので、自力で作った型やライブラリから引っ張ってきた型はハッシュ値を計算する関数を自作して３個目の引数にわたす。このとき関数テンプレートで渡す。operator()はこのハッシュ値を計算する関数テンプレートをオーバーロードしてると思われる。std::size_t型は、32bit環境なら32bit,64bit環境なら64bitの符号なし整数型で、配列やコンテナのサイズ、シード値に使われる。
        std::hash<std::remove_cvref_t<decltype(key[0])>> hasher;
        //テンプレートの引数にconstを入れると、コンパイルの最適化の時におかしなことが起きるらしい。remove_cvref_tはconstとかを剥ぎ取れる(swiftでやったアンラップ)。また、型がはっきりしてないと困るので、テンプレートの入れ子構造みたいなことはできないはずだが、decltype()は()の中の型を動的に取得でき、正常に作動する。
        std::size_t seed = 0;
        for(int i=0;i<N;++i){
            seed ^= hasher(key[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

struct StructHasher2d{
    template <typename T>
    std::size_t operator()(T const& key)const{
        auto const& [elem0,elem1] = key;
        std::hash<std::remove_cvref_t<decltype(elem0)>> hasher0;
        std::hash<std::remove_cvref_t<decltype(elem1)>> hasher1;
        std::size_t seed = 0;
        seed ^= hasher0(elem0) + 0x9e3779b9 +(seed << 6) + (seed >> 2);
        seed ^= hasher1(elem1) + 0x9e3779b9 +(seed << 6) + (seed >> 2);
        return seed;
    }
};
struct StructHasher3d{
    template <typename T>
    std::size_t operator()(T const& key)const{
        auto const& [elem0,elem1,elem2] = key;
        std::hash<std::remove_cvref_t<decltype(elem0)>> hasher0;
        std::hash<std::remove_cvref_t<decltype(elem1)>> hasher1;
        std::hash<std::remove_cvref_t<decltype(elem2)>> hasher2;
        std::size_t seed = 0;
        seed ^= hasher0(elem0) + 0x9e3779b9 +(seed << 6) + (seed >> 2);
        seed ^= hasher1(elem1) + 0x9e3779b9 +(seed << 6) + (seed >> 2);
        seed ^= hasher1(elem2) + 0x9e3779b9 +(seed << 6) + (seed >> 2);
        return seed;
    }
};
#endif /* Hasher_h */
