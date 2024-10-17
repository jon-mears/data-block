#ifndef DATA_BLOCK_HPP
#define DATA_BLOCK_HPP

#include <concepts>
#include <initializer_list>
#include <iostream>
#include <type_traits>

#include "index-range.hpp"
#include "operation-constraints.hpp"

using namespace OpTraits;

template <typename T, size_t N>
struct RawVector;

template <typename T, size_t R, size_t C>
struct RawMatrix;

template <typename T, template <typename, size_t...> typename TEMPLATE>
inline constexpr bool is_instance_of_v = false;

template <template <typename, size_t...> typename TEMPLATE, typename T, size_t... ARGS>
inline constexpr bool is_instance_of_v<TEMPLATE<T, ARGS...>, TEMPLATE> = true;

template <typename T>
struct FillInfo {
    T const& fill;
};

template <typename T1, typename T2>
concept IsExplicitlyConvertible = requires (T1 a) {
    static_cast<T2>(a);
};

template <typename T1, typename T2>
concept IsImplicitlyConvertible = (std::is_convertible_v<T1, T2>);

// datablocks and datablockviews

// base
template <size_t D1, size_t... Ds> requires (sizeof...(Ds) == 0)
size_t Block2FlatIndex(size_t index) {
    return index;
}

// recursive
template <size_t D1, size_t... Ds> requires (sizeof...(Ds) != 0)
size_t Block2FlatIndex(size_t first, std::convertible_to<size_t> auto const... rest) {
    return first * (Ds * ...) + Block2FlatIndex<Ds...>(rest...);
}

template <typename T, size_t... Ds>
struct DataBlock;

template <typename T, size_t... Ds>
struct DataBlockView;

template <typename T, size_t... Ds>
using ConstDataBlockView = DataBlockView<const T, Ds...>;

template <typename T, typename DB, size_t ARGUMENT_COUNT, size_t UNSPECIFIED_COUNT, typename IS = std::make_index_sequence<UNSPECIFIED_COUNT>>
struct data_block_view_type;

template <typename T, typename DB, size_t ARGUMENT_COUNT, size_t UNSPECIFIED_COUNT, size_t... Is>
struct data_block_view_type<T, DB, ARGUMENT_COUNT, UNSPECIFIED_COUNT, std::index_sequence<Is...>> {
    typedef DataBlockView<T, DB::shape[ARGUMENT_COUNT + Is]...> type;
};

template <size_t ARGUMENT_COUNT, size_t UNSPECIFIED_COUNT, typename T, size_t... Ds>
using data_block_view_t = typename data_block_view_type<T, DataBlock<T, Ds...>, ARGUMENT_COUNT, UNSPECIFIED_COUNT>::type;

// multiple potential DataBlockView types

template <typename T, size_t... Ds>
struct DataBlockView {
public:
    using type = DataBlockView<T, Ds...>;
    using value_type = T;

    using view_type = DataBlockView<T, Ds...>;
    using non_view_type = DataBlock<T, Ds...>;
private:

    T* mpData;
    int mStep;

public:

    static constexpr size_t NElems() {
        return DataBlock<T, Ds...>::nelems;
    }

    static constexpr size_t NDims() {
        return DataBlock<T, Ds...>::ndims;
    }

    static constexpr size_t Shape(size_t idx) {
        return DataBlock<T, Ds...>::shape[idx];
    }

    template <typename... Ts>
    requires (sizeof...(Ts) >= sizeof...(Ds))
    inline T& operator()(Ts... indices) {
        return mpData[Block2FlatIndex<Ds...>(indices...)];
    }

    template <typename... Ts>
    requires (sizeof...(Ts) >= sizeof...(Ds))
    inline T const& operator()(Ts... indices) const {
        return mpData[Block2FlatIndex<Ds...>(indices...)];
    }

    template <typename... Ts>
    requires (sizeof...(Ts) < sizeof...(Ds))
    inline auto View(Ts... indices) {
        return data_block_view_t<sizeof...(Ts), sizeof...(Ds) - sizeof...(Ts), T, Ds...>(*this, indices...);
    }

    template <typename... Ts>
    requires (sizeof...(Ts) < sizeof...(Ds))
    inline const auto View(Ts... indices) const {
        return data_block_view_t<sizeof...(Ts), sizeof...(Ds) - sizeof...(Ts), T, Ds...>(*this, indices...);
    }

    inline T& operator[](size_t i) {
        return mpData[i * mStep];
    }

    inline T const& operator[](size_t i) const {
        return mpData[i * mStep];
    }

    template <size_t... D2s>
    DataBlockView(DataBlock<T, D2s...>& rDB, SkipIndexRange const& crSIR) :
        mpData{ &rDB[crSIR.start] }, mStep(crSIR.skip) { }

    template <size_t... D2s>
    DataBlockView(DataBlock<T, D2s...>& rDB, IndexRange const& crIR) :
        mpData{ &rDB[crIR.start] }, mStep{ 1 } { }

    template <size_t... D2s>
    DataBlockView(DataBlock<std::remove_const_t<T>, D2s...> const& crDB, SkipIndexRange const& crSIR) :
        mpData{ &crDB[crSIR.start] }, mStep{ crSIR.skip } { }

    template <size_t... D2s>
    DataBlockView(DataBlock<std::remove_const_t<T>, D2s...> const& crDB, IndexRange const& crIR) :
        mpData{ &crDB[crIR.start] }, mStep{ 1 } { }

    template <size_t... D2s>
    DataBlockView(DataBlock<std::remove_const_t<T>, D2s...>&& rrDB, SkipIndexRange const& crSIR) = delete;

    template <size_t... D2s>
    DataBlockView(DataBlock<std::remove_const_t<T>, D2s...>&& rrDB, IndexRange const& crIR) = delete;

    DataBlockView(T* pData, int step = 1) : mpData{ pData }, mStep{ step } { }

    operator DataBlockView<std::add_const_t<T>, Ds...>() {
        return DataBlockView<std::add_const_t<T>, Ds...>(mpData, mStep);
    }

    // assignment operators

    template <size_t N> requires (N == NElems())
    DataBlockView& operator=(T const (&crArr)[N]) {
        for (int i = 0; i < NElems(); ++i) {
            operator[](i) = crArr[i];
        }

        return *this;
    }

    DataBlockView& operator=(DataBlock<T, Ds...> const& crDB) {
        for (int i = 0; i < NElems(); ++i) {
            operator[](i) = crDB[i];
        }

        return *this;
    }

    DataBlockView& operator=(DataBlockView<T, Ds...> const& crDBV) {
        for (int i = 0; i < NElems(); ++i) {
            operator[](i) = crDBV[i];
        }

        return *this;
    }
};

template <typename T>
concept NotADataBlock = (!is_instance_of_v<T, DataBlock> && !is_instance_of_v<T, DataBlockView>);

template <typename T>
concept ADataBlock = (is_instance_of_v<T, DataBlock> || is_instance_of_v<T, DataBlockView>);

template <typename T>
constexpr size_t NElems();

template <typename T, size_t... Ds>
struct DataBlock {
public:
    typedef T value_type;
protected:

    T mData[(Ds * ...)];

    // introduce _IPairwiseExpressions (w/ a function) and _IBroadcastableExpressions (w/ a function)

    class _IMultiExpression {
    private:
        DataBlock<T, Ds...>& mrDB;
    public:
        _IMultiExpression(DataBlock<T, Ds...>& rDB) : mrDB{ rDB } { }
    };

    friend class _IMultiExpression;

public:

    static constexpr size_t NElems() {
        return nelems;
    }

    static constexpr size_t NDims() {
        return ndims;
    }

    static constexpr size_t Shape(size_t idx) {
        return shape[idx];
    }

    DataBlock(std::convertible_to<T> auto... data)
        requires (sizeof...(data) == DataBlock<T, Ds...>::nelems) : mData(data...) { }

    DataBlock(FillInfo<T> const& crFillInfo) {
        for (int i = 0; i < nelems; ++i) {
            mData[i] = crFillInfo.fill;
        }
    }

    DataBlock() { }

    DataBlock(DataBlock<T, Ds...> const& crDB) {
        memcpy(mData, &crDB[0], sizeof(mData));
    }

    template <std::convertible_to<T> U> requires (!std::is_same_v<T, U>)
    DataBlock(DataBlock<U, Ds...> const& crDB) {
        for (int i = 0; i < nelems; ++i) {
            mData[i] = crDB[i];
        }
    }

    DataBlock(DataBlockView<T, Ds...> const& crDBV) {
        memcpy(mData, &crDBV[0], sizeof(T) * nelems);
    }

    template <std::convertible_to<T> U> requires (!std::is_same_v<T, U>)
    DataBlock(DataBlockView<U, Ds...> const& crDBV) {
        for (int i = 0; i < nelems; ++i) {
            mData[i] = crDBV[i];
        }
    }

    // constructor that takes data, at least one of which is another DataBlock that is expanded
    template <typename... Ts>
    DataBlock(Ts const&... crDBs) requires ((::NElems<Ts>() + ...) <= NElems() && 
        ((std::convertible_to<Ts, T> || ADataBlock<Ts>) && ...) && (ADataBlock<Ts> || ...)) : mData{} {

        auto Helper = []<typename U>(DataBlock<T, Ds...>& rDest, U const& crSource, int dest_idx) {

            if constexpr (ADataBlock<U>) {
                for (int source_idx = 0; source_idx < crSource.NElems(); ++source_idx) {
                    rDest[dest_idx++] = crSource[source_idx];
                }
            }

            else {
                rDest[dest_idx++] = crSource;
            }

            return dest_idx;
        };
        
        int dest_idx = 0;

        ((dest_idx = Helper(*this, crDBs, dest_idx)), ...);
    }

    DataBlock<T, Ds...>& operator=(DataBlock<T, Ds...> const& crDB) {
        memcpy(mData, crDB.mData, sizeof(crDB.mData));
        return *this;
    }

    template <std::convertible_to<T> U> requires (!std::is_same_v<T, U>)
    DataBlock<T, Ds...>& operator=(DataBlock<U, Ds...> const& crDB) {
        for (int i = 0; i < nelems; ++i) {
            mData[i] = crDB[i];
        }

        return *this;
    }

    DataBlock<T, Ds...>& operator=(DataBlockView<T, Ds...> const& crDBV) {
        memcpy(mData, &crDBV((Ds, 0)...), sizeof(T) * crDBV.nelems);

        return *this;
    }

    template <std::convertible_to<T> U> requires (!std::is_same_v<T, U>)
    DataBlock<T, Ds...>& operator=(DataBlockView<U, Ds...> const& crDBV) {
        for (int i = 0; i < nelems; ++i) {
            mData[i] = crDBV[i];
        }

        return *this;
    }

    _IMultiExpression IMultiExpression() { return _IMultiExpression{ *this }; }

    static constexpr size_t ndims{ sizeof...(Ds) };
    static constexpr size_t nelems{ (Ds * ...) };
    static constexpr size_t shape[]{ Ds... };

    template <typename... Ts>
        requires (sizeof...(Ts) >= sizeof...(Ds))
    inline T& operator()(Ts... indices) {
        return mData[Block2FlatIndex<Ds...>(indices...)];
    }

    template <typename... Ts>
        requires (sizeof...(Ts) >= sizeof...(Ds))
    inline T const& operator()(Ts... indices) const {
        return mData[Block2FlatIndex<Ds...>(indices...)];
    }

    template <typename... Ts>
        requires (sizeof...(Ts) < sizeof...(Ds))
    auto View(Ts... indices) {
        return data_block_view_t<sizeof...(Ts), sizeof...(Ds) - sizeof...(Ts), T, Ds...>(*this, indices...);
    }

    template <typename... Ts>
        requires (sizeof...(Ts) < sizeof...(Ds))
    auto View(Ts... indices) const {
        return data_block_view_t<sizeof...(Ts), sizeof...(Ds) - sizeof...(Ts), T, Ds...>(*this, indices...);
    }

    inline T& operator[](size_t idx) {
        return mData[idx];
    }

    inline T const& operator[](size_t idx) const {
        return mData[idx];
    }

    inline DataBlockView<T, Ds...> operator[](SkipIndexRange const& crSIR) {
        return DataBlockView<T, Ds...>(*this, crSIR);
    }

    inline DataBlockView<T, Ds...> operator[](IndexRange const& crIR) {
        return DataBlockView<T, Ds...>(*this, crIR);
    }

    template <int START, int STOP, int STEP> requires ValidIndexRangeV<START, STOP, STEP>
    DataBlockView<T, (STOP - START) / STEP + ((STOP - START) % STEP ? 1 : 0)> operator[]
    (TemplatedIndexRange<START, STOP, STEP> tip) {
        return DataBlockView<T, (STOP - START) / STEP + ((STOP - START) % STEP ? 1 : 0)>(
            &mData[START], STEP
        );
    }

    template <int START, int STOP, int STEP> requires ValidIndexRangeV<START, STOP, STEP>
    ConstDataBlockView<T, (STOP - START) / STEP + ((STOP - START) % STEP ? 1 : 0)> operator[]
    (TemplatedIndexRange<START, STOP, STEP> tip) const {
        return ConstDataBlockView<T, (STOP - START) / STEP + ((STOP - START) % STEP ? 1 : 0)>(
            &mData[START], STEP
        );
    }

    explicit operator bool() requires (IsExplicitlyConvertible<T, bool>) {
        for (int i = 0; i < nelems; ++i) {
            if (static_cast<bool>(mData[i])) {}
            else { return false; }
        }

        return true;
    }

    template <typename U>
    operator DataBlock<U, Ds...>() requires (IsImplicitlyConvertible<T, U>) {
        DataBlock<U, Ds...> result;

        for (int i = 0; i < nelems; ++i) {
            result[i] = mData[i];
        }

        return result;
    }

    template <typename U>
    explicit operator DataBlock<U, Ds...>() requires (!IsImplicitlyConvertible<T, U>&& IsExplicitlyConvertible<T, U>) {
        DataBlock<U, Ds...> result;

        for (int i = 0; i < nelems; ++i) {
            result[i] = static_cast<U>(mData[i]);
        }

        return result;
    }

    // Matrix

    auto Row(size_t i) requires (sizeof...(Ds) == 2) {
        return DataBlockView<T, shape[1]>(*this, { .start = i * shape[1], .end = (i + 1) * shape[1],
            .skip = 1 });
    }

    auto Column(size_t i) requires (sizeof...(Ds) == 2) {
        return DataBlockView<T, shape[0]>(*this, { .start = i,
            .end = (shape[0] - 1) * shape[1] + i + 1, .skip = shape[1] });
    }

    // Square Matrix

    auto Diagonal() requires (sizeof...(Ds) == 2 && (Ds == ...)) {
        return DataBlockView<T, shape[0]>(*this, 
            { .start = 0, .end = nelems, .skip = shape[1] + 1 });
    }

    static DataBlock<T, Ds...> Identity() requires (sizeof...(Ds) == 2 && (Ds == ...)) {
        DataBlock<T, Ds...> result({ .fill = 0 });

        for (int i = 0; i < nelems; i += shape[0] + 1) {
            result[i] = 1;
        }

        return result;
    }
};

template <typename T>
constexpr size_t NElems() {
    if constexpr (ADataBlock<T>) {
        return T::NElems();
    }

    else {
        return 1;
    }
}

template <size_t F = 0, size_t... R>
struct size_t_pack;

template <typename T1, typename T2>
struct ConsistentDims;

template <size_t CPTY, size_t SIZE, typename SZP, typename IS = std::make_index_sequence<CPTY - SIZE>>
struct Size_T_Arr;

template <size_t CPTY, size_t SIZE, size_t... DATA, size_t... Is>
struct Size_T_Arr<CPTY, SIZE, size_t_pack<DATA...>, std::index_sequence<Is...>> {
    typedef size_t_pack<DATA..., (Is, 0)...> arr;
};

template <>
struct ConsistentDims<size_t_pack<0>, size_t_pack<0>> {
    static constexpr bool value = true;
};

template <size_t F1, size_t F2, size_t... R1, size_t... R2>
struct ConsistentDims<size_t_pack<F1, R1...>, size_t_pack<F2, R2...>> {
    static constexpr bool value = ConsistentDims<size_t_pack<R1...>, size_t_pack<R2...>>::value && (F1 == F2 || (F1 == 1 || F2 == 1) || (F1 == 0 || F2 == 0));
};

template <typename T1, typename T2>
struct wider {
};

template <typename T>
struct wider<T, T> {
    typedef T type;
};

template <>
struct wider<short, int> {
    typedef int type;
};

template <typename T1, typename T2>
using wider_t = wider<T1, T2>::type;

template <typename T1, typename T2, typename SZP1, typename SZP2>
struct broadcast_type_info;

template <typename T1, typename T2, size_t... D1, size_t... D2>
struct broadcast_type_info<T1, T2, size_t_pack<D1...>, size_t_pack<D2...>> {
    typedef DataBlock<wider_t<T1, T2>, (D1 > D2 ? D1 : D2)...> type;
};

template <typename T1, typename T2>
struct broadcast_info;

template <typename T1, typename T2, size_t... R1, size_t... R2>
struct broadcast_info<DataBlock<T1, R1...>, DataBlock<T2, R2...>> {
    typedef typename broadcast_type_info<T1, T2, typename Size_T_Arr<(sizeof...(R1) > sizeof...(R2) ? sizeof...(R1) : sizeof...(R2)), sizeof...(R1),
        size_t_pack<R1...>>::arr, typename Size_T_Arr<(sizeof...(R1) > sizeof...(R2) ? sizeof...(R1) : sizeof...(R2)), sizeof...(R2),
        size_t_pack<R2...>>::arr>::type type;
};

template <typename T1, typename T2>
using broadcast_t = broadcast_info<T1, T2>::type;

struct AnyType { };

template <typename T, typename VT = AnyType, size_t... Ds>
constexpr bool IsADataBlockV = false;

template <typename VT, size_t... Ds>
constexpr bool IsADataBlockV<DataBlock<VT, Ds...>, VT, Ds...> = true;

template <typename VT, size_t... Ds>
constexpr bool IsADataBlockV<DataBlock<VT, Ds...>, AnyType, Ds...> = true;

template <typename VT, size_t... Ds>
constexpr bool IsADataBlockV<DataBlock<VT, Ds...>, AnyType> = true;

template <typename VT, size_t... Ds>
constexpr bool IsADataBlockV<DataBlockView<VT, Ds...>, VT, Ds...> = true;

template <typename VT, size_t... Ds>
constexpr bool IsADataBlockV<DataBlockView<VT, Ds...>, AnyType, Ds...> = true;

template <typename VT, size_t... Ds>
constexpr bool IsADataBlockV<DataBlockView<VT, Ds...>, AnyType> = true;

template <typename T, typename VT = AnyType, size_t... Ds>
concept IsADataBlockC = IsADataBlockV<T, VT, Ds...>;

template <IsADataBlockC T1, IsADataBlockC T2>
constexpr bool HasSameDimensionsV = false;

template <template <typename, size_t...> typename DB1, template <typename, size_t...> typename DB2, typename VT1, typename VT2, size_t... Ds>
constexpr bool HasSameDimensionsV<DB1<VT1, Ds...>, DB2<VT2, Ds...>> = true;

constexpr bool HasSameDimensions(IsADataBlockC auto const& crLHS, IsADataBlockC auto const& crRHS) {
    using T1 = std::remove_cvref_t<decltype(crLHS)>;
    using T2 = std::remove_cvref_t<decltype(crRHS)>;

    return HasSameDimensionsV<T1, T2>;
}

template <IsADataBlockC T1, IsADataBlockC T2>
constexpr bool HasSameNumberOfDimensionsV = false;

template <template <typename, size_t...> typename DB1, template <typename, size_t...> typename DB2, typename VT1, typename VT2, size_t... D1s, size_t... D2s>
constexpr bool HasSameNumberOfDimensionsV<DB1<VT1, D1s...>, DB2<VT2, D2s...>> = sizeof...(D1s) == sizeof...(D2s);

constexpr bool HasSameNumberOfDimensions(IsADataBlockC auto const& crLHS, IsADataBlockC auto const& crRHS) {
    using T1 = std::remove_cvref_t<decltype(crLHS)>;
    using T2 = std::remove_cvref_t<decltype(crRHS)>;

    return HasSameNumberOfDimensionsV<T1, T2>;
}

template <typename T>
inline constexpr bool IsATypeV = false;

template <template <typename> typename TMP, typename T>
inline constexpr bool IsATypeV<TMP<T>> = std::is_base_of_v<Type<T>, TMP<T>> || std::is_same_v<Type<T>, TMP<T>>;

template <typename T>
concept IsATypeC = IsATypeV<T>;

namespace Detail {
    template <typename T, typename U>
    struct DataBlockWithHelper;

    template <IsATypeC T, size_t... Ds>
    struct DataBlockWithHelper<T, Dimensions<Ds...>> {
        using type = DataBlock<typename T::type, Ds...>;
    };
}

template <auto X, auto Y>
using DataBlockWith = Detail::DataBlockWithHelper<std::remove_const_t<decltype(X)>, std::remove_const_t<decltype(Y)>>::type;

template <typename T>
concept HasValueTypeC = requires {
    typename T::value_type;
};

// unary plus
auto operator+(IsADataBlockC auto const& crDB) requires (CanUnaryPlus(ValueType(crDB))) {
    DataBlockWith<UnaryPlusType(ValueType(crDB)), Dimensions(crDB)> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = +crDB[i]; 
    }

    return result;
}

// unary minus
auto operator-(IsADataBlockC auto const& crDB) requires (CanUnaryMinus(ValueType(crDB))) {
    DataBlockWith<ValueType(crDB), Dimensions(crDB)> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = -crDB[i];
    }

    return result;
}

// dereference
template <typename VT, template <typename, size_t...> typename DB, size_t... Ds>
auto operator*(DB<VT, Ds...> const& crDB)
    requires (ADataBlock<DB<VT, Ds...>>) {
    static_assert(CanDereference<VT>, "Cannot dereference a DataBlock whose elements cannot be dereferenced.");

    DataBlock<dereference_t<VT>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = *crDB[i];
    }

    return result;
}

// default address-of; see IMultiExpression for constructing a new DataBlock containing the addresses of all elements of the
// original DataBlock

template <CanBitwiseNot VT, template <typename, size_t...> typename DB, size_t... Ds>
DataBlock<bitwise_not_t<VT>, Ds...> operator~(DB<VT, Ds...> const& crDB)
    requires (ADataBlock<DB<VT, Ds...>>) {
    //static_assert(CanBitwiseNot<VT>, "Cannot apply bitwise-not to a DataBlock whose elements cannot have bitwise-not applied to them.");

    DataBlock<bitwise_not_t<VT>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = ~crDB[i];
    }

    return result;
}

template <typename VT, template <typename, size_t...> typename DB, size_t... Ds>
bool operator!(DB<VT, Ds...> const& crDB)
    requires (ADataBlock<DB<VT, Ds...>>) {
    static_assert(CanLogicalNot<VT>, "Cannot apply logical-not to a DataBlock whose elements cannot have logical-not applied to them.");
    static_assert(IsExplicitlyConvertible<logical_not_t<VT>, bool>, "Cannot apply logical-not to a DataBlock whose elements do not yield a"
        " type explicilty convertible to bool when logical-not is applied to them.");

    for (int i = 0; i < crDB.nelems; ++i) {
        if (!crDB[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT, size_t... Ds>
auto& operator++(DataBlock<VT, Ds...>& rDB) {
    static_assert(CanPrefixIncrement<VT>, "Cannot prefix increment a DataBlock whose elements cannot be prefix incremented.");

    for (int i = 0; i < rDB.nelems; ++i) {
        ++rDB[i];
    }

    return rDB;
}

template <typename VT, size_t... Ds>
auto& operator++(DataBlockView<VT, Ds...>&& rDBV) {
    static_assert(CanPrefixIncrement<VT>, "Cannot prefix increment a DataBlock whose elements cannot be prefix incremented.");

    for (int i = 0; i < rDBV.nelems; ++i) {
        ++rDBV[i];
    }

    return rDBV;
}

template <typename VT, size_t... Ds>
auto operator++(DataBlock<VT, Ds...>& rDB, int) {
    static_assert(CanPostfixIncrement<VT>, "Cannot postfix increment a DataBlock whose elements cannot be postfix incremented.");

    DataBlock<VT, Ds...> db = rDB;

    for (int i = 0; i < rDB.nelems; ++i) {
        rDB[i]++;
    }

    return db;
}

template <typename VT, size_t... Ds>
auto operator++(DataBlockView<VT, Ds...>&& rDBV, int) {
    static_assert(CanPostfixIncrement<VT>, "Cannot postfix increment a DataBlock whose elements cannot be postfix incremented.");

    DataBlock<VT, Ds...> db = rDBV;

    for (int i = 0; i < rDBV.nelems; ++i) {
        rDBV[i]++;
    }

    return db;
}

template <typename VT, size_t... Ds>
auto& operator--(DataBlock<VT, Ds...>& rDB) {
    static_assert(CanPrefixDecrement<VT>, "Cannot prefix decrement a DataBlock whose elements cannot be prefix decremented.");

    for (int i = 0; i < rDB.nelems; ++i) {
        --rDB[i];
    }

    return rDB;
}

template <typename VT, size_t... Ds>
auto& operator--(DataBlockView<VT, Ds...>&& rDBV) {
    static_assert(CanPrefixDecrement<VT>, "Cannot prefix decrement a DataBlock whose elements cannot be prefix decremented.");

    for (int i = 0; i < rDBV.nelems; ++i) {
        --rDBV[i];
    }

    return rDBV;
}

template <typename VT, size_t... Ds>
auto operator--(DataBlock<VT, Ds...>& rDB, int) {
    static_assert(CanPostfixDecrement<VT>, "Cannot postfix decrement a DataBlock whose elements cannot be postfix decremented.");

    DataBlock<VT, Ds...> db = rDB;

    for (int i = 0; i < rDB.nelems; ++i) {
        rDB[i]--;
    }

    return db;
}

template <typename VT, size_t... Ds>
auto operator--(DataBlockView<VT, Ds...>&& rDBV, int) {
    static_assert(CanPostfixDecrement<VT>, "Cannot postfix decrement a DataBlock whose elements cannot be postfix decremented.");

    DataBlock<VT, Ds...> db = rDBV;

    for (int i = 0; i < rDBV.nelems; ++i) {
        rDBV[i]--;
    }

    return db;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
DataBlock<sum_t<VT1, VT2>, Ds...> operator+(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanSum<VT1, VT2>, "Cannot sum DataBlocks whose element types cannot be summed.");

    DataBlock<sum_t<VT1, VT2>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] + crRHS[i];
    }

    return result;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
DB<sum_t<VT, ST>, Ds...> operator+(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanSum<VT, ST>, "Cannot sum a DataBlock and some other type when the elements of"
        " the other type and the elements of the DataBlock cannot be summed.");

    DB<sum_t<VT, ST>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] + crRHS;
    }

    return result;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
auto operator+(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanSum<ST, VT>, "Cannot sum a DataBlock and some other type when the elements of"
        " the other type and the elements of the DataBlock cannot be summed.");

    DataBlock<sum_t<ST, VT>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS + crRHS[i];
    }

    return result;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator+=(DB1<VT1, Ds...>&& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<sum_t<VT1, VT2>, VT1>, "Cannot apply compound-sum-assignment to two DataBlock when "
        "the sum type of the element types of the DataBlocks are not implicitly convertable to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] += crRHS[i];
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator+=(DB1<VT1, Ds...>& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<sum_t<VT1, VT2>, VT1>, "Cannot apply compound-sum-assignment to two DataBlock when "
        "the sum type of the element types of the DataBlocks are not implicitly convertable to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] += crRHS[i];
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator+=(DB<VT, Ds...>&& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<sum_t<VT, ST>, VT>, "Cannot apply compound-sum-assignment to a DataBlock and another type when "
        "the sum type of the right-hand type and the element type of the DataBlock is not implicitly convertable to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] += crRHS;
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator+=(DB<VT, Ds...>& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<sum_t<VT, ST>, VT>, "Cannot apply compound-sum-assignment to a DataBlock and another type when "
        "the sum type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] += crRHS;
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto operator-(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanSubtract<VT1, VT2>, "Cannot subtract DataBlocks whose element types cannot be subtracted.");

    DataBlock<difference_t<VT1, VT2>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] - crRHS[i];
    }

    return result;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
auto operator-(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanSubtract<VT, ST>, "Cannot subtract a DataBlock and some other type when"
        " the other type and the elements of the DataBlock cannot be subtracted.");

    DataBlock<difference_t<VT, ST>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] - crRHS;
    }

    return result;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
auto operator-(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanSum<ST, VT>, "Cannot subtract a DataBlock and some other type when"
        " the other type and the elements of the DataBlock cannot be subtracted.");

    DataBlock<difference_t<ST, VT>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS - crRHS[i];
    }

    return result;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator-=(DB1<VT1, Ds...>&& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<difference_t<VT1, VT2>, VT1>, "Cannot apply compound-subtraction-assignment to two DataBlock when "
        "the difference type of the element types of the DataBlocks are not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] -= crRHS[i];
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator-=(DB1<VT1, Ds...>& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<difference_t<VT1, VT2>, VT1>, "Cannot apply compound-subtraction-assignment to two DataBlocks when "
        "the difference type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] -= crRHS[i];
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator-=(DB<VT, Ds...>&& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<difference_t<VT, ST>, VT>, "Cannot apply compound-subtraction-assignment to a DataBlock and another type when "
        "the difference type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] -= crRHS;
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator-=(DB<VT, Ds...>& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<difference_t<VT, ST>, VT>, "Cannot apply compound-subtraction-assignment to a DataBlock and another type when "
        "the difference type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] -= crRHS;
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
DataBlock<product_t<VT1, VT2>, Ds...> operator*(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanMultiply<VT1, VT2>, "Cannot multiply DataBlocks whose element types cannot be multiplied.");

    DataBlock<product_t<VT1, VT2>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] * crRHS[i];
    }

    return result;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
auto operator*(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanMultiply<VT, ST>, "Cannot multiply a DataBlock and some other type when"
        " the other type and the elements of the DataBlock cannot be multiplied.");

    DataBlock<product_t<VT, ST>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] * crRHS;
    }

    return result;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
auto operator*(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanMultiply<ST, VT>, "Cannot multiply a DataBlock and some other type when"
        " the other type and the elements of the DataBlock cannot be multiplied.");

    DataBlock<product_t<ST, VT>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS * crRHS[i];
    }

    return result;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator*=(DB1<VT1, Ds...>&& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<product_t<VT1, VT2>, VT1>, "Cannot apply compound-multiplication-assignment to two DataBlocks when "
        "the product type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] *= crRHS[i];
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator*=(DB1<VT1, Ds...>& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>> && ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<product_t<VT1, VT2>, VT1>, "Cannot apply compound-multiplication-assignment to two DataBlocks when "
        "the product type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] *= crRHS[i];
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator*=(DB<VT, Ds...>&& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<product_t<VT, ST>, VT>, "Cannot apply compound-multiplication-assignment to a DataBlock and another type when "
        "the product type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] *= crRHS;
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator*=(DB<VT, Ds...>& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<product_t<VT, ST>, VT>, "Cannot apply compound-multiplication-assignment to a DataBlock and another type when "
        "the product type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] *= crRHS;
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto operator/(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanDivide<VT1, VT2>, "Cannot divide DataBlocks whose element types cannot be divided.");

    DataBlock<quotient_t<VT1, VT2>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] / crRHS[i];
    }

    return result;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
auto operator/(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanDivide<VT, ST>, "Cannot divide a DataBlock and some other type when"
        " the other type and the elements of the DataBlock cannot be divided.");

    DataBlock<quotient_t<VT, ST>, Ds...> result;
    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] / crRHS;
    }

    return result;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
auto operator/(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanDivide<ST, VT>, "Cannot divide a DataBlock and some other type when"
        " the other type and the elements of the DataBlock cannot be divided.");

    DataBlock<quotient_t<ST, VT>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS / crRHS[i];
    }

    return result;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator/=(DB1<VT1, Ds...>&& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<quotient_t<VT1, VT2>, VT1>, "Cannot apply compound-division-assignment to two DataBlocks when "
        "the quotient type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] /= crRHS[i];
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator/=(DB1<VT1, Ds...>& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<quotient_t<VT1, VT2>, VT1>, "Cannot apply compound-division-assignment to two DataBlocks when "
        "the quotient type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] /= crRHS[i];
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator/=(DB<VT, Ds...>&& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<quotient_t<VT, ST>, VT>, "Cannot apply compound-division-assignment to a DataBlock and another type when "
        "the quotient type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] /= crRHS;
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator/=(DB<VT, Ds...>& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<quotient_t<VT, ST>, VT>, "Cannot apply compound-division-assignment to a DataBlock and another type when "
        "the quotient type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] /= crRHS;
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto operator%(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanModulus<VT1, VT2>, "Cannot take the modulus of DataBlocks for which the modulus of their element types cannot be taken.");

    DataBlock<modulus_t<VT1, VT2>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] % crRHS[i];
    }

    return result;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
auto operator%(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanModulus<VT, ST>, "Cannot take the modulus of a DataBlock and some other type when"
        " the modulus of the other type and the elements of the DataBlock cannot be taken.");

    DataBlock<modulus_t<VT, ST>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] % crRHS;
    }

    return result;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
auto operator%(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanModulus<ST, VT>, "Cannot take the modulus of a DataBlock and some other type when"
        " the modulus of the other type and the elements of the DataBlock cannot be taken.");

    DataBlock<modulus_t<ST, VT>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS % crRHS[i];
    }

    return result;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator%=(DB1<VT1, Ds...>&& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<modulus_t<VT1, VT2>, VT1>, "Cannot apply compound-modulus-assignment to two DataBlocks when "
        "the modulus type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] %= crRHS[i];
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator%=(DB1<VT1, Ds...>& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<modulus_t<VT1, VT2>, VT1>, "Cannot apply compound-modulus-assignment to two DataBlocks when "
        "the modulus type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] %= crRHS[i];
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator%=(DB<VT, Ds...>&& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<modulus_t<VT, ST>, VT>, "Cannot apply compound-modulus-assignment to a DataBlock and another type when "
        "the modulus type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] %= crRHS;
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator%=(DB<VT, Ds...>& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<modulus_t<VT, ST>, VT>, "Cannot apply compound-modulus-assignment to a DataBlock and another type when "
        "the modulus type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] %= crRHS;
    }

    return rLHS;
}


template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto operator^(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanBitwiseXor<VT1, VT2>, "Cannot apply bitwise xor to DataBlocks for which bitwise xor cannot be applied to their element types.");

    DataBlock<bitwise_xor_t<VT1, VT2>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] ^ crRHS[i];
    }

    return result;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
auto operator^(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanBitwiseXor<VT, ST>, "Cannot apply bitwise xor to a DataBlock and some other type when"
        " bitwise xor cannot be applied to the other type and the elements of the DataBlock.");

    DataBlock<bitwise_xor_t<VT, ST>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] ^ crRHS;
    }

    return result;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
auto operator^(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanBitwiseXor<ST, VT>, "Cannot apply bitwise xor to a DataBlock and some other type when"
        " bitwise xor cannot be applied to the other type and the elements of the DataBlock.");

    DataBlock<bitwise_xor_t<ST, VT>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS ^ crRHS[i];
    }

    return result;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator^=(DB1<VT1, Ds...>&& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<bitwise_xor_t<VT1, VT2>, VT1>, "Cannot apply compound-bitwise-xor-assignment to two DataBlocks when "
        "the bitwise-xor type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] ^= crRHS[i];
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator^=(DB1<VT1, Ds...>& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<bitwise_xor_t<VT1, VT2>, VT1>, "Cannot apply compound-bitwise-xor-assignment to two DataBlocks when "
        "the bitwise-xor type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] ^= crRHS[i];
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator^=(DB<VT, Ds...>&& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<bitwise_xor_t<VT, ST>, VT>, "Cannot apply compound-bitwise-xor-assignment to a DataBlock and another type when "
        "the bitwise-xor type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] ^= crRHS;
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator^=(DB<VT, Ds...>& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<bitwise_xor_t<VT, ST>, VT>, "Cannot apply compound-bitwise-xor-assignment to a DataBlock and another type when "
        "the bitwise-xor type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] ^= crRHS;
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto operator&(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanBitwiseAnd<VT1, VT2>, "Cannot apply bitwise and to DataBlocks for which bitwise and cannot be applied to their element types.");

    DataBlock<bitwise_and_t<VT1, VT2>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] & crRHS[i];
    }

    return result;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
auto operator&(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanBitwiseAnd<VT, ST>, "Cannot apply bitwise and to a DataBlock and some other type when"
        " bitwise and cannot be applied to the other type and the elements of the DataBlock.");

    DataBlock<bitwise_and_t<VT, ST>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] & crRHS;
    }

    return result;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
auto operator&(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanBitwiseAnd<ST, VT>, "Cannot apply bitwise and to a DataBlock and some other type when"
        " bitwise and cannot be applied to the other type and the elements of the DataBlock.");

    DataBlock<bitwise_and_t<ST, VT>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS & crRHS[i];
    }

    return result;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator&=(DB1<VT1, Ds...>&& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<bitwise_and_t<VT1, VT2>, VT1>, "Cannot apply compound-bitwise-and-assignment to two DataBlocks when "
        "the bitwise-and type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] &= crRHS[i];
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator&=(DB1<VT1, Ds...>& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<bitwise_and_t<VT1, VT2>, VT1>, "Cannot apply compound-bitwise-and-assignment to two DataBlocks when "
        "the bitwise-and type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] &= crRHS[i];
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator&=(DB<VT, Ds...>&& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<bitwise_and_t<VT, ST>, VT>, "Cannot apply compound-bitwise-and-assignment to a DataBlock and another type when "
        "the bitwise-and type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] &= crRHS;
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator&=(DB<VT, Ds...>& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<bitwise_and_t<VT, ST>, VT>, "Cannot apply compound-bitwise-and-assignment to a DataBlock and another type when "
        "the bitwise-and type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] &= crRHS;
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto operator|(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanBitwiseOr<VT1, VT2>, "Cannot apply bitwise-or to DataBlocks for which bitwise-or cannot be applied to their element types.");

    DataBlock<bitwise_or_t<VT1, VT2>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] | crRHS[i];
    }

    return result;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
auto operator|(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanBitwiseOr<VT, ST>, "Cannot apply bitwise-or to a DataBlock and some other type when"
        " bitwise-or cannot be applied to the other type and the elements of the DataBlock.");

    DataBlock<bitwise_or_t<VT, ST>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] | crRHS;
    }

    return result;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
auto operator|(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanBitwiseOr<ST, VT>, "Cannot apply bitwise-or to a DataBlock and some other type when"
        " bitwise-or cannot be applied to the other type and the elements of the DataBlock.");

    DataBlock<bitwise_or_t<ST, VT>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS | crRHS[i];
    }

    return result;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator|=(DB1<VT1, Ds...>&& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<bitwise_or_t<VT1, VT2>, VT1>, "Cannot apply compound-bitwise-or-assignment to two DataBlocks when "
        "the bitwise-or type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] |= crRHS[i];
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator|=(DB1<VT1, Ds...>& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<bitwise_or_t<VT1, VT2>, VT1>, "Cannot apply compound-bitwise-or-assignment to two DataBlocks when "
        "the bitwise-or type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] |= crRHS[i];
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator|=(DB<VT, Ds...>&& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<bitwise_or_t<VT, ST>, VT>, "Cannot apply compound-bitwise-or-assignment to a DataBlock and another type when "
        "the bitwise-or type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] |= crRHS;
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator|=(DB<VT, Ds...>& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<bitwise_or_t<VT, ST>, VT>, "Cannot apply compound-bitwise-or-assignment to a DataBlock and another type when "
        "the bitwise-or type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] |= crRHS;
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
bool operator<(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanLessThan<VT1, VT2>, "Cannot apply a less-than comparison to DataBlocks whose element types cannot be tested in less-than comparisons.");
    static_assert(IsExplicitlyConvertible<less_than_t<VT1, VT2>, bool>, "Cannot apply a less-than comparison to DataBlocks whose element types yield a type that "
        "cannot be explicitly converted to bool under less-than comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] < crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
bool operator<(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanLessThan<VT, ST>, "Cannot apply a less-than comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in less-than comparisons.");
    static_assert(IsExplicitlyConvertible<less_than_t<VT, ST>, bool>, "Cannot apply a less-than comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under less-than comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] < crRHS) {}
        else { return false; }
    }

    return true;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
bool operator<(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanLessThan<ST, VT>, "Cannot apply a less-than comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in less-than comparisons.");
    static_assert(IsExplicitlyConvertible<less_than_t<ST, VT>, bool>, "Cannot apply a less-than comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under less-than comparison.");

    for (int i = 0; i < crRHS.nelems; ++i) {
        if (crLHS < crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
bool operator>(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanGreaterThan<VT1, VT2>, "Cannot apply a greater-than comparison to DataBlocks whose element types cannot be tested in greater-than comparisons.");
    static_assert(IsExplicitlyConvertible<greater_than_t<VT1, VT2>, bool>, "Cannot apply a greater-than comparison to DataBlocks whose element types yield a type that "
        "cannot be explicitly converted to bool under greater-than comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] > crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
bool operator>(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanGreaterThan<VT, ST>, "Cannot apply a greater-than comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in greater-than comparisons.");
    static_assert(IsExplicitlyConvertible<greater_than_t<VT, ST>, bool>, "Cannot apply a greater-than comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under greater-than comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] > crRHS) {}
        else { return false; }
    }

    return true;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
bool operator>(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanGreaterThan<ST, VT>, "Cannot apply a greater-than comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in greater-than comparisons.");
    static_assert(IsExplicitlyConvertible<greater_than_t<ST, VT>, bool>, "Cannot apply a greater-than comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under greater-than comparison.");

    for (int i = 0; i < crRHS.nelems; ++i) {
        if (crLHS > crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto operator<<(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanLeftShift<VT1, VT2>, "Cannot left-shift DataBlocks whose element types cannot be left-shifted.");

    DataBlock<left_shift_t<VT1, VT2>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] << crRHS[i];
    }

    return result;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
auto operator<<(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanLeftShift<VT, ST>, "Cannot left-shift a DataBlock and some other type when"
        " the other type and the elements of the DataBlock cannot be left-shifted.");

    DataBlock<left_shift_t<VT, ST>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] << crRHS;
    }

    return result;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
auto operator<<(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanLeftShift<ST, VT>, "Cannot left-shift a DataBlock and some other type when"
        " the other type and the elements of the DataBlock cannot be left-shifted.");

    DataBlock<left_shift_t<ST, VT>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS << crRHS[i];
    }

    return result;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator<<=(DB1<VT1, Ds...>&& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<left_shift_t<VT1, VT2>, VT1>, "Cannot apply compound-left-shift-assignment to two DataBlocks when "
        "the left-shift type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] <<= crRHS[i];
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator<<=(DB1<VT1, Ds...>& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<left_shift_t<VT1, VT2>, VT1>, "Cannot apply compound-left-shift-assignment to two DataBlocks when "
        "the left-shift type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] <<= crRHS[i];
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator<<=(DB<VT, Ds...>&& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<left_shift_t<VT, ST>, VT>, "Cannot apply compound-left-shift-assignment to a DataBlock and another type when "
        "the left-shift type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] <<= crRHS;
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator<<=(DB<VT, Ds...>& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<left_shift_t<VT, ST>, VT>, "Cannot apply compound-left-shift-assignment to a DataBlock and another type when "
        "the left-shift type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] <<= crRHS;
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto operator>>(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanRightShift<VT1, VT2>, "Cannot right-shift DataBlocks whose element types cannot be right-shifted.");

    DataBlock<right_shift_t<VT1, VT2>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] >> crRHS[i];
    }

    return result;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
auto operator>>(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanRightShift<VT, ST>, "Cannot right-shift a DataBlock and some other type when"
        " the other type and the elements of the DataBlock cannot be right-shifted.");

    DataBlock<right_shift_t<VT, ST>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS[i] >> crRHS;
    }

    return result;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
auto operator>>(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanRightShift<ST, VT>, "Cannot right-shift a DataBlock and some other type when"
        " the other type and the elements of the DataBlock cannot be right-shifted.");

    DataBlock<right_shift_t<ST, VT>, Ds...> result;

    for (int i = 0; i < result.nelems; ++i) {
        result[i] = crLHS >> crRHS[i];
    }

    return result;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator>>=(DB1<VT1, Ds...>&& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<right_shift_t<VT1, VT2>, VT1>, "Cannot apply compound-right-shift-assignment to two DataBlocks when "
        "the right-shift type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] >>= crRHS[i];
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
auto& operator>>=(DB1<VT1, Ds...>& rLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(std::convertible_to<right_shift_t<VT1, VT2>, VT1>, "Cannot apply compound-right-shift-assignment to two DataBlocks when "
        "the right-shift type of the element types of the DataBlocks is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] >>= crRHS[i];
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator>>=(DB<VT, Ds...>&& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<right_shift_t<VT, ST>, VT>, "Cannot apply compound-right-shift-assignment to a DataBlock and another type when "
        "the right-shift type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] >>= crRHS;
    }

    return rLHS;
}

template <template <typename, size_t...> typename DB, typename VT, NotADataBlock ST, size_t... Ds>
auto& operator>>=(DB<VT, Ds...>& rLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(std::convertible_to<right_shift_t<VT, ST>, VT>, "Cannot apply compound-right-shift-assignment to a DataBlock and another type when "
        "the right-shift type of the right-hand type and the element type of the DataBlock is not implicitly convertible to the element type of the"
        " assigned-to DataBlock.");

    for (int i = 0; i < rLHS.nelems; ++i) {
        rLHS[i] >>= crRHS;
    }

    return rLHS;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
bool operator==(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanLogicalEquality<VT1, VT2>, "Cannot apply a logical-equality comparison to DataBlocks whose element types cannot be tested in logical-equality comparisons.");
    static_assert(IsExplicitlyConvertible<logical_equality_t<VT1, VT2>, bool>, "Cannot apply a logical-equality comparison to DataBlocks whose element types yield a type that "
        "cannot be explicitly converted to bool under logical-equality comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] == crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
bool operator==(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanLogicalEquality<VT, ST>, "Cannot apply a logical-equality comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in logical-equality comparisons.");
    static_assert(IsExplicitlyConvertible<logical_equality_t<VT, ST>, bool>, "Cannot apply a logical-equality comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under logical-equality comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] == crRHS) {}
        else { return false; }
    }

    return true;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
bool operator==(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanLogicalEquality<ST, VT>, "Cannot apply a logical-equality comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in logical-equality comparisons.");
    static_assert(IsExplicitlyConvertible<logical_equality_t<ST, VT>, bool>, "Cannot apply a logical-equality comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under logical-equality comparison.");

    for (int i = 0; i < crRHS.nelems; ++i) {
        if (crLHS == crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
bool operator!=(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanLogicalInequality<VT1, VT2>, "Cannot apply a logical-inequality comparison to DataBlocks whose element types cannot be tested in logical-inequality comparisons.");
    static_assert(IsExplicitlyConvertible<logical_inequality_t<VT1, VT2>, bool>, "Cannot apply a logical-inequality comparison to DataBlocks whose element types yield a type that "
        "cannot be explicitly converted to bool under logical-inequality comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] != crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
bool operator!=(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanLogicalInequality<VT, ST>, "Cannot apply a logical-inequality comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in logical-inequality comparisons.");
    static_assert(IsExplicitlyConvertible<logical_inequality_t<VT, ST>, bool>, "Cannot apply a logical-inequality comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under logical-inequality comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] != crRHS) {}
        else { return false; }
    }

    return true;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
bool operator!=(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanLogicalInequality<ST, VT>, "Cannot apply a logical-inequality comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in logical-inequality comparisons.");
    static_assert(IsExplicitlyConvertible<logical_inequality_t<ST, VT>, bool>, "Cannot apply a logical-inequality comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under logical-inequality comparison.");

    for (int i = 0; i < crRHS.nelems; ++i) {
        if (crLHS != crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
bool operator<=(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanLessThanOrEqual<VT1, VT2>, "Cannot apply a less-than-or-equal comparison to DataBlocks whose element types cannot be tested in less-than-or-equal comparisons.");
    static_assert(IsExplicitlyConvertible<less_than_or_equal_t<VT1, VT2>, bool>, "Cannot apply a less-than-or-equal comparison to DataBlocks whose element types yield a type that "
        "cannot be explicitly converted to bool under less-than-or-equal comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] <= crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
bool operator<=(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanLessThanOrEqual<VT, ST>, "Cannot apply a less-than-or-equal comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in less-than-or-equal comparisons.");
    static_assert(IsExplicitlyConvertible<less_than_or_equal_t<VT, ST>, bool>, "Cannot apply a less-than-or-equal comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under less-than-or-equal comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] <= crRHS) {}
        else { return false; }
    }

    return true;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
bool operator<=(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanLessThanOrEqual<ST, VT>, "Cannot apply a less-than-or-equal comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in less-than-or-equal comparisons.");
    static_assert(IsExplicitlyConvertible<less_than_or_equal_t<ST, VT>, bool>, "Cannot apply a less-than-or-equal comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under less-than-or-equal comparison.");

    for (int i = 0; i < crRHS.nelems; ++i) {
        if (crLHS <= crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
bool operator>=(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanGreaterThanOrEqual<VT1, VT2>, "Cannot apply a greater-than-or-equal comparison to DataBlocks whose element types cannot be tested in greater-than-or-equal comparisons.");
    static_assert(IsExplicitlyConvertible<greater_than_or_equal_t<VT1, VT2>, bool>, "Cannot apply a greater-than-or-equal comparison to DataBlocks whose element types yield a type that "
        "cannot be explicitly converted to bool under greater-than-or-equal comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] >= crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
bool operator>=(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanGreaterThanOrEqual<VT, ST>, "Cannot apply a greater-than-or-equal comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in greater-than-or-equal comparisons.");
    static_assert(IsExplicitlyConvertible<greater_than_or_equal_t<VT, ST>, bool>, "Cannot apply a greater-than-or-equal comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under greater-than-or-equal comparison.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        if (crLHS[i] >= crRHS) {}
        else { return false; }
    }

    return true;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
bool operator>=(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanGreaterThanOrEqual<ST, VT>, "Cannot apply a greater-than-or-equal comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in greater-than-or-equal comparisons.");
    static_assert(IsExplicitlyConvertible<greater_than_or_equal_t<ST, VT>, bool>, "Cannot apply a greater-than-or-equal comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that cannot be explicitly converted"
        " to bool under greater-than-or-equal comparison.");

    for (int i = 0; i < crRHS.nelems; ++i) {
        if (crLHS >= crRHS[i]) {}
        else { return false; }
    }

    return true;
}

template <typename VT1, typename VT2, template <typename, size_t...> typename DB1,
    template <typename, size_t...> typename DB2, size_t... Ds>
three_way_compare_t<VT1, VT2> operator<=>(DB1<VT1, Ds...> const& crLHS, DB2<VT2, Ds...> const& crRHS)
    requires (ADataBlock<DB1<VT1, Ds...>>&& ADataBlock<DB2<VT2, Ds...>>) {

    static_assert(CanThreeWayCompare<VT1, VT2>, "Cannot apply a three-way comparison to DataBlocks whose element types cannot be tested in three-way comparisons.");
    static_assert(IsValidThreeWayCompareReturn<three_way_compare_t<VT1, VT2>>, "Cannot apply a three-way comparison to DataBlocks whose element types yield a type that "
        "is not a canonical three-way comparison return type.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        three_way_compare_t<VT1, VT2> result = crLHS[i] <=> crRHS[i];
        if (result > 0) return result;
        if (result < 0) return result;
    }

    return crLHS[crLHS.nelems - 1] <=> crRHS[crRHS.nelems - 1];
}

template <typename VT, NotADataBlock ST, template <typename, size_t...> typename DB, size_t... Ds>
three_way_compare_t<VT, ST> operator<=>(DB<VT, Ds...> const& crLHS, ST const& crRHS)
    requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanThreeWayCompare<VT, ST>, "Cannot apply a three-way comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in three-way comparisons.");
    static_assert(IsValidThreeWayCompareReturn<three_way_compare_t<VT, ST>>, "Cannot apply a three-way comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that is not a canonical "
        "three-way comparison return type.");

    for (int i = 0; i < crLHS.nelems; ++i) {
        three_way_compare_t<VT, ST> result = crLHS[i] <=> crRHS;
        if (result > 0) return result;
        if (result < 0) return result;
    }

    return crLHS[crLHS.nelems - 1] <=> crRHS;
}

template <NotADataBlock ST, typename VT, template <typename, size_t...> typename DB, size_t... Ds>
three_way_compare_t<ST, VT> operator<=>(ST const& crLHS, DB<VT, Ds...> const& crRHS) requires (ADataBlock<DB<VT, Ds...>>) {

    static_assert(CanThreeWayCompare<ST, VT>, "Cannot apply a three-way comparison to a DataBlock and some other type when"
        " the other type and the element type of the DataBlock cannot be tested in three-way comparisons.");
    static_assert(IsValidThreeWayCompareReturn<three_way_compare_t<ST, VT>>, "Cannot apply a three-way comparison to a DataBlock"
        " and some other type when the other type and the element type of the DataBlock yield a type that is not a canonical three-way"
        " comparison return type.");

    for (int i = 0; i < crRHS.nelems; ++i) {
        three_way_compare_t<ST, VT> result = crLHS <=> crRHS[i];
        if (result > 0) return result;
        if (result < 0) return result;
    }

    return crLHS <=> crRHS[crRHS.nelems - 1];
}

template <typename T>
concept Printable = requires (T a, std::ostream& rOS) {
    rOS << a;
};

size_t PrintLen(Printable auto const& crData) {
    std::ostringstream stream;

    stream << crData;

    std::string str = stream.str();
    return str.size();
}

template <typename T>
class Width {

private:
    T const& data;
    size_t const width;

public:

    Width(T const& crData, size_t width) : data{ crData }, width{ width } { }

private:
    friend std::ostream& operator<<(std::ostream& rOS, Width const& crWidth) {
        std::ostringstream stream;

        stream << crWidth.data;

        std::string str = stream.str();
        const char* cpStr = str.c_str();

        size_t len = strlen(cpStr);

        if (len > crWidth.width) {
            for (int i = 0; i < crWidth.width; ++i) {
                rOS << cpStr[i];
            }
        }

        else {
            int i;

            for (i = 0; i < len; ++i) {
                rOS << cpStr[i];
            }

            for (; i < crWidth.width; ++i) {
                rOS << ' ';
            }
        }

        return rOS;
    }
};

template <typename C, typename RETURN, typename... Ts, typename... As>
auto PartialBind(RETURN(C::* pMember)(Ts...), C* pInstance, As... bound_args) {
    return [=](auto&&... args) {return (pInstance->*pMember)(bound_args..., std::forward<decltype(args)>(args)...);};
}

template <typename F, typename... As>
auto PartialBind(F& rFunc, As... bound_args) {
    return [=](auto&&... args) {return rFunc(bound_args..., std::forward<decltype(args)>(args)...);};
}

template <typename DB>
std::string ToString(DB const& crDB, size_t width = -1) {

    if (width == -1) {
        size_t len = 0, max_len = 0;
        for (int i = 0; i < crDB.NElems(); ++i) {
            len = PrintLen(crDB[i]);

            if (len > max_len) max_len = len;
        }

        width = max_len;
    }

    std::ostringstream stream;

    size_t row_len = crDB.Shape(crDB.NDims() - 1);

    for (int flat_idx = 0; flat_idx < crDB.NElems(); flat_idx += row_len) {
        stream << "[ ";
        for (int row_idx = 0; row_idx < row_len; ++row_idx) {
            stream << Width(crDB[flat_idx + row_idx], width) << ' ';
        }
        stream << ']';

        size_t new_flat_idx = flat_idx + row_len;
        size_t boundary_idx = 1;

        for (int dim_idx = crDB.NDims() - 1; dim_idx > 0; --dim_idx) {
            boundary_idx *= crDB.Shape(dim_idx);

            if (new_flat_idx % boundary_idx == 0) stream << '\n';
        }
    }

    return stream.str();
}

template <typename T, size_t... Ds>
std::string ToString(DataBlockView<T, Ds...> const& crDBV, size_t width = 3) {
    std::ostringstream stream;

    stream << "[ ";
    size_t nelems = crDBV.Size();

    for (size_t i = 0; i < nelems; ++i) {
        stream << Width(crDBV[i], width) << ' ';
    }

    stream << ']';

    return stream.str();
}

template <typename T, size_t... Ds>
T const* ValuePointer(DataBlock<T, Ds...>&& crDB) = delete;

template <typename T, size_t... Ds>
T const* ValuePointer(DataBlock<T, Ds...> const& crDB) {
    return &crDB[0];
}

void foo(ADataBlock auto const&) { }

#endif