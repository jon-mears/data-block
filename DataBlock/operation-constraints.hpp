#ifndef OPERATION_CONSTRAINTS_HPP
#define OPERATION_CONSTRAINTS_HPP

#include <concepts>

template <typename T>
struct Type { };

template <typename T>
struct ValueType : Type<T> {
    template <typename U>
    constexpr ValueType(U) : Type<T>() { }
};

template <typename U>
ValueType(U) -> ValueType<typename U::value_type>;

template <typename T, T... Ts>
struct HomogeneousPack { };

template <size_t... Ds>
struct Dimensions : HomogeneousPack<size_t, Ds...> {
    template <typename T, template <typename, size_t...> typename DB>
    constexpr Dimensions(DB<T, Ds...>) : HomogeneousPack<size_t, Ds...>() { }
};

template <typename T, size_t... Ds, template <typename, size_t...> typename DB>
Dimensions(DB<T, Ds...>) -> Dimensions<Ds...>;

namespace OpTraits {
    template <typename T1, typename T2>
    concept CanAssignC = requires (T1 a, T2 b) {
        a = b;
    };

    template <typename T1, typename T2>
    constexpr bool CanAssignV = CanAssignC<T1, T2>;

    template <typename T1, typename T2>
    constexpr bool CanAssign(Type<T1>, Type<T2>) {
        return CanAssignV<T1, T2>;
    }

    template <typename T1, typename T2>
    constexpr bool CanAssign(T1, T2) {
        return CanAssignC<T1, T2>;
    }

    template <typename T1, typename T2> requires (CanAssignC<T1, T2>)
    struct AssignTraits {
        typedef decltype(std::declval<T1>() = std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using AssignT = AssignTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    constexpr Type<AssignT<T1, T2>> AssignType(Type<T1>, Type<T2>) {
        return Type<AssignT<T1, T2>>();
    }

    template <typename T1, typename T2>
    constexpr Type<AssignT<T1, T2>> AssignType(T1, T2) {
        return Type<AssignT<T1, T2>>();
    }

    template <typename T>
    concept CanArrowC = requires (T a) {
        a.operator->();
    };

    template <typename T>
    constexpr bool CanArrowV = CanArrowC<T>;

    template <typename T>
    constexpr bool CanArrow(Type<T>) {
        return CanArrowC<T>;
    } 

    template <typename T>
    constexpr bool CanArrow(T) {
        return CanArrowC<T>;
    }

    template <typename T> requires (CanArrowC<T>)
    struct ArrowTraits {
        typedef decltype(std::declval<T>().operator->()) result_type;
    };

    template <typename T>
    using ArrowT = ArrowTraits<T>::result_type;

    template <typename T>
    constexpr Type<ArrowT<T>> ArrowType(Type<T>) {
        return Type<ArrowT<T>>();
    }

    template <typename T>
    constexpr Type<ArrowT<T>> ArrowType(T) {
        return Type<ArrowT<T>>();
    }

    template <typename T>
    concept CanonicalThreeWayComparisonReturnC = requires (T a) {
        { a < 0 }  -> std::convertible_to<bool>;
        { a > 0 }  -> std::convertible_to<bool>;
        { a == 0 } -> std::convertible_to<bool>;
    };

    template <typename T>
    concept CanUnaryPlusC = requires (T a) {
        +a;
    };

    template <typename T>
    constexpr bool CanUnaryPlusV = CanUnaryPlusC<T>;

    template <typename T>
    constexpr bool CanUnaryPlus(Type<T>) {
        return CanUnaryPlusV<T>;
    }

    template <typename T>
    constexpr bool CanUnaryPlus(T) {
        return CanUnaryPlusV<T>;
    }

    template <typename T> requires (CanUnaryPlusC<T>)
    struct UnaryPlusTraits {
        typedef decltype(+std::declval<T>()) result_type;
    };

    template <typename T>
    using UnaryPlusT = UnaryPlusTraits<T>::result_type;

    template <typename T>
    constexpr Type<UnaryPlusT<T>> UnaryPlusType(Type<T>) {
        return Type<UnaryPlusT<T>>();
    }

    template <typename T>
    constexpr Type<UnaryPlusT<T>> UnaryPlusType(T) {
        return Type<UnaryPlusT<T>>();
    }

    template <typename T>
    concept CanDereferenceC = requires (T a) {
        *a;
    };

    template <typename T>
    constexpr bool CanDereferenceV = CanDereferenceC<T>;

    template <typename T> requires (CanDereferenceC<T>)
    struct DereferenceTraits {
        typedef decltype(*std::declval<T>()) result_type;
    };

    template <typename T>
    using DereferenceT = DereferenceTraits<T>::result_type;

    template <typename T>
    concept CanAddressOfC = requires (T a) {
        &a;
    };

    template <typename T> requires (CanAddressOfC<T>)
        struct AddressOfTraits {
        typedef decltype(&std::declval<T>()) result_type;
    };

    template <typename T>
    using AddressOfT = AddressOfTraits<T>::result_type;

    template <typename T>
    concept CanPrefixIncrementC = requires (T a) {
        ++a;
    };

    template <typename T> requires (CanPrefixIncrementC<T>)
    struct PrefixIncrementTraits {
        typedef decltype(++std::declval<T>()) result_type;
    };

    template <typename T>
    using PrefixIncrementT = PrefixIncrementTraits<T>::result_type;

    template <typename T>
    concept CanPostfixIncrementC = requires (T a) {
        a++;
    };

    template <typename T> requires (CanPostfixIncrementC<T>)
    struct PostfixIncrementTraits {
        typedef decltype(std::declval<T>()++) result_type;
    };

    template <typename T>
    using PostfixIncrementT = PostfixIncrementTraits<T>::result_type;

    template <typename T>
    concept CanPrefixDecrementC = requires (T a) {
        --a;
    };

    template <typename T> requires (CanPrefixDecrementC<T>)
    struct PrefixDecrementTraits {
        typedef decltype(--std::declval<T>()) result_type;
    };

    template <typename T>
    using PrefixDecrementT = PrefixDecrementTraits<T>::result_type;

    template <typename T>
    concept CanPostfixDecrementC = requires (T a) {
        a--;
    };

    template <typename T> requires (CanPostfixDecrementC<T>)
    struct PostfixDecrementTraits {
        typedef decltype(std::declval<T>()--) result_type;
    };

    template <typename T>
    using PostfixDecrementT = PostfixDecrementTraits<T>::result_type;

    template <typename T1, typename T2>
    concept CanThreeWayCompareC = requires (T1 a, T2 b) {
        a <=> b;
    };

    template <typename T1, typename T2> requires (CanThreeWayCompareC<T1, T2>)
    struct ThreeWayCompareTraits {
        typedef decltype(std::declval<T1>() <=> std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using ThreeWayCompareT = ThreeWayCompareTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanSumC = requires (T1 a, T2 b) {
        a + b;
    };

    template <typename T1, typename T2> requires (CanSumC<T1, T2>)
    struct SumTraits {
        typedef decltype(std::declval<T1>() + std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using SumT = SumTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanSumAssignC = requires (T1 a, T2 b) {
        a += b;
    };

    template <typename T1, typename T2> requires (CanSumAssignC<T1, T2>)
    struct SumAssignTraits {
        typedef decltype(std::declval<T1>() += std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using SumAssignT = SumAssignTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanSubtractC = requires (T1 a, T2 b) {
        a - b;
    };

    template <typename T1, typename T2> requires (CanSubtractC<T1, T2>)
    struct SubtractTraits {
        typedef decltype(std::declval<T1>() - std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using SubtractT = SubtractTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanSubtractAssignC = requires (T1 a, T2 b) {
        a -= b;
    };

    template <typename T1, typename T2> requires (CanSubtractAssignC<T1, T2>)
    struct SubtractAssignTraits {
        typedef decltype(std::declval<T1>() -= std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using SubtractAssignT = SubtractAssignTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanMultiplyC = requires (T1 a, T2 b) {
        a* b;
    };

    template <typename T1, typename T2> requires (CanMultiplyC<T1, T2>)
    struct MultiplyTraits {
        typedef decltype(std::declval<T1>()* std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using MultiplyT = MultiplyTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanMultiplyAssignC = requires (T1 a, T2 b) {
        a *= b;
    };

    template <typename T1, typename T2> requires (CanMultiplyAssignC<T1, T2>)
    struct MultiplyAssignTraits {
        typedef decltype(std::declval<T1>() *= std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using MultiplyAssignT = MultiplyAssignTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanDivideC = requires (T1 a, T2 b) {
        a / b;
    };

    template <typename T1, typename T2> requires (CanDivideC<T1, T2>)
    struct DivideTraits {
        typedef decltype(std::declval<T1>() / std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using DivideT = DivideTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanDivideAssignC = requires (T1 a, T2 b) {
        a /= b;
    };

    template <typename T1, typename T2> requires (CanDivideAssignC<T1, T2>)
    struct DivideAssignTraits {
        typedef decltype(std::declval<T1>() /= std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using DivideAssignT = DivideAssignTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanRemainderC = requires (T1 a, T2 b) {
        a % b;
    };

    template <typename T1, typename T2> requires (CanRemainderC<T1, T2>)
    struct RemainderTraits {
        typedef decltype(std::declval<T1>() % std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using RemainderT = RemainderTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanRemainderAssignC = requires (T1 a, T2 b) {
        a %= b;
    };

    template <typename T1, typename T2> requires (CanRemainderAssignC<T1, T2>)
    struct RemainderAssignTraits {
        typedef decltype(std::declval<T1>() %= std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using RemainderAssignT = RemainderAssignTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanBitwiseAndC = requires (T1 a, T2 b) {
        a& b;
    };

    template <typename T1, typename T2> requires (CanBitwiseAndC<T1, T2>)
    struct BitwiseAndTraits {
        typedef decltype(std::declval<T1>()& std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using BitwiseAndT = BitwiseAndTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanBitwiseAndAssignC = requires (T1 a, T2 b) {
        a &= b;
    };

    template <typename T1, typename T2> requires (CanBitwiseAndAssignC<T1, T2>)
    struct BitwiseAndAssignTraits {
        typedef decltype(std::declval<T1>() &= std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using BitwiseAndAssignT = BitwiseAndAssignTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanBitwiseOrC = requires (T1 a, T2 b) {
        a | b;
    };

    template <typename T1, typename T2> requires (CanBitwiseOrC<T1, T2>)
    struct BitwiseOrTraits {
        typedef decltype(std::declval<T1>() | std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using BitwiseOrT = BitwiseOrTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanBitwiseOrAssignC = requires (T1 a, T2 b) {
        a |= b;
    };

    template <typename T1, typename T2> requires (CanBitwiseOrAssignC<T1, T2>)
    struct BitwiseOrAssignTraits {
        typedef decltype(std::declval<T1>() |= std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using BitwiseOrAssignT = BitwiseOrAssignTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanBitwiseLeftShiftC = requires (T1 a, T2 b) {
        a << b;
    };

    template <typename T1, typename T2> requires (CanBitwiseLeftShiftC<T1, T2>)
    struct BitwiseLeftShiftTraits {
        typedef decltype(std::declval<T1>() << std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using BitwiseLeftShiftT = BitwiseLeftShiftTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanBitwiseLeftShiftAssignC = requires (T1 a, T2 b) {
        a <<= b;
    };

    template <typename T1, typename T2> requires (CanBitwiseLeftShiftAssignC<T1, T2>)
    struct BitwiseLeftShiftAssignTraits {
        typedef decltype(std::declval<T1>() <<= std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using BitwiseLeftShiftAssignT = BitwiseLeftShiftAssignTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanBitwiseRightShiftC = requires (T1 a, T2 b) {
        a >> b;
    };

    template <typename T1, typename T2> requires (CanBitwiseRightShiftC<T1, T2>)
    struct BitwiseRightShiftTraits {
        typedef decltype(std::declval<T1>() >> std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using BitwiseRightShiftT = BitwiseRightShiftTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanBitwiseRightShiftAssignC = requires (T1 a, T2 b) {
        a >>= b;
    };

    template <typename T1, typename T2> requires (CanBitwiseRightShiftAssignC<T1, T2>)
    struct BitwiseRightShiftAssignTraits {
        typedef decltype(std::declval<T1>() >>= std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using BitwiseRightShiftAssignT = BitwiseRightShiftAssignTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanLogicalAndC = requires (T1 a, T2 b) {
        a && b;
    };

    template <typename T1, typename T2> requires (CanLogicalAndC<T1, T2>)
    struct LogicalAndTraits {
        typedef decltype(std::declval<T1>() && std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using LogicalAndT = LogicalAndTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanLogicalOrC = requires (T1 a, T2 b) {
        a || b;
    };

    template <typename T1, typename T2> requires (CanLogicalOrC<T1, T2>)
    struct LogicalOrTraits {
        typedef decltype(std::declval<T1>() || std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using LogicalOrT = LogicalOrTraits<T1, T2>::result_type;

    template <typename T>
    concept CanUnaryMinusC = requires (T a) {
        -a;
    };

    template <typename T> requires (CanUnaryMinusC<T>)
    struct UnaryMinusTraits {
        typedef decltype(-std::declval<T>()) result_type;
    };

    template <typename T>
    using UnaryMinusT = UnaryMinusTraits<T>::result_type;

    template <typename T>
    concept CanLogicalNotC = requires (T a) {
        !a;
    };

    template <typename T> requires (CanLogicalNotC<T>)
    struct LogicalNotTraits {
        typedef decltype(!std::declval<T>()) result_type;
    };

    template <typename T>
    using LogicalNotT = LogicalNotTraits<T>::result_type;

    template <typename T>
    concept CanBitwiseNotC = requires (T a) {
        ~a;
    };

    template <typename T> requires (CanBitwiseNotC<T>)
    struct BitwiseNotTraits {
        typedef decltype(~std::declval<T>()) result_type;
    };

    template <typename T>
    using BitwiseNotT = BitwiseNotTraits<T>::result_type;

    template <typename F, typename... ARGS>
    concept CanFunctionCallC = requires (F function, ARGS... arguments) {
        function(arguments...);
    };

    template <typename F, typename... ARGS> requires (CanFunctionCallC<F, ARGS...>)
    struct FunctionCallTraits {
        typedef decltype(std::declval<F>()(std::declval<ARGS>()...)) result_type;
    };

    template <typename F, typename... ARGS>
    using FunctionCallT = FunctionCallTraits<F, ARGS...>::result_type;

    template <typename F>
    concept CanCanonicalSubscriptC = requires (F a, size_t idx) {
        a[idx];
    };

    template <typename F> requires (CanCanonicalSubscriptC<F>)
    struct CanonicalSubscriptTraits {
        typedef decltype(std::declval<F>()[std::declval<size_t>()]) result_type;
    };

    template <typename F>
    using CanonicalSubscriptT = CanonicalSubscriptTraits<F>::result_type;

    template <typename T1, typename T2>
    concept CanLogicalLessThanC = requires (T1 a, T2 b) {
        a < b;
    };

    template <typename T1, typename T2> requires (CanLogicalLessThanC<T1, T2>)
    struct less_than {
        typedef decltype(std::declval<T1>() < std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using LogicalLessThanT = less_than<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanLogicalGreaterThanC = requires (T1 a, T2 b) {
        a > b;
    };

    template <typename T1, typename T2> requires (CanLogicalGreaterThanC<T1, T2>)
    struct LogicalGreaterThanTraits {
        typedef decltype(std::declval<T1>() > std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using LogicalGreaterThanT = LogicalGreaterThanTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanBitwiseXorC = requires (T1 a, T2 b) {
        a ^ b;
    };

    template <typename T1, typename T2> requires (CanBitwiseXorC<T1, T2>)
    struct BitwiseXorTraits {
        typedef decltype(std::declval<T1>() ^ std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using BitwiseXorT = BitwiseXorTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanBitwiseXorAssignC = requires (T1 a, T2 b) {
        a ^= b;
    };

    template <typename T1, typename T2> requires (CanBitwiseXorAssignC<T1, T2>)
    struct BitwiseXorTraits {
        typedef decltype(std::declval<T1>() ^= std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using BitwiseXorAssignT = BitwiseXorAssignTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanLogicalEqualsC = requires (T1 a, T2 b) {
        a == b;
    };

    template <typename T1, typename T2> requires (CanLogicalEqualsC<T1, T2>)
    struct LogicalEqualsTraits {
        typedef decltype(std::declval<T1>() == std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using LogicalEqualsT = LogicalEqualsTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanLogicalNotEqualsC = requires (T1 a, T2 b) {
        a != b;
    };

    template <typename T1, typename T2> requires (CanLogicalNotEqualsC<T1, T2>)
    struct LogicalNotEqualsTraits {
        typedef decltype(std::declval<T1>() != std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using LogicalNotEqualsT = LogicalNotEqualsTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanLogicalLessThanOrEqualsC = requires (T1 a, T2 b) {
        a <= b;
    };

    template <typename T1, typename T2> requires (CanLogicalLessThanOrEqualsC<T1, T2>)
    struct LogicalLessThanOrEqualsTraits {
        typedef decltype(std::declval<T1>() <= std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using LogicalLessThanOrEqualsT = LogicalLessThanOrEqualsTraits<T1, T2>::result_type;

    template <typename T1, typename T2>
    concept CanLogicalGreaterThanOrEqualsC = requires (T1 a, T2 b) {
        a >= b;
    };

    template <typename T1, typename T2> requires (CanLogicalGreaterThanOrEqualsC<T1, T2>)
    struct LogicalGreaterThanOrEqualsTraits {
        typedef decltype(std::declval<T1>() >= std::declval<T2>()) result_type;
    };

    template <typename T1, typename T2>
    using LogicalGreaterThanOrEqualsT = LogicalGreaterThanOrEqualsTraits<T1, T2>::result_type;
}

#endif