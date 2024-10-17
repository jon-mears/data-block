#ifndef INDEX_RANGE_HPP
#define INDEX_RANGE_HPP

struct Index;
struct IndexRange;
struct SkipIndexRange;

struct SkipIndexRange {
    size_t start;
    size_t end;
    size_t skip;
};

struct IndexRange {
    size_t start;
    size_t end;

    SkipIndexRange operator|(Index rhs);
};

struct Index {
    size_t index;

    IndexRange operator-(Index rhs) {
        return { .start = index, .end = rhs.index };
    }
};

SkipIndexRange IndexRange::operator|(Index rhs) {
    return { .start = start, .end = end, .skip = rhs.index };
}

template <int IDX>
struct TemplatedIndex { };

template <int START, int STOP, int STEP>
inline constexpr bool ValidIndexRangeV = (
    START >= 0 &&
    STOP >= -1 &&
    STEP != 0 &&
    (START < STOP && STEP > 0) || (START > STOP && STEP < 0)
    );

template <int START, int STOP, int STEP = 1>
struct TemplatedIndexRange { };

template <int IDX>
TemplatedIndex<-IDX> operator-(TemplatedIndex<IDX> const& crIndex) {
    return { };
}

template <int START, int STOP>
TemplatedIndexRange<START, STOP> operator>>(TemplatedIndex<START> const& crLHS, TemplatedIndex<STOP> const& crRHS) {
    return { };
}

template <int START, int STOP, int STEP>
TemplatedIndexRange<START, STOP, STEP> operator|(TemplatedIndexRange<START, STOP> const& crLHS, TemplatedIndex<STEP> const& crRHS) {
    return { };
}

template <int V>
constexpr int make_index() {
    return V;
}

template <int V, char C, char... Cs>
constexpr int make_index() {
    return make_index<V * 10 + C - '0', Cs...>();
}

template <char... Cs>
constexpr TemplatedIndex<make_index<0, Cs...>()> operator""_idx() {
    return { };
}

#endif