#include "lib/posit.h"
#include "lib/util.h"
#include "lib/pack.h"
#include "lib/op1.h"
#include "lib/op2.h"
#include "lib/posit_types.h"
#include <cstdio>
#include <cmath>

using namespace std;

Posit::Posit(POSIT_UTYPE bits, int nbits, int es) :
    mBits(bits),
    mNbits(nbits),
    mEs(es)
{
}

Posit::Posit(int nbits, int es) :
    Posit(POSIT_ZERO, nbits, es)
{
}

bool Posit::isZero() const
{
    return util_is_zero(mBits);
}

bool Posit::isNar() const
{
    return util_is_nar(mBits);
}

bool Posit::isNeg() const
{
    return util_is_neg(mBits);
}

int Posit::nbits() const
{
    return mNbits;
}

int Posit::ss() const
{
    return util_ss();
}

int Posit::rs() const
{
    return util_rs(mBits, mNbits);
}

int Posit::es() const
{
    return util_es(mBits, mNbits, mEs);
}

int Posit::fs() const
{
    return util_fs(mBits, mNbits, mEs);
}

Posit Posit::zero() const
{
    return Posit(POSIT_ZERO, mNbits, mEs);
}

Posit Posit::one() const
{
    return Posit(POSIT_ONE, mNbits, mEs);
}

Posit Posit::nar() const
{
    return Posit(POSIT_NAR, mNbits, mEs);
}

Posit Posit::neg() const
{
    if (isNar()) {
        return nar();
    }

    return Posit(util_neg(mBits, mNbits), mNbits, mEs);
}

Posit Posit::rec() const
{
    if (isNar() || isZero()) {
        return nar();
    }

    return one().div(*this);
}

Posit Posit::sqrt() const
{
    if (isNar() || isNeg()) {
        return nar();
    } else if (isZero()) {
        return zero();
    }

    unpacked_t aup = unpack_posit(mBits, mNbits, mEs);
    unpacked_t up = op1_sqrt(aup);

    return Posit(pack_posit(up, mNbits, mEs), mNbits, mEs);
}

Posit Posit::add(const Posit& p) const
{
    if (isNar() || p.isNar()) {
        return nar();
    } else if (isZero()) {
        return p;
    } else if (p.isZero()) {
        return *this;
    } else if (neg().eq(p)) {
        return zero();
    }

    unpacked_t aup = unpack_posit(mBits, mNbits, mEs);
    unpacked_t bup = unpack_posit(p.mBits, p.mNbits, p.mEs);
    unpacked_t up = op2_add(aup, bup);

    return Posit(pack_posit(up, mNbits, mEs), mNbits, mEs);
}

Posit Posit::sub(const Posit& p) const
{
    if (isNar() || p.isNar()) {
        return nar();
    } else if (isZero()) {
        return p.neg();
    } else if (p.isZero()) {
        return *this;
    } else if (eq(p)) {
        return zero();
    }

    unpacked_t aup = unpack_posit(mBits, mNbits, mEs);
    unpacked_t bup = unpack_posit(p.mBits, p.mNbits, p.mEs);
    unpacked_t up = op2_sub(aup, bup);

    return Posit(pack_posit(up, mNbits, mEs), mNbits, mEs);
}

Posit Posit::mul(const Posit& p) const
{
    if (isNar() || p.isNar()) {
        return nar();
    } else if (isZero() || p.isZero()) {
        return zero();
    }

    unpacked_t aup = unpack_posit(mBits, mNbits, mEs);
    unpacked_t bup = unpack_posit(p.mBits, p.mNbits, p.mEs);
    unpacked_t up = op2_mul(aup, bup);

    return Posit(pack_posit(up, mNbits, mEs), mNbits, mEs);
}

Posit Posit::div(const Posit& p) const
{
    if (isNar() || p.isNar() || p.isZero()) {
        return nar();
    } else if (isZero()) {
        return zero();
    }

    unpacked_t aup = unpack_posit(mBits, mNbits, mEs);
    unpacked_t bup = unpack_posit(p.mBits, p.mNbits, p.mEs);
    unpacked_t up = op2_div(aup, bup);

    return Posit(pack_posit(up, mNbits, mEs), mNbits, mEs);
}

bool Posit::eq(const Posit& p) const
{
    return mBits == p.mBits;
}

bool Posit::gt(const Posit& p) const
{
    return (POSIT_STYPE)mBits > (POSIT_STYPE)p.mBits;
}

bool Posit::ge(const Posit& p) const
{
    return gt(p) || eq(p);
}

bool Posit::lt(const Posit& p) const
{
    return !gt(p) && !eq(p);
}

bool Posit::le(const Posit& p) const
{
    return !gt(p);
}

void Posit::set(Posit p)
{
    mBits = pack_posit(unpack_posit(p.mBits, p.mNbits, p.mEs), mNbits, mEs);
}

void Posit::set(float n)
{
    switch (fpclassify(n)) {
    case FP_INFINITE:
    case FP_NAN:
        mBits = POSIT_NAR;
        break;
    case FP_ZERO:
        mBits = POSIT_ZERO;
        break;
    default:
        mBits = pack_posit(unpack_float(n), mNbits, mEs);
        break;
    }
}

void Posit::set(double n)
{
    switch (fpclassify(n)) {
    case FP_INFINITE:
    case FP_NAN:
        mBits = POSIT_NAR;
        break;
    case FP_ZERO:
        mBits = POSIT_ZERO;
        break;
    default:
        mBits = pack_posit(unpack_double(n), mNbits, mEs);
        break;
    }
}

float Posit::getFloat() const
{
    if (isZero()) {
        return 0.f;
    } else if (isNar()) {
        return 0.f / 0.f;
    }

    return pack_float(unpack_posit(mBits, mNbits, mEs));
}

double Posit::getDouble() const
{
    if (isZero()) {
        return 0.0;
    } else if (isNar()) {
        return 0.0 / 0.0;
    }

    return pack_double(unpack_posit(mBits, mNbits, mEs));
}

void Posit::setBits(POSIT_UTYPE bits)
{
    mBits = LSHIFT(bits, POSIT_WIDTH - mNbits);
}

POSIT_UTYPE Posit::getBits()
{
    return RSHIFT(mBits, POSIT_WIDTH - mNbits);
}

void Posit::print()
{
    Posit p = (isNeg() ? neg() : *this);

    printf("{%d, %d} ", mNbits, mEs);

    if (isNar()) {
        printf("NaR\n");
        return;
    }

    for (int i = POSIT_WIDTH - 1; i >= POSIT_WIDTH - mNbits; i--) {
        printf("%d", RSHIFT(mBits, i) & 1);
    }

    printf(" -> ");
    printf(isNeg() ? "-" : "+");

    for (int i = POSIT_WIDTH - ss() - 1; i >= POSIT_WIDTH - mNbits; i--) {
        printf("%d", RSHIFT(p.mBits, i) & 1);

        if (i != POSIT_WIDTH - mNbits &&
            ((i == POSIT_WIDTH - ss() - p.rs()) ||
             (i == POSIT_WIDTH - ss() - p.rs() - mEs))) {
            printf(" ");
        }
    }

    printf(" = %lg\n", getDouble());
}

Posit8::Posit8() :
    Posit(8, 0)
{

}

Posit8::Posit8(Posit v) :
    Posit8()
{
    set(v);
}

Posit8::Posit8(float v) :
    Posit8()
{
    set(v);
}

Posit8::Posit8(double v) :
    Posit8()
{
    set(v);
}

Posit16::Posit16() :
    Posit(16, 1)
{

}

Posit16::Posit16(Posit v) :
    Posit16()
{
    set(v);
}

Posit16::Posit16(float v) :
    Posit16()
{
    set(v);
}

Posit16::Posit16(double v) :
    Posit16()
{
    set(v);
}

Posit32::Posit32() :
    Posit(32, 2)
{

}

Posit32::Posit32(Posit v) :
    Posit32()
{
    set(v);
}

Posit32::Posit32(float v) :
    Posit32()
{
    set(v);
}

Posit32::Posit32(double v) :
    Posit32()
{
    set(v);
}

Posit operator+(const Posit& a, const Posit& b)
{
    return a.add(b);
}

Posit operator-(const Posit& a, const Posit& b)
{
    return a.sub(b);
}

Posit operator*(const Posit& a, const Posit& b)
{
    return a.mul(b);
}

Posit operator/(const Posit& a, const Posit& b)
{
    return a.div(b);
}

Posit operator-(const Posit& a)
{
    return a.neg();
}

bool operator<(const Posit&a , const Posit& b)
{
    return a.lt(b);
}

bool operator<=(const Posit&a , const Posit& b)
{
    return a.le(b);
}

bool operator>(const Posit&a , const Posit& b)
{
    return a.gt(b);
}

bool operator>=(const Posit&a , const Posit& b)
{
    return a.ge(b);
}
bool operator==(const Posit&a , const Posit& b)
{
    return a.eq(b);
}

bool operator!=(const Posit&a , const Posit& b)
{
    return !a.eq(b);
}

POSIT_UTYPE pack_posit(struct unpacked_t up, int nbits, int es)
{
    POSIT_UTYPE p;
    POSIT_UTYPE regbits;
    POSIT_UTYPE expbits;

    // handle underflow and overflow.
    // in either case, exponent and fraction bits will disappear.
    int maxexp = POW2(es) * (nbits - 2);
    if (up.exp < -maxexp) {
        up.exp = -maxexp;
    } else if (up.exp > maxexp) {
        up.exp = maxexp;
    }

    int reg = FLOORDIV(up.exp, POW2(es));
    int ss = util_ss();
    int rs = MAX(-reg + 1, reg + 2);

    // FIXME: round exponent up if needed
    if (ss + rs + es >= nbits && up.frac >= POSIT_MSB) {
        up.exp++;
        reg = FLOORDIV(up.exp, POW2(es));
        rs = MAX(-reg + 1, reg + 2);
    }

    POSIT_UTYPE exp = up.exp - POW2(es) * reg;

    if (reg < 0) {
        regbits = RSHIFT(POSIT_MSB, -reg);
    } else {
        regbits = LMASK(POSIT_MASK, reg + 1);
    }
    expbits = LMASK(LSHIFT(exp, POSIT_WIDTH - es), es);

    p = up.frac;
    p = expbits | RSHIFT(p, es);
    p = regbits | RSHIFT(p, rs);
    p = RSHIFT(p, ss);

    if (up.neg) {
        return util_neg(p, nbits);
    } else {
        return LMASK(p, nbits);
    }
}

float pack_float(struct unpacked_t up)
{
    int fexp = up.exp + 127;

    // left aligned
    uint32_t fexpbits;
    uint32_t ffracbits;

    if (fexp > 254) {
        // overflow, set maximum value
        fexpbits = LSHIFT(254, 24);
        ffracbits = -1;
    } else if (fexp < 1) {
        // underflow, pack as denormal
        fexpbits = 0;
#if POSIT_WIDTH <= 32
        ffracbits = LSHIFT((uint32_t)(POSIT_MSB | RSHIFT(up.frac, 1)), 32 - POSIT_WIDTH);
#else
        ffracbits = RSHIFT(POSIT_MSB | RSHIFT(up.frac, 1), POSIT_WIDTH - 32);
#endif
        ffracbits = RSHIFT(ffracbits, -fexp);
    } else {
        fexpbits = LSHIFT(fexp & 0xFF, 24);
#if POSIT_WIDTH <= 32
        ffracbits = LSHIFT((uint32_t)up.frac, 32 - POSIT_WIDTH);
#else
        ffracbits = RSHIFT(up.frac, POSIT_WIDTH - 32);
#endif
    }

    union {
        float f;
        uint32_t u;
    } un;

    un.u = ffracbits;
    un.u = fexpbits | RSHIFT(un.u, 8);
    un.u = LSHIFT((uint32_t)up.neg, 31) | RSHIFT(un.u, 1);

    // don't underflow to zero
    if (LSHIFT(un.u, 1) == 0) {
        un.u++;
    }

    return un.f;
}

double pack_double(struct unpacked_t up)
{
    int fexp = up.exp + 1023;

    // left aligned
    uint64_t fexpbits;
    uint64_t ffracbits;

    if (fexp > 2046) {
        // overflow, set maximum value
        fexpbits = LSHIFT((uint64_t)2046, 53);
        ffracbits = -1;
    } else if (fexp < 1) {
        // underflow, pack as denormal
        fexpbits = 0;
#if POSIT_WIDTH <= 64
        ffracbits = LSHIFT((uint64_t)(POSIT_MSB | RSHIFT(up.frac, 1)), 64 - POSIT_WIDTH);
#else
        ffracbits = RSHIFT(POSIT_MSB | RSHIFT(up.frac, 1), POSIT_WIDTH - 64);
#endif
        ffracbits = RSHIFT(ffracbits, -fexp);
    } else {
        fexpbits = LSHIFT((uint64_t)(fexp & 0x7FF), 53);
#if POSIT_WIDTH <= 64
        ffracbits = LSHIFT((uint64_t)up.frac, 64 - POSIT_WIDTH);
#else
        ffracbits = RSHIFT(up.frac, POSIT_WIDTH - 64);
#endif
    }

    union {
        double f;
        uint64_t u;
    } un;

    un.u = ffracbits;
    un.u = fexpbits | RSHIFT(un.u, 11);
    un.u = LSHIFT((uint64_t)up.neg, 63) | RSHIFT(un.u, 1);

    // don't underflow to zero
    if (LSHIFT(un.u, 1) == 0) {
        un.u++;
    }

    return un.f;
}

struct unpacked_t unpack_posit(POSIT_UTYPE p, int nbits, int es)
{
    struct unpacked_t up;

    bool neg = util_is_neg(p);
    if (neg) {
        p = util_neg(p, nbits);
    }

    int ss = util_ss();
    int rs = util_rs(p, nbits);

    int lz = CLZ(LSHIFT(p, ss));
    int lo = CLZ(LSHIFT(~p, ss) | 1); // add LSB to compensate for sign bit

    int reg = (lz == 0 ? lo - 1 : -lz);
    POSIT_UTYPE exp = RSHIFT(LSHIFT(p, ss + rs), POSIT_WIDTH - es);

    up.neg = neg;
    up.exp = POW2(es) * reg + exp;
    up.frac = LSHIFT(p, ss + rs + es);

    return up;
}

struct unpacked_t unpack_float(float f)
{
    struct unpacked_t up;
    int bias = 127;

    union {
        float f;
        uint32_t u;
    } un;

    un.f = f;

    up.neg = RSHIFT(un.u, 31);
    up.exp = (RSHIFT(un.u, 23) & 0xFF) - bias;
#if POSIT_WIDTH <= 32
    up.frac = RSHIFT(LSHIFT(un.u, 9), 32 - POSIT_WIDTH);
#else
    up.frac = LSHIFT((POSIT_UTYPE)un.u, POSIT_WIDTH - 32 + 9);
#endif

    if (up.exp == -bias) {
        // normalize
        // FIXME: some precision is lost if frac was downcasted
        up.exp -= CLZ(up.frac);
        up.frac = LSHIFT(up.frac, CLZ(up.frac) + 1);
    }

    return up;
}

struct unpacked_t unpack_double(double f)
{
    struct unpacked_t up;
    int bias = 1023;

    union {
        double f;
        uint64_t u;
    } un;

    un.f = f;

    up.neg = RSHIFT(un.u, 63);
    up.exp = (RSHIFT(un.u, 52) & 0x7FF) - bias;
#if POSIT_WIDTH <= 64
    up.frac = RSHIFT(LSHIFT(un.u, 12), 64 - POSIT_WIDTH);
#else
    up.frac = LSHIFT((POSIT_UTYPE)un.u, POSIT_WIDTH - 64 + 12);
#endif

    if (up.exp == -bias) {
        // normalize
        // FIXME: some precision is lost if frac was downcasted
        up.exp -= CLZ(up.frac);
        up.frac = LSHIFT(up.frac, CLZ(up.frac) + 1);
    }

    return up;
}

static struct unpacked_t add(struct unpacked_t a, struct unpacked_t b, bool neg)
{
    struct unpacked_t r;

    POSIT_LUTYPE afrac = HIDDEN_BIT(a.frac);
    POSIT_LUTYPE bfrac = HIDDEN_BIT(b.frac);
    POSIT_LUTYPE frac;

    if (a.exp > b.exp) {
        r.exp = a.exp;
        bfrac = RSHIFT(bfrac, a.exp - b.exp);
    } else {
        r.exp = b.exp;
        afrac = RSHIFT(afrac, b.exp - a.exp);
    }

    frac = afrac + bfrac;
    if (RSHIFT(frac, POSIT_WIDTH) != 0) {
        r.exp++;
        frac = RSHIFT(frac, 1);
    }

    r.neg = neg;
    r.frac = LSHIFT(frac, 1);

    return r;
}

static struct unpacked_t sub(struct unpacked_t a, struct unpacked_t b, bool neg)
{
    struct unpacked_t r;

    POSIT_UTYPE afrac = HIDDEN_BIT(a.frac);
    POSIT_UTYPE bfrac = HIDDEN_BIT(b.frac);
    POSIT_UTYPE frac;

    if (a.exp > b.exp || (a.exp == b.exp && a.frac > b.frac)) {
        r.exp = a.exp;
        bfrac = RSHIFT(bfrac, a.exp - b.exp);
        frac = afrac - bfrac;
    } else {
        neg = !neg;
        r.exp = b.exp;
        afrac = RSHIFT(afrac, b.exp - a.exp);
        frac = bfrac - afrac;
    }

    r.neg = neg;
    r.exp -= CLZ(frac);
    r.frac = LSHIFT(frac, CLZ(frac) + 1);

    return r;
}

struct unpacked_t op2_mul(struct unpacked_t a, struct unpacked_t b)
{
    struct unpacked_t r;

    POSIT_LUTYPE afrac = HIDDEN_BIT(a.frac);
    POSIT_LUTYPE bfrac = HIDDEN_BIT(b.frac);
    POSIT_UTYPE frac = RSHIFT(afrac * bfrac, POSIT_WIDTH);
    POSIT_STYPE exp = a.exp + b.exp + 1;

    if ((frac & POSIT_MSB) == 0) {
        exp--;
        frac = LSHIFT(frac, 1);
    }

    r.neg = a.neg ^ b.neg;
    r.exp = exp;
    r.frac = LSHIFT(frac, 1);

    return r;
}

struct unpacked_t op2_div(struct unpacked_t a, struct unpacked_t b)
{
    struct unpacked_t r;

    POSIT_LUTYPE afrac = HIDDEN_BIT(a.frac);
    POSIT_LUTYPE bfrac = HIDDEN_BIT(b.frac);
    POSIT_STYPE exp = a.exp - b.exp;

    if (afrac < bfrac) {
        exp--;
        bfrac = RSHIFT(bfrac, 1);
    }

    r.neg = a.neg ^ b.neg;
    r.exp = exp;
    r.frac = LSHIFT(afrac, POSIT_WIDTH) / bfrac;

    return r;
}

struct unpacked_t op2_add(struct unpacked_t a, struct unpacked_t b)
{
    if (a.neg == b.neg) {
        return add(a, b, a.neg);
    } else {
        return sub(a, b, a.neg);
    }
}

struct unpacked_t op2_sub(struct unpacked_t a, struct unpacked_t b)
{
    if (a.neg == b.neg) {
        return sub(a, b, a.neg);
    } else {
        return add(a, b, a.neg);
    }
}

static struct unpacked_t half(struct unpacked_t a)
{
    struct unpacked_t r = a;

    r.exp--;

    return r;
}

struct unpacked_t op1_sqrt(struct unpacked_t a)
{
    struct unpacked_t r = a;

    // initial guess: half exponent is the sqrt if we ignore fraction bits
    r.exp /= 2;

    for (int i = 0; i < 100; i++) {
        struct unpacked_t rn;

        // Newton-Raphson: rn = r - (r^2 - a) / (2 * r) = (r + a / r) / 2
        rn = half(op2_add(r, op2_div(a, r)));

        if (rn.exp == r.exp && rn.frac == r.frac) {
            break;
        }

        r = rn;
    }

    return r;
}

bool util_is_zero(POSIT_UTYPE p)
{
    return p == POSIT_ZERO;
}

bool util_is_nar(POSIT_UTYPE p)
{
    return p == POSIT_NAR;
}

bool util_is_neg(POSIT_UTYPE p)
{
    return (POSIT_STYPE)p < 0 && !util_is_nar(p);
}

int util_ss()
{
    return 1;
}

int util_rs(POSIT_UTYPE p, int nbits)
{
    int ss = util_ss();
    int lz = CLZ(LSHIFT(p, ss));
    int lo = CLZ(LSHIFT(~p, ss));
    int rs = MAX(lz, lo) + 1;

    return MIN(rs, nbits - ss);
}

int util_es(POSIT_UTYPE p, int nbits, int es)
{
    int ss = util_ss();
    int rs = util_rs(p, nbits);

    return MIN(MAX(nbits - ss - rs, 0), es);
}

int util_fs(POSIT_UTYPE p, int nbits, int es)
{
    int ss = util_ss();
    int rs = util_rs(p, nbits);

    return MAX(nbits - ss - rs - es, 0);
}

POSIT_UTYPE util_neg(POSIT_UTYPE p, int nbits)
{
    // reverse all bits and add one
    return LMASK(-LMASK(p, nbits), nbits);
}

int main(int argc, char *argv[])
{
    auto a = Posit32(0.1);
    auto b = Posit32(0.2);
    auto c = a+b;
    c.print(); //build with command emcc test.cpp -s WASM=1 -o test.html -std=c++11
    return 0;
}
