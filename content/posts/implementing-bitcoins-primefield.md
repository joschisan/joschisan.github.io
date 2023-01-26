---
title: "Implementing secp256k1's Primefield"
date: 2023-01-04T19:01:41+01:00
draft: false
---

In the following we assume familiarity with the mathematical concept of a Primefield. 
Bitcoins choosen prime is
$$p = 2^{256} - 2^{32} - 977$$
hence 256 bits would be sufficient to represent any field element. Therefore we could 
encode a field element in only four u64 unsigned integers, called limbs. This, however,
would force us too handle possible overflow on any additon of two limbs as we would 
have to use their full range, which would slow down the addition of two field elements
substantially. Furthermore, the product of two limbs is at most a 128 bit integer. 
Hence it can be represented by Rust's u128 integer type, however, the sum of two of such
products might be a 129 bit integer, again requiring us to handle possible overflow. 
This would slow down the multiplication of two field elements which requires us to multiply
their limbs and summ over the products. To prevent this we therefore represent a field
element $x \in F_p$ with five u64 limbs $l_0,l_1,l_2,l_3,l_4$. We call these limbs
the representation of $x$ and define the representations value as
$$\sum_{k \le 4} l_k \cdot 2^{52k}$$
For the representation of 
$$p_0 = FFFFEFFFFFC2F_{16}$$
$$p_1 = FFFFFFFFFFFFF_{16}$$
$$p_2 = FFFFFFFFFFFFF_{16}$$
$$p_3 = FFFFFFFFFFFFF_{16}$$
$$p_4 = 0FFFFFFFFFFFF_{16}$$
in base 16 we get the value of our prime $p$ itself. Notice that we can now have multiple
possible representations with the same value as the highest 12 bits of $l_k$ overlap with
the lowest 12 bits of $l_{k + 1}$.  The value of a representation may be any 272 bit 
unsigned int. Hence, the map from representations to field elements defined by 
$$(l_k) \mapsto \sum_{k \le 4} l_k \cdot 2^{52k} \in F_p$$
is a surjection. Using five limbs we can now represent up to 272 bit
unsigned integers using only the lowest 52 bits of limbs 1 to 3 and 48 bits of limbs 4.
Therefore, we can compute the sum of any two elements without any overflow if 
we choose representations with limbs of at least one leading zero. Furthermore the product 
of two 53 bit limbs would be at most a 104 bit integer, leaving us plenty of headroom too add 
several of those products when using the u128 integer type. To assure that we maintain sufficient
headroom to prevent overflow we track the the maximal number of bits of a representations limbs, 
which we call the representations magnitude. Finally, we define a field element as

```rust
pub struct Element {
    pub limbs: [u64; 5],
    pub magnitude: u64
}

```
and implement the method

```rust
pub fn verify(self) {
        assert!((self.limbs[0] >> (self.magnitude + 52)) == 0);
        assert!((self.limbs[1] >> (self.magnitude + 52)) == 0);
        assert!((self.limbs[2] >> (self.magnitude + 52)) == 0);
        assert!((self.limbs[3] >> (self.magnitude + 52)) == 0);
        assert!((self.limbs[4] >> (self.magnitude + 49)) == 0);
    }
```
to verify that we track our headroom correctly. We can now implement 
additon without any overflow by overloading the + operator:
```rust
 fn add(self, rhs: Self) -> Self {
        assert!(self.magnitude < 12);
        assert!(rhs.magnitude < 12);
        self.verify();
        rhs.verify();

        Self {
            limbs: [0, 1, 2, 3, 4].map(|i| self.limbs[i] + rhs.limbs[i]),
            magnitude: cmp::max(self.magnitude, rhs.magnitude) + 1,
        }
    }

```
Furthermore we implemement a method to calculate the additive inverse of an element:
```rust
 pub fn negative(self) -> Self {
        const P: [u64; 5] = [
            0xFFFFEFFFFFC2F,
            0xFFFFFFFFFFFFF,
            0xFFFFFFFFFFFFF,
            0xFFFFFFFFFFFFF,
            0xFFFFFFFFFFFF,
        ];

        assert!(self.magnitude < 12);
        self.verify();

        Self {
            limbs: [0, 1, 2, 3, 4].map(|i| (P[i] << (self.magnitude + 1)) - self.limbs[i]),
            magnitude: self.magnitude + 1,
        }
    }
```

When we run out of headroom after several operations such as additions or negation we 
then subtract a multiple of p from our representation to reduce the maximal number of 
bits in each limb. We recall that $p = 2^{256} - r$ for $r = 2^{32} + 977$. Hence, 
when we subtract a $m2^{256}$ for $m = l_k >> 48$ from our representation by setting
the 16 most significant bits of $l_4$ to zero, we need to add $mr$ to our representation
to maintain congruency modulo p. This is where the choice of prime becomes important for 
the effifiency of this implementation as we can fit the product $mr$ of at most 45bits
in a u64, which can now be added to $l_0$. Therefore the Prime of secp256k1 has been 
specifically choosen such that $2^{256} - p$ has a low number of bits. To prevent overflow
we require at least a single bit of headroom in all limbs when calling the reduce method:

```rust
pub fn reduce(self) -> Self {
        const MASK_52: u64 = 0xFFFFFFFFFFFFF;
        const MASK_48: u64 = 0xFFFFFFFFFFFF;
        const R: u64 = 0x1000003D1;

        assert!(self.magnitude < 12);
        self.verify();

        let mut limbs = self.limbs;

        limbs[0] += (limbs[4] >> 48) * R;
        limbs[4] &= MASK_48;

        limbs[1] += limbs[0] >> 52;
        limbs[2] += limbs[1] >> 52;
        limbs[3] += limbs[2] >> 52;
        limbs[4] += limbs[3] >> 52;

        limbs[0] &= MASK_52;
        limbs[1] &= MASK_52;
        limbs[2] &= MASK_52;
        limbs[3] &= MASK_52;

        Self {
            limbs,
            magnitude: 0,
        }
    }

```
By masking $l_4$ we obtain $l_4 \le 2^{48}$. After we have added a 12 bit carry from 
$l_3$ we can derive 
$$l_4 \le 2^{48} + 2^{12} \le 2^{49} - 2 = 2(2^{48} - 1) = 2p_4$$
where $p_4$ is the most significant limbs of our primes representation as defined earlier.
Therefore the value of a reduced representation is smaller the 2p. However, this bound
is not sufficient to compute a unique representation for any field element. To achieve
this we need to find a äquivalent representation of value is smaller then p, which we 
call normalized. Hence a representation represents the zero element if and only if its
corresponding normalized representations value is zero. We will now implement
a method to compute the normal representation of a given representation.
```rust
pub fn normalize(self) -> Self {
        const MASK_52: u64 = 0xFFFFFFFFFFFFF;
        const MASK_48: u64 = 0xFFFFFFFFFFFF;
        const R: u64 = 0x1000003D1;
        const P: [u64;5] = [
            0xFFFFEFFFFFC2F,
            0xFFFFFFFFFFFFF,
            0xFFFFFFFFFFFFF,
            0xFFFFFFFFFFFFF,
            0xFFFFFFFFFFFF
        ];

        // after reduction the value of of a representation is less than 2p
        // since l_4 <= MASK_48 + 2^12 < 2 * MASK_48 and the bits of l_3
        // that overlap with l_4 are zero
        let mut limbs = self.reduce().limbs;


        // Since there is no overlap in the limbs we can compare the arrays in
        // colexicographical order to determine if the representation is normalized
        let normalized = limbs.iter().rev().lt(P.iter().rev());

        if !normalized {
            // Since a reduced representations value is less then 2p we
            // have to subtract P only once to achieve a normalized representation.
            limbs[0] += R;
            limbs[1] += limbs[0] >> 52;
            limbs[2] += limbs[1] >> 52;
            limbs[3] += limbs[2] >> 52;
            limbs[4] += limbs[3] >> 52;

            limbs[0] &= MASK_52;
            limbs[1] &= MASK_52;
            limbs[2] &= MASK_52;
            limbs[3] &= MASK_52;

            // Since the value of a representaion is larger or equal to P
            // and P + R = 2^256 we carry to bit 49 of limb 4 if it was not
            // set to 1 already. Furthermore the representations value is 
            // smaller then 2p + R = 2^256 + P <= 2^257
            assert!(limbs[4] >> 48 == 1);

            // by zeroing bit 49 we subtract 2^256 and maintain congruecy modulo P
            limbs[4] &= MASK_48;
        };

        // verify the representation is normalized
        assert!(limbs.iter().rev().lt(P.iter().rev()));

        Self {
            limbs,
            magnitude: 0,
        }
 ```
Now two elements are equal if and only if the normal representation of their difference 
is zero and we can implement the PartialEq trait for `Element`.
```rust
impl PartialEq for Element {
    fn eq(&self, rhs: &Self) -> bool {
        (*self - *rhs).normalize().limbs == [0,0,0,0,0]
    }
}
```
Let $a,b$ be two field elements
represented by the limbs $a_i,b_i$. In the following we define the notation
$$[c_n ... c_0] = \sum_{k \leq n} c_k2^{52} \in F_p$$
to track the computation in the implementation. Notice that 
$$[c_5 ... c_0] = [c_4 ... c_0 + c_5r]$$
for $p = 2^{256} - r$. Once again, the efficiency of the implementation
depends on $r$ being a small number and we can use the headroom of zeros
in our limbs to prevent overflow. Now let
$$p_k = \sum_{i + j = k} a_i b_i$$
then we can can compute the product of $a$ and $b$ by
$$ [p_7 ... p_0] \in F_p$$
Finally, we implement the multiplication. 

```rust
pub fn multiply(a: field::Element, b: field::Element) -> field::Element {
    const MASK_48: u128 = 0xFFFFFFFFFFFF;
    const MASK_52: u128 = 0xFFFFFFFFFFFFF;
    const MASK_64: u128 = 0xFFFFFFFFFFFFFFFF;
    const R: u128 = 0x1000003D10;

    assert!(a.magnitude < 5);
    assert!(b.magnitude < 5);
    a.verify();
    b.verify();

    let a0 = a.limbs[0] as u128;
    let a1 = a.limbs[1] as u128;
    let a2 = a.limbs[2] as u128;
    let a3 = a.limbs[3] as u128;
    let a4 = a.limbs[4] as u128;

    let b0 = b.limbs[0] as u128;
    let b1 = b.limbs[1] as u128;
    let b2 = b.limbs[2] as u128;
    let b3 = b.limbs[3] as u128;
    let b4 = b.limbs[4] as u128;

    let mut d = a0 * b3 + a1 * b2 + a2 * b1 + a3 * b0;
    assert!(d >> 114 == 0);
    // [d 0 0 0] = [p3 0 0 0]

    let mut c = a4 * b4;
    assert!(c >> 112 == 0);
    // [c 0 0 0 0 d 0 0 0] = [p8 0 0 0 0 p3 0 0 0]

    d += (c & MASK_64) * R;
    c >>= 64;
    assert!(d >> 115 == 0);
    assert!(c >> 48 == 0);
    // [(c<<12) 0 0 0 0 0 d 0 0 0] = [p8 0 0 0 0 p3 0 0 0]

    let t3 = d & MASK_52;
    d >>= 52;
    assert!(t3 >> 52 == 0);
    assert!(d >> 63 == 0);
    // [(c<<12) 0 0 0 0 d t3 0 0 0] = [p8 0 0 0 0 p3 0 0 0]

    d += a0 * b4 + a1 * b3 + a2 * b2 + a3 * b1 + a4 * b0;
    assert!(d >> 115 == 0);
    // [(c<<12) 0 0 0 0 d t3 0 0 0] = [p8 0 0 0 p4 p3 0 0 0]

    d += (c & MASK_64) * (R << 12);
    assert!(d >> 116 == 0);
    // [d t3 0 0 0] = [p8 0 0 0 p4 p3 0 0 0]

    let mut t4 = d & MASK_52;
    d >>= 52;
    assert!(t4 >> 52 == 0);
    assert!(d >> 64 == 0);
    // [d t4 t3 0 0 0] = [p8 0 0 0 p4 p3 0 0 0]

    let tx = t4 >> 48;
    t4 &= MASK_48;
    assert!(tx >> 4 == 0);
    assert!(t4 >> 48 == 0);
    // [d (t4 + (tx << 48)) t3 0 0 0] = [p8 0 0 0 p4 p3 0 0 0] */
    //
    c = a0 * b0;
    assert!(c >> 112 == 0);
    // [d t4+(tx<<48) t3 0 0 c] = [p8 0 0 0 p4 p3 0 0 p0]

    d += a1 * b4 + a2 * b3 + a3 * b2 + a4 * b1;
    assert!(d >> 115 == 0);
    // [d (t4 + (tx << 48)) t3 0 0 c] = [p8 0 0 p5 p4 p3 0 0 p0]

    let mut u0 = d & MASK_52;
    d >>= 52;
    assert!(u0 >> 52 == 0);
    assert!(d >> 63 == 0);
    // [d u0 (t4 + (tx << 48)) t3 0 0 c] = [p8 0 0 p5 p4 p3 0 0 p0]
    // [d 0 (t4 + (tx << 48) + (u0 << 52)) t3 0 0 c] = [p8 0 0 p5 p4 p3 0 0 p0]

    u0 = (u0 << 4) | tx;
    assert!(u0 >> 56 == 0);
    // [d 0 t4+(u0<<48) t3 0 0 c] = [p8 0 0 p5 p4 p3 0 0 p0]

    c += u0 * (R >> 4);
    assert!(c >> 115 == 0);
    // [d 0 t4 t3 0 0 c] = [p8 0 0 p5 p4 p3 0 0 p0]

    let r0 = c & MASK_52;
    c >>= 52;
    assert!(r0 >> 52 == 0);
    assert!(c >> 61 == 0);
    // [d 0 t4 t3 0 c r0] = [p8 0 0 p5 p4 p3 0 0 p0]

    c += a0 * b1 + a1 * b0;
    assert!(c >> 114 == 0);
    // [d 0 t4 t3 0 c r0] = [p8 0 0 p5 p4 p3 0 p1 p0]

    d += a2 * b4 + a3 * b3 + a4 * b2;
    assert!(d >> 114 == 0);
    // [d 0 t4 t3 0 c r0] = [p8 0 p6 p5 p4 p3 0 p1 p0]

    c += (d & MASK_52) * R;
    d >>= 52;
    assert!(c >> 115 == 0);
    assert!(d >> 62 == 0);
    // [d 0 0 t4 t3 0 c r0] = [p8 0 p6 p5 p4 p3 0 p1 p0]

    let r1 = c & MASK_52;
    c >>= 52;
    assert!(r1 >> 52 == 0);
    assert!(c >> 63 == 0);
    // [d 0 0 t4 t3 c r1 r0] = [p8 0 p6 p5 p4 p3 0 p1 p0]

    c += a0 * b2 + a1 * b1 + a2 * b0;
    assert!(c >> 114 == 0);
    // [d 0 0 t4 t3 c r1 r0] = [p8 0 p6 p5 p4 p3 p2 p1 p0]

    d += a3 * b4 + a4 * b3;
    assert!(d >> 114 == 0);
    // [d 0 0 t4 t3 c t1 r0] = [p8 p7 p6 p5 p4 p3 p2 p1 p0]

    c += (d & MASK_64) * R;
    d >>= 64;
    assert!(c >> 115 == 0);
    assert!(d >> 50 == 0);
    // [(d<<12) 0 0 0 t4 t3 c r1 r0] = [p8 p7 p6 p5 p4 p3 p2 p1 p0]

    let r2 = c & MASK_52;
    c >>= 52;
    assert!(r2 >> 52 == 0);
    assert!(c >> 63 == 0);
    // [(d<<12) 0 0 0 t4 (t3 + c) r2 r1 r0] = [p8 p7 p6 p5 p4 p3 p2 p1 p0]

    c += (d & MASK_64) * (R << 12);
    c += t3;
    assert!(c >> 100 == 0);
    // [t4 c r2 r1 r0] = [p8 p7 p6 p5 p4 p3 p2 p1 p0]

    let r3 = c & MASK_52;
    c >>= 52;
    assert!(r3 >> 52 == 0);
    assert!(c >> 48 == 0);
    // [t4+c r3 r2 r1 r0] = [p8 p7 p6 p5 p4 p3 p2 p1 p0]

    let r4 = (c & MASK_64) + t4;
    assert!(r4 >> 49 == 0);
    // [r4 r3 r2 r1 r0] = [p8 p7 p6 p5 p4 p3 p2 p1 p0]

    field::Element {
        limbs: [r0 as u64, r1 as u64, r2 as u64, r3 as u64, r4 as u64],
        magnitude: 0,
    }
}

```
