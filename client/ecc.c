#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct{

	mpz_t x;
	mpz_t y;

} point_t;

typedef struct{

	point_t *generator;
	mpz_t a;
	mpz_t b;
	mpz_t p;
	mpz_t cofactor;
} curve_t;

/*
 ----------------- Modular arifmetic -----------------
*/


/*
*   Function solves a first degree modulo comparison
*   over a prime field Fp:
*
*   +++++++++++++++++++++++++++++++++++++++++++++++++
*   +                                               +
*   +    1)  a*x = b (mod p)                        +
*   +       (a, p) = 1, because p - primary         +
*   +                                               +
*   +    2)  a*u + p*v = 1 - Co-prime criterion,    +
*   +        u and v - Bezout odds.                 +
*   +                                               +
*   +    3)  a*u = 1 (mod p) | * b                  +
*   +                                               +
*   +    4)  a*(u*b) = b (mod p)                    +
*   +                                               +
*   +    5)  x = (u*b) (mod p)                      +
*   +                                               +
*   +++++++++++++++++++++++++++++++++++++++++++++++++
*
*   Function time: Extended Euclidean algorithm.
*
*   Function writing result to lymda.
*/

void
solution(mpz_t a, mpz_t b, mpz_t p, mpz_t lymda){

	/*extended Euclidean elgoritm*/
    mpz_t s, old_s, t, old_t, quotient, temp, r, old_r;

    mpz_init(s);
    mpz_init(old_s);
    mpz_init(t);
    mpz_init(old_t);
    mpz_init(r);
    mpz_init(old_r);
    mpz_init(quotient);
    mpz_init(temp);

    mpz_set_si(s, 0);
    mpz_set_si(old_s, 1);
    mpz_set_si(t, 1);
    mpz_set_si(old_t, 0);
    mpz_set(r, p);
    mpz_set(old_r, a);
    mpz_set_si(quotient, 0);
    mpz_set_si(temp, 0);

repeat:

    mpz_fdiv_q(quotient,old_r,r);

    mpz_set(temp,old_r);

    mpz_set(old_r,r);

    mpz_mul(r,quotient,r);
    mpz_sub(r,temp,r);

    mpz_set(temp,old_s);

    mpz_set(old_s,s);

    mpz_mul(s,quotient,s);
    mpz_sub(s,temp,s);

    mpz_set(temp,old_t);

    mpz_set(old_t,t);

    mpz_mul(t,quotient,t);
    mpz_sub(t,temp,t);

	if(mpz_cmp_si(r, 0))
        goto repeat;

	mpz_mul(temp,old_s,b);
    mpz_mod(lymda,temp,p);

    mpz_clear(s);
    mpz_clear(old_s);
    mpz_clear(t);
    mpz_clear(old_t);
    mpz_clear(r);
    mpz_clear(old_r);
    mpz_clear(quotient);
    mpz_clear(temp);
}

/*
 ----------------- Functions that work with points on curve -----------------
*/

/*
*   This functions check necessary condition:
*   y^2 = x^3 + ax + b (entry in Wieierstrass form)
*
*   Return value if successful is integer <> 0.
*/

int
is_on_curve(point_t *point, mpz_t a, mpz_t b, mpz_t p){

	if(point == NULL)
		return 1;

    mpz_t x, y;

    mpz_init(x);
    mpz_init(y);

    mpz_set(x, point->x);
    mpz_set(y, point->y);

    mpz_mul(y,y,y);
    mpz_mul(x,x,x);
    mpz_mul(x,x,x);
    mpz_sub(y,y,x);
    mpz_sub(y,y,a);
    mpz_sub(y,y,b);
    mpz_fdiv_q(y,y,p);

    return mpz_cmp_si(y, 0);
}

/*
*   Where the points P and Q are coincident (at the same coordinates), addition is similar,
*   except that there is no well-defined straight line through P,
*   so the operation is closed using a limiting case, the tangent to the curve, E, at P.
*
*   This is calculated as above, taking derivatives (dE/dx)/(dE/dy):
*   Lymda = 3x1^2 + a / 2y1,
*   where a is from the defining equation of the curve, E, above.
*
*   x = lymda^2 - 2x1
*   y = lymda(x1 - x) - y1
*   Function return struct point_t, that countains calculated point(x,y).
*/

point_t
*point_doubling(point_t *point, mpz_t a, mpz_t p){

    mpz_t x, y, lymda, old_x, old_y, old_lymda;

    mpz_init(x);
    mpz_init(y);
    mpz_init(lymda);
    mpz_init(old_x);
    mpz_init(old_y);
    mpz_init(old_lymda);

    mpz_set(old_x, point->x);
    mpz_set(old_y, point->y);

    mpz_mul(x,old_x,old_x);
    mpz_mul_si(x,x,3);
    mpz_add(x,x,a);
    mpz_mul_si(y,old_y,2);
    solution(y,x,p,old_lymda);

    mpz_mul(lymda,old_lymda,old_lymda);
    mpz_mul_si(x,old_x,2);
    mpz_sub(x,lymda,x);
    mpz_mod(point->x,x,p);

    mpz_sub(x,old_x,x);
    mpz_mul(lymda,old_lymda,x);
    mpz_sub(y,lymda,old_y);
    mpz_mod(point->y,y,p);

    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(lymda);
    mpz_clear(old_x);
    mpz_clear(old_y);
    mpz_clear(old_lymda);

    return point;
}


/*
*   With 2 distinct points, P and Q,
*   addition is defined as the negation
*   of the point resulting from the intersection of the curve, E,
*   and the straight line defined by the points P and Q, giving the point, R.
*
*   P + Q = R | point1 + point2 = return point
*   (Xp , Yp) + (Xq , Yq) = (Xr , Yr) | (x1 , y1) + (x2 , y2) = (x3 , y3)
*
*   Assuming the elliptic curve, E,
*   is given by y2 = x3 + ax + b, this can be calculated as:
*   Lymda = y1-y2 / x1-x2
*   Xr = Lymda^2 - x1 - x2
*   Yr = lymda(x1 - x3) - y1
*
*   These equations are correct when neither point is the point at infinity, O,
*   and if the points have different x coordinates (they're not mutual inverses).
*   This is important for the ECDSA verification algorithm where the hash value could be zero.
*
*   Function return struct point_t, that countains calculated point(x,y).
*/

point_t
*point_add(point_t *point1, point_t *point2, mpz_t p){

    mpz_t x, y, x1, x2, y1, y2, lymda, old_lymda;

    mpz_init(x1);
    mpz_init(x2);
    mpz_init(y1);
    mpz_init(y2);
    mpz_init(lymda);
    mpz_init(old_lymda);
    mpz_init(x);
    mpz_init(y);

    mpz_set(x1, point1->x);
    mpz_set(x2, point2->x);
    mpz_set(y1, point1->y);
    mpz_set(y2, point2->y);

    mpz_sub(x,x2,x1);
    mpz_sub(y,y2,y1);
    solution(x,y,p,old_lymda);

    mpz_mul(lymda,old_lymda,old_lymda);
    mpz_sub(x,lymda,x1);
    mpz_sub(x,x,x2);
    mpz_mod(point1->x,x,p);

    mpz_sub(x,x1,point1->x);
    mpz_mul(lymda,old_lymda,x);
    mpz_sub(y,lymda,y1);
    mpz_mod(point1->y,y,p);

    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(x1);
    mpz_clear(x2);
    mpz_clear(y1);
    mpz_clear(y2);
    mpz_clear(lymda);
    mpz_clear(old_lymda);
    return point1;

}


/*
*   Function check the points and selects the method:
*
*       1) point_doubling
*
*       2) point_add
*
*   Function return struct point that countains x and y coordinates.
*/

point_t
*point_addition(point_t *point1,point_t *point2,
                mpz_t a, mpz_t b, mpz_t p){

//	if(!is_on_curve(point1, a, b, p)){}
//	if(!is_on_curve(point2, a, b, p)){}

	if(point1 == NULL)
		return point2;

	if(point2 == NULL)
		return point1;

	if(mpz_cmp(point1->x, point2->x) == 0)
        if(mpz_cmp(point1->y, point2->y) != 0)
            return NULL;

    if(mpz_cmp(point1->x, point2->x) != 0)
        return point_add(point1, point2, p);


	return point_doubling(point1, a, p);
}

// return -point
point_t
*point_neg(point_t *point, mpz_t p){

    mpz_mul_si(point->y, point->y, -1);
    mpz_fdiv_q(point->y, point->y, p);

    return point;
}


/*
 -------------------- Calculate ECC point --------------------
*/

/*
*   Function calculate a public key on EC with prime field Fp:
*
*   +++++++++++++++++++++++++
*   +                       +
*   +    P = G^K (mod p)    +
*   +                       +
*   +++++++++++++++++++++++++
*
*   Function time: O(log_2(k))
*   Returns k*point computed using the point_doubling and
*   point_add algorithm.
*/

point_t
*scalar_multiply(point_t *point, mpz_t k, mpz_t a, mpz_t b, mpz_t p){

//	if(!is_on_curve(point, a, b, p))
//		return NULL;

	if(point == NULL)
		return NULL;

	if(mpz_cmp_si(k, 0) < 0){
		mpz_mul_si(k, k, -1);
        return scalar_multiply(point_neg(point, p), k, a, b, p);
    }

	point_t *addend = point, *result = NULL;

repeat:

    if(mpz_odd_p(k))
        if(result != NULL){
            result = point_addition(result, addend, a, b, p);
        } else {
            result = (point_t*)malloc(sizeof(point_t));
            mpz_init(result->x);
            mpz_init(result->y);
            mpz_set(result->x, addend->x);
            mpz_set(result->y, addend->y);
        }

    addend = point_doubling(addend, a, p);
    mpz_divexact_ui(k, k, 2);

    if(mpz_cmp_si(k, 0))
        goto repeat;

//	if(!is_on_curve(result, a, b, p))
//		return NULL;

	return result;
}

/*int
main(int argc, char *argv[], char **env){


    curve_t curve;
    curve.generator = (point_t*)malloc(sizeof(point_t));
    mpz_t secure_key;
    mpz_init(curve.generator->x); mpz_init(curve.generator->y);
    mpz_init(curve.a); mpz_init(curve.b); mpz_init(curve.p); mpz_init(curve.cofactor);
    mpz_init(secure_key);

    const char *generator_x = "602046282375688656758213480587526111916698976636884684818";
    const char *generator_y = "174050332293622031404857552280219410364023488927386650641";
    const char *a = "6277101735386680763835789423207666416083908700390324961276";
    const char *b = "2455155546008943817740293915197451784769108058161191238065";
    const char *p = "6277101735386680763835789423207666416083908700390324961279";
    mpz_set_str(curve.generator->x, generator_x, 10);
    mpz_set_str(curve.generator->y, generator_y, 10);
    mpz_set_str(curve.a, a, 10);
    mpz_set_str(curve.b, b, 10), mpz_set_str(curve.p, p, 10);


    const char *key = argv[1];
    mpz_set_str(secure_key, key, 10);


    point_t *public_key = scalar_multiply(curve.generator, secure_key, curve.a, curve.b, curve.p);


    gmp_printf("%Zd;%Zd", public_key->x, public_key->y);

    return 0;
}*/
