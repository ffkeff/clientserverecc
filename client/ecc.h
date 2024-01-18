#ifndef ECC_H
#define ECC_H

#include <gmp.h>

typedef struct{

	mpz_t x;
	mpz_t y;
} point_t;

void solution(mpz_t a, mpz_t b, mpz_t p, mpz_t lymda);

int is_on_curve(point_t *point, mpz_t a, mpz_t b, mpz_t p);

point_t *point_doubling(point_t *point, mpz_t a, mpz_t p);

point_t *point_add(point_t *point1, point_t *point2, mpz_t p);

point_t *point_addition(point_t *point1,point_t *point2, mpz_t a, mpz_t b, mpz_t p);

point_t *point_neg(point_t *point, mpz_t p);

point_t *scalar_multiply(point_t *point, mpz_t k, mpz_t a, mpz_t b, mpz_t p);

#endif
