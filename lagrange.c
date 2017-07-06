/* Compile and run: gcc -W -Wall -o lagrange lagrange.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef enum {false, true} bool;

struct point
{
	double x;
	double y;
};

struct points
{
	uint32_t n;
	struct point *p;
};

struct polynomial
{
	uint32_t n;
	double *pol;
};

void ps_free(struct points *ps)
{
	free(ps->p);
	free(ps);
}

void pols_free(struct polynomial *pols, uint32_t n)
{
	uint32_t i;

	for (i = 0; i < n; ++i) {
		free(pols[i].pol);
	}
	free(pols);
}

struct points *parse_file(const char *file)
{
	FILE *f = fopen(file, "r");
	uint32_t n;

	fscanf(f, "%u", &n);
	if (n == 0) {
		return NULL;
	}
	struct points *ret = malloc(sizeof(struct points));
	ret->n = n;
	ret->p = malloc(ret->n * sizeof(struct point));

	size_t i;
	for (i = 0; i < n; ++i) {
		fscanf(f, "%lf", &ret->p[i].x);
		fscanf(f, "%lf", &ret->p[i].y);
	}
	return ret;
}

void ps_to_string(const struct points *ps)
{
	size_t i;

	for (i = 0; i < ps->n; ++i) {
		printf("(%lf, %lf)\n", ps->p[i].x, ps->p[i].y);
	}
	printf("\n");
}

void pol_to_string(const struct polynomial *pol)
{
	int32_t i;

	for (i = pol->n - 1; i >= 0; --i) {
		if (i > 1) {
			printf("%2.12fx^%d + ", pol->pol[i], i);
		} else if (i == 1) {
			printf("%2.12fx + ", pol->pol[i]);
		} else {
			printf("%2.12f", pol->pol[i]);
		}
	}
	printf("\n");
}

void pol_clear(const struct polynomial *pol)
{
	uint32_t i;

	for (i = 0; i < pol->n; ++i) {
		pol->pol[i] = 0.0;
	}
}

uint32_t pol_mult_size(const uint32_t n1, const uint32_t n2)
{
	return n1 + n2 - 1;
}

void pol_mult(const struct polynomial *a, const struct polynomial *b, struct polynomial *ret)
{
	uint32_t i, j;

	if (!b) {
		for (i = 0; i < a->n; ++i) {
			ret->pol[i] = a->pol[i];
		}
		return;
	}

	struct polynomial a_cpy;
	a_cpy.n = a->n;
	a_cpy.pol = malloc(sizeof(double) * a->n);

	for (size_t i = 0; i < a->n; ++i) {
		a_cpy.pol[i] = a->pol[i];
		ret->pol[i] = 0.0;
	}

	for (i = 0; i < a->n; ++i) {
		for (j = 0; j < b->n; ++j) {
			if (i + j < ret->n) {
				ret->pol[i + j] += a_cpy.pol[i] * b->pol[j];
			}
		}
	}

	free(a_cpy.pol);
}

void pol_const_mult(const struct polynomial *ret, const double c)
{
	uint32_t i;

	for (i = 0; i < ret->n; ++i) {
		ret->pol[i] *= c;
	}
}

void pol_sum(const struct polynomial *a, const struct polynomial *b, struct polynomial *ret)
{
	uint32_t i;

	for (i = 0; i < ret->n; ++i) {
		ret->pol[i] = a->pol[i] + b->pol[i];
	}
}

void lagrange(const struct points *ps)
{
	struct polynomial *pols = malloc(sizeof(struct polynomial) * ps->n);
	double *consts = malloc(sizeof(double) * ps->n);
	uint32_t i, k;

	/* Iterate through polynomials */
	for (i = 0; i < ps->n; ++i) {
		pols[i].pol = malloc(sizeof(double) * ps->n);
		pols[i].n = ps->n;
		pol_clear(&pols[i]);
		struct polynomial *tmp_pols = malloc(sizeof(struct polynomial) * (ps->n - 1));
		int32_t ii = (i + 1) % ps->n;
		/* Iterate through sub polynomials */
		for (k = 0; k < ps->n - 1; ++k) {
			tmp_pols[k].pol = malloc(sizeof(double) * 2);
			tmp_pols[k].n = 2;
			tmp_pols[k].pol[0] = -ps->p[ii].x;
			tmp_pols[k].pol[1] = 1.0;
			//pol_to_string(&pols[i]);
			if (k == 0) {
				pol_mult(&tmp_pols[k], NULL, &pols[i]);
				consts[i] = (ps->p[i].x - ps->p[ii].x);
			} else {
				pol_mult(&pols[i], &tmp_pols[k], &pols[i]);
				consts[i] *= (ps->p[i].x - ps->p[ii].x);
			}
			ii = (ii + 1) % ps->n;
		}
		printf("L%d = ", i);
		pol_to_string(&pols[i]);
		printf("c%d = %2.10f\n\n", i, consts[i]);
		pol_const_mult(&pols[i], ps->p[i].y / consts[i]);
		if (i > 0) {
			pol_sum(&pols[0], &pols[i], &pols[0]);
		}
		pols_free(tmp_pols, ps->n - 1);
	}

	printf("P(x) = ");
	pol_to_string(&pols[0]);
	free(consts);
	pols_free(pols, ps->n);
}

int main(int argc, char *argv[])
{
	if (argc != 2) {
		printf("Wrong number of parameters!\n");
		return 1;
	}

	struct points *ps = parse_file(argv[1]);

	if (!ps) {
		printf("Input error\n");
		return 2;
	}
	ps_to_string(ps);
	lagrange(ps);
	ps_free(ps);

	return 0;
}
