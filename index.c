#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 1000

int compare(const void *a, const void *b) {
    const double *da = (const double *)a;
    const double *db = (const double *)b;
    return (*da > *db) - (*da < *db);
}

double norm_rasp(double mean, double stddev) {
    double u1, u2, v1, v2, s;
    do {
        u1 = ((double)rand() / RAND_MAX) * 2 - 1;
        u2 = ((double)rand() / RAND_MAX) * 2 - 1;
        s = u1 * u1 + u2 * u2;
    } while (s >= 1 || s == 0);

    v1 = sqrt(-2 * log(s) / s) * u1;
    v2 = sqrt(-2 * log(s) / s) * u2;

    return mean + stddev * v1;
}

double norm_cdf(double x, double mean, double stddev) {
    return 0.5 * erfc(-(x - mean) / (stddev * sqrt(2)));
}

double chisqinv(double p, int df) {
    double epsilon = 0.0001; // Погрешность
    double x_low = 0.0;      // Нижняя граница диапазона
    double x_high = df;      // Верхняя граница диапазона
    double x_mid = 0.0;      // Середина диапазона
    double p_mid = 0.0;      // Значение функции распределения в середине диапазона

    while (x_high - x_low > epsilon) {
        x_mid = (x_low + x_high) / 2.0;
        p_mid = 0.5 * (1 + erf((x_mid - df) / (sqrt(2 * df))));

        if (p_mid < p) {
            x_low = x_mid;
        } else {
            x_high = x_mid;
        }
    }

    return x_mid;
}

int main() {
  double m = 3.0;
  double sigma = 1.0;
  double alfa = 0.01;
  double sample[N];

  srand(time(NULL));
  for (int i = 0; i < N; i++) {
    sample[i] = m + sigma * sqrt(-2 * log((double)rand() / RAND_MAX)) * cos(2 * M_PI * ((double)rand() / RAND_MAX));
  }
  qsort(sample, N, sizeof(double), compare);

  double border = 0.0;
  int count = 0;
  double pockets[15] = {0.0}; //границы карманов
  int pockets_sample[15] = {0}; //содержимое карманов
  int pockets_count = 0; //количество карманов

  for (int i = 0; i < N; i++) {
    if (sample[i] >= border) {
      if (count >= 5) {
        pockets_sample[pockets_count] = count;
        pockets[pockets_count] = border;
        count = 0;
        pockets_count++;
      }
      border += 0.6;
      if (border > 6.6) {
        pockets_sample[pockets_count] = count;
        pockets[pockets_count] = 6.6;
        pockets_count++;
        break;
      }
      i--;
      continue;
    }    
    count++;
  }

  if (count < 5) {
    pockets_sample[pockets_count - 1] += count;
  } else {
    pockets_sample[pockets_count] = count;
    pockets[pockets_count] = 6.6;
    pockets_count++;
  }

  int rit = 0;
  for (int i = 0; i <= pockets_count; i++) {
    rit += pockets_sample[i];
  }

  pockets[pockets_count - 1] = 100;
  double *ri_sample = malloc(pockets_count * sizeof(double));
  ri_sample[0] = (norm_cdf(pockets[0], m, sigma) - norm_cdf(-100.0, m, sigma)) * N;

  for (int i = 1; i < pockets_count; i++) {
    ri_sample[i] = (norm_cdf(pockets[i], m, sigma) - norm_cdf(pockets[i-1], m, sigma)) * N;
  }

  double ri = 0.0;
  for (int i = 0; i < pockets_count; i++) {
    ri += ri_sample[i];
  }

  double S = (ri - rit) * (ri - rit) / ri;
  double chi2_inv = chisqinv(alfa, pockets_count - 1);

  if (S < chi2_inv) {
    printf("true\n");
  } else {
    printf("false\n");
  }

  free(ri_sample);
  return 0;
}

