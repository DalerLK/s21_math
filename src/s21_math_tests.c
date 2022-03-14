#include <check.h>
#include <string.h>
#include "s21_math.h"

#define DESTROY(VALUE) \
  if (VALUE) {         \
    free(VALUE);       \
  }

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int return_double(double *number, int count) {
    switch (count) {
        case 0:*number = -10.0; break;
        case 1:*number = -1.0; break;
        case 2:*number = -0.5; break;
        case 3:*number = 0.0; break;
        case 4:*number = 0.5; break;
        case 5:*number = 1.0; break;
        case 6:*number = 10.0; break;
        case 7:*number = s21_INFINITY; break;
        case 8:*number = - s21_INFINITY; break;
        case 9:*number = s21_NAN; break;
        case 10:*number = -5.55; break;
        case 11:*number = -0.67; break;
        case 12:*number = 12.34; break;
        case 13:*number = 99.98; break;
        case 14:*number = 0.0001; break;
        case 15:*number = 1000; break;
        default: *number = 0; count = -1;
    }
    return count;
}


START_TEST(s21_math_sin_test) {
  double step = 0.001, step_counter = 0, start = -4 * s21_PI, stop = 4 * s21_PI,
         x;
  x = start;
  while (x < stop + 0.01) {
    long double s21_sin_x, sin_x;
    s21_sin_x = s21_sin(x);
    sin_x = sin(x);
    ck_assert_double_lt(fabs((double)s21_sin_x - sin_x), 0.0000001);
    x = start + s21_PI * step * (long double)step_counter;
    step_counter++;
  }
}
END_TEST

START_TEST(s21_math_cos_test) {
  double step = 0.001, step_counter = 0, start = -4 * s21_PI, stop = 4 * s21_PI,
         x;
  x = start;
  while (x < stop + 0.01) {
    long double s21_cos_x, cos_x;
    s21_cos_x = s21_cos(x);
    cos_x = cos(x);
    ck_assert_double_lt(fabs((double)s21_cos_x - cos_x), 0.0000001);
    x = start + s21_PI * step * (long double)step_counter;
    step_counter++;
  }
}
END_TEST

START_TEST(s21_math_tan_test) {
  double step = 0.01, step_counter = 0, start = -4 * s21_PI, stop = 4 * s21_PI,
         x;
  long double max, prescision;
  x = start;
  while (x < stop + 0.01) {
    long double diff, s21_tan_x, tan_x;
    s21_tan_x = s21_tan(x);
    tan_x = tan(x);
    if (fabs(tan_x) > fabs((double)s21_tan_x)) {
      max = fabs(tan_x);
    } else {
      max = fabs((double)s21_tan_x);
    }
    diff = fabs(tan_x - (double)s21_tan_x);
    if (max > 1) {
      prescision = diff / max;
    } else {
      prescision = fabs((double)s21_tan_x - tan_x);
    }
    if (max < 400000000000000) {
      ck_assert_double_lt(prescision, 0.00000009);
    }
    x = start + s21_PI * step * (long double)step_counter;
    step_counter++;
  }
}
END_TEST

START_TEST(s21_sqrt_test) {
  double step = 0.01, step_counter = 0, start = 0, stop = 1.2, x;
  x = start;
  while (x < stop + 0.01) {
    long double s21_sqrt_x = s21_sqrt(x), diff, max, prescision;
    double sqrt_x = sqrt(x);
    if (fabs(sqrt_x) > fabs((double)s21_sqrt_x)) {
      max = fabs(sqrt_x);
    } else {
      max = fabs((double)s21_sqrt_x);
    }
    diff = fabs(sqrt_x - (double)s21_sqrt_x);
    if (max > 1) {
      prescision = diff / max;
    } else {
      prescision = fabs((double)s21_sqrt_x - sqrt_x);
    }
    ck_assert_double_lt(prescision, 0.00000009);
    x = start + step * (long double)step_counter;
    step_counter++;
  }
}
END_TEST

START_TEST(s21_log_test) {
  double res2 = s21_INFINITY;
  double res3 = s21_NAN;
  int nantest1 = isnan(log(res3));
  int nantest2 = isnan(s21_log(res3));
  ck_assert_double_eq((double)s21_log(res2), log(res2));
  ck_assert_int_eq(nantest1, nantest2);
}
END_TEST

START_TEST(s21_math_log_test) {
    long double step =  0.1, start = -10, stop = 10, x;
    x = start;
    int step_counter = 0;
    while (x < stop + 0.01) {
        long double s21_log_x = s21_log((double)x), max, diff, prescision;
        long double log_x = (long double) log(x);
        if (fabsl(log_x) > fabsl(s21_log_x)) {
            max = fabsl(log_x);
        } else {
            max = fabsl(s21_log_x);
        }
        diff = fabsl(log_x - s21_log_x);
        if (max > 1) {
            prescision = diff/max;
        } else {
            prescision = fabsl(s21_log_x - log_x);
        }
        if (prescision > 0.0000009) {
            ck_assert_double_lt(prescision, 0.00000009);
        }
        x = start + step * (long double) step_counter;
        step_counter++;
    }
}
END_TEST

START_TEST(s21_abs_test) {
  double res1 = 2.1;
  ck_assert_int_eq(s21_abs(res1), abs(res1));
}
END_TEST

START_TEST(s21_fmod_test) {
  double res1 = 2.1;
  double res2 = 2.1345;
  ck_assert_int_eq(s21_fmod(res1, res2), fmod(res1, res2));
}
END_TEST

START_TEST(s21_exp_test) {
  double res1 = 22.1345;
  ck_assert_int_eq(s21_exp(res1), exp(res1));
}
END_TEST

START_TEST(s21_math_exp_test) {
    long double step =  0.01, start = -10, stop = 10, x;
    x = start;
    int step_counter = 0;
    while (x < stop + 0.01) {
        long double s21_exp_x = s21_exp((double)x), max, diff, prescision;
        long double exp_x = (long double) exp((double)x);
        if (fabsl(exp_x) > fabsl(s21_exp_x)) {
            max = fabsl(exp_x);
        } else {
            max = fabsl(s21_exp_x);
        }
        diff = fabsl(exp_x - s21_exp_x);
        if (max > 1) {
            prescision = diff/max;
        } else {
            prescision = fabsl(s21_exp_x - exp_x);
        }
        if (prescision > 0.0000009) {
            ck_assert_double_lt(prescision, 0.00000009);
        }
        x = start + step * (long double) step_counter;
        step_counter++;
    }
}
END_TEST

START_TEST(s21_math_exp_test_stage_2) {
    long double step =  1000, start = 0, stop = 10000000, x;
    x = start;
    int step_counter = 0;
    while (x < stop + 0.01) {
        long double s21_exp_x = s21_exp((double)x), max, diff, prescision;
        long double exp_x = (long double) exp((double)x);
        if (fabsl(exp_x) > fabsl(s21_exp_x)) {
            max = fabsl(exp_x);
        } else {
            max = fabsl(s21_exp_x);
        }
        diff = fabsl(exp_x - s21_exp_x);
        if (max > 1) {
            prescision = diff/max;
        } else {
            prescision = fabsl(s21_exp_x - exp_x);
        }
        if (prescision > 0.0000009) {
            ck_assert_double_lt(prescision, 0.00000009);
        }
        x = start + step * (long double) step_counter;
        step_counter++;
    }
}
END_TEST

START_TEST(s21_pow_test) {
  double res1 = 44.89;
  double res2 = 2.005;
  ck_assert_double_eq_tol(s21_pow(res1, res2), pow(res1, res2), 6);
}
END_TEST

START_TEST(s21_math_pow_test) {
    long double step =  0.5, start = -5, stop = 5, x;
    long double exp_step = 0.5, exp_start = -5, exp_stop = 5, exp;
    x = start;
    exp = exp_start;
    int step_counter = 0, exp_step_counter = 0;
    while (exp < exp_stop + 0.01) {
        while (x < stop + 0.01) {
            long double s21_pow_x = s21_pow((double)x, (double)exp), max, diff, prescision;
            long double pow_x = (long double) pow((double)x, (double)exp);
            if (fabsl(pow_x) > fabsl(s21_pow_x)) {
                max = fabsl(pow_x);
            } else {
                max = fabsl(s21_pow_x);
            }
            diff = fabsl(pow_x - s21_pow_x);
            if (max > 1) {
                prescision = diff/max;
            } else {
                prescision = fabsl(s21_pow_x - pow_x);
            }
            if (prescision > 0.0000009) {
               ck_assert_double_lt(prescision, 0.00000009);
            }
            x = start + step * (long double) step_counter;
            step_counter++;
        }
        exp = exp_start + exp_step * (long double) exp_step_counter;
        exp_step_counter++;
        step_counter = 0;
        x = start;
    }
}
END_TEST

START_TEST(s21_math_pow_test_stage_2) {
    int base_counter = 0, exp_counter = 0;
    while (base_counter > -1) {
        double x;
        base_counter = return_double(&x, base_counter);
        while (exp_counter > -1) {
            double exp;
            exp_counter = return_double(&exp, exp_counter);
            long double s21_pow_x = s21_pow(x, exp), max, diff, prescision;
            long double pow_x = (long double) pow(x, exp);
            if (fabsl(pow_x) > fabsl(s21_pow_x)) {
                max = fabsl(pow_x);
            } else {
                max = fabsl(s21_pow_x);
            }
            diff = fabsl(pow_x - s21_pow_x);
            if (max > 1) {
                prescision = diff/max;
            } else {
                prescision = fabsl(s21_pow_x - pow_x);
            }
            if (prescision > 0.0000009) {
               ck_assert_double_lt(prescision, 0.00000009);
            }
            if (exp_counter != -1) exp_counter++;
        }
        if (base_counter != -1) base_counter++;
        exp_counter = 0;
    }
}
END_TEST



START_TEST(s21_fabs_test) {
  double fabstest1 = 0;
  double fabstest2 = -2.123456;
  double fabstest3 = 4.654321;
  double fabstest4 = 3;
  ck_assert_double_eq(s21_fabs(fabstest1), fabs(fabstest1));
  ck_assert_double_eq(s21_fabs(fabstest2), fabs(fabstest2));
  ck_assert_double_eq(s21_fabs(fabstest3), fabs(fabstest3));
  ck_assert_double_eq(s21_fabs(fabstest4), fabs(fabstest4));
}
END_TEST

START_TEST(s21_ceil_test) {
  double test1 = -1.5;
  double test2 = -5.351;
  double test3 = 0.5;
  double test4 = 0;
  double test5 = s21_NAN;
  double test6 = s21_INFINITY;
  int nantest1 = isnan(ceil(test5));
  int nantest2 = isnan(s21_ceil(test5));
  ck_assert_int_eq(nantest1, nantest2);
  ck_assert_int_eq(s21_ceil(test1), ceil(test1));
  ck_assert_int_eq(s21_ceil(test2), ceil(test2));
  ck_assert_int_eq(s21_ceil(test3), ceil(test3));
  ck_assert_int_eq(s21_ceil(test4), ceil(test4));
  ck_assert_int_eq(s21_ceil(test5), ceil(test5));
  ck_assert_int_eq(s21_ceil(test6), ceil(test6));
}
END_TEST

START_TEST(s21_floor_test) {
  double test1 = -1.5;
  double test2 = -5.351;
  double test3 = 0.5;
  double test4 = 0;
  double test5 = s21_NAN;
  double test6 = s21_INFINITY;
  int nantest1 = isnan(floor(test5));
  int nantest2 = isnan(s21_floor(test5));
  ck_assert_int_eq(nantest1, nantest2);
  ck_assert_int_eq(s21_floor(test1), floor(test1));
  ck_assert_int_eq(s21_floor(test2), floor(test2));
  ck_assert_int_eq(s21_floor(test3), floor(test3));
  ck_assert_int_eq(s21_floor(test4), floor(test4));
  ck_assert_int_eq(s21_floor(test5), floor(test5));
  ck_assert_int_eq(s21_floor(test6), floor(test6));
}
END_TEST

START_TEST(s21_math_atan_test) {
  double step = 0.01, step_counter = 0, start = -4, stop = 4, x;
  x = start;
  while (x < stop + 0.01) {
    long double s21_atan_x = s21_atan(x), diff, max, prescision;
    double atan_x = atan(x);
    if (fabs(atan_x) > fabs((double)s21_atan_x)) {
      max = fabs(atan_x);
    } else {
      max = fabs((double)s21_atan_x);
    }
    diff = fabs(atan_x - (double)s21_atan_x);
    if (max > 1) {
      prescision = diff / max;
    } else {
      prescision = fabs((double)s21_atan_x - atan_x);
    }
    ck_assert_double_lt(prescision, 0.00000009);
    x = start + step * (long double)step_counter;
    step_counter++;
  }
}
END_TEST

START_TEST(s21_math_acos_test) {
  double step = 0.01, step_counter = 0, start = -1.2, stop = 1.2, x;
  x = start;
  while (x < stop + 0.01) {
    long double s21_acos_x = s21_acos(x), diff, max, prescision;
    double acos_x = acos(x);
    if (fabs(acos_x) > fabs((double)s21_acos_x)) {
      max = fabs(acos_x);
    } else {
      max = fabs((double)s21_acos_x);
    }
    diff = fabs(acos_x - (double)s21_acos_x);
    if (max > 1) {
      prescision = diff / max;
    } else {
      prescision = fabs((double)s21_acos_x - acos_x);
    }
    if (isnan(acos_x)) {
      ck_assert_ldouble_nan(s21_acos_x);
    } else {
      ck_assert_double_lt(prescision, 0.00000009);
    }
    x = start + step * (long double)step_counter;
    step_counter++;
  }
}
END_TEST

START_TEST(s21_math_asin_test) {
  double step = 0.01, step_counter = 0, start = -1.2, stop = 1.2, x;
  x = start;
  while (x < stop + 0.01) {
    long double s21_asin_x = s21_asin(x), diff, max, prescision;
    double asin_x = asin(x);
    if (fabs(asin_x) > fabs((double)s21_asin_x)) {
      max = fabs(asin_x);
    } else {
      max = fabs((double)s21_asin_x);
    }
    diff = fabs(asin_x - (double)s21_asin_x);
    if (max > 1) {
      prescision = diff / max;
    } else {
      prescision = fabs((double)s21_asin_x - asin_x);
    }
    if (isnan(asin_x)) {
      ck_assert_ldouble_nan(s21_asin_x);
    } else {
      ck_assert_double_lt(prescision, 0.00000009);
    }
    x = start + step * (long double)step_counter;
    step_counter++;
  }
}
END_TEST

START_TEST(s21_math_inf_test) {
  char s21_str[10] = " ", origin_str[10] = " ";
  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_sin(s21_INFINITY));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)sin(s21_INFINITY));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_sin(-s21_INFINITY));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)sin(-s21_INFINITY));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_cos(s21_INFINITY));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)cos(s21_INFINITY));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_cos(-s21_INFINITY));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)cos(-s21_INFINITY));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_tan(s21_INFINITY));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)tan(s21_INFINITY));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_tan(-s21_INFINITY));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)tan(-s21_INFINITY));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_atan(s21_INFINITY));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)atan(s21_INFINITY));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_atan(-s21_INFINITY));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)atan(-s21_INFINITY));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_asin(s21_INFINITY));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)asin(s21_INFINITY));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_asin(-s21_INFINITY));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)asin(-s21_INFINITY));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_acos(s21_INFINITY));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)acos(s21_INFINITY));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_acos(-s21_INFINITY));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)acos(-s21_INFINITY));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_sqrt(s21_INFINITY));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)sqrt(s21_INFINITY));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_sqrt(-s21_INFINITY));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)sqrt(-s21_INFINITY));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_log(s21_INFINITY));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)log(s21_INFINITY));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_log(-s21_INFINITY));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)log(-s21_INFINITY));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_exp(s21_INFINITY));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)exp(s21_INFINITY));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_exp(-s21_INFINITY));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)exp(-s21_INFINITY));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_pow(s21_INFINITY, 1));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)pow(s21_INFINITY, 1));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_pow(-s21_INFINITY, 1));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)pow(-s21_INFINITY, 1));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_pow(1, s21_INFINITY));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)pow(1, s21_INFINITY));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_pow(1, -s21_INFINITY));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)pow(1, -s21_INFINITY));
  ck_assert_str_eq(s21_str, origin_str);
}
END_TEST

START_TEST(s21_math_nan_test) {
  char s21_str[10] = " ", origin_str[10] = " ";

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_sin(s21_NAN));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)sin(s21_NAN));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_sin(-s21_NAN));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)sin(-s21_NAN));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_cos(s21_NAN));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)cos(s21_NAN));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_cos(-s21_NAN));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)cos(-s21_NAN));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_tan(s21_NAN));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)tan(s21_NAN));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_tan(-s21_NAN));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)tan(-s21_NAN));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_atan(s21_NAN));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)atan(s21_NAN));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_atan(-s21_NAN));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)atan(-s21_NAN));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_asin(s21_NAN));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)asin(s21_NAN));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_asin(-s21_NAN));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)asin(-s21_NAN));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_acos(s21_NAN));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)acos(s21_NAN));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_acos(-s21_NAN));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)acos(-s21_NAN));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_sqrt(s21_NAN));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)sqrt(s21_NAN));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_sqrt(-s21_NAN));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)sqrt(-s21_NAN));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_log(s21_NAN));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)log(s21_NAN));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_log(-s21_NAN));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)log(-s21_NAN));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_exp(s21_NAN));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)exp(s21_NAN));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_exp(-s21_NAN));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)exp(-s21_NAN));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_pow(s21_NAN, 1));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)pow(s21_NAN, 1));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_pow(-s21_NAN, 1));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)pow(-s21_NAN, 1));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_pow(1, s21_NAN));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)pow(1, s21_NAN));
  ck_assert_str_eq(s21_str, origin_str);

  snprintf(s21_str, sizeof(s21_str), " ");
  snprintf(origin_str, sizeof(origin_str), " ");
  snprintf(s21_str, sizeof(s21_str), "%f", (double)s21_pow(1, -s21_NAN));
  snprintf(origin_str, sizeof(origin_str), "%f", (double)pow(1, -s21_NAN));
  ck_assert_str_eq(s21_str, origin_str);
}
END_TEST

int main(void) {
  Suite *s1 = suite_create("s21_math.h");
  SRunner *sr = srunner_create(s1);
  int nf;

  TCase *tc1_1 = tcase_create("s21_sin");
  tcase_add_test(tc1_1, s21_math_sin_test);
  suite_add_tcase(s1, tc1_1);

  TCase *tc2_1 = tcase_create("s21_cos");
  tcase_add_test(tc2_1, s21_math_cos_test);
  suite_add_tcase(s1, tc2_1);
  TCase *tc3_1 = tcase_create("s21_tan");
  tcase_add_test(tc3_1, s21_math_tan_test);
  suite_add_tcase(s1, tc3_1);
  TCase *tc6_1 = tcase_create("s21_sqrt");
  tcase_add_test(tc6_1, s21_sqrt_test);
  suite_add_tcase(s1, tc6_1);
  TCase *tc7_1 = tcase_create("s21_log");
  tcase_add_test(tc7_1, s21_log_test);
  suite_add_tcase(s1, tc7_1);
  TCase *tc8_1 = tcase_create("s21_abs");
  tcase_add_test(tc8_1, s21_abs_test);
  suite_add_tcase(s1, tc8_1);
  TCase *tc9_1 = tcase_create("s21_fmod");
  tcase_add_test(tc9_1, s21_fmod_test);
  suite_add_tcase(s1, tc9_1);
  TCase *tc10_1 = tcase_create("s21_exp");
  tcase_add_test(tc10_1, s21_exp_test);
  suite_add_tcase(s1, tc10_1);
  TCase *tc11_1 = tcase_create("s21_pow");
  tcase_add_test(tc11_1, s21_pow_test);
  suite_add_tcase(s1, tc11_1);
  TCase *tc12_1 = tcase_create("s21_fabs");
  tcase_add_test(tc12_1, s21_fabs_test);
  suite_add_tcase(s1, tc12_1);
  TCase *tc13_1 = tcase_create("s21_ceil");
  tcase_add_test(tc13_1, s21_ceil_test);
  suite_add_tcase(s1, tc13_1);
  TCase *tc14_1 = tcase_create("s21_floor");
  tcase_add_test(tc14_1, s21_floor_test);
  suite_add_tcase(s1, tc14_1);
  TCase *tc15_1 = tcase_create("s21_asin");
  tcase_add_test(tc15_1, s21_math_asin_test);
  suite_add_tcase(s1, tc15_1);
  TCase *tc16_1 = tcase_create("s21_acos");
  tcase_add_test(tc16_1, s21_math_acos_test);
  suite_add_tcase(s1, tc16_1);
  TCase *tc17_1 = tcase_create("s21_atan");
  tcase_add_test(tc17_1, s21_math_atan_test);
  suite_add_tcase(s1, tc17_1);
  TCase *tc18_1 = tcase_create("s21_inf");
  tcase_add_test(tc18_1, s21_math_inf_test);
  suite_add_tcase(s1, tc18_1);
  TCase *tc19_1 = tcase_create("s21_inf");
  tcase_add_test(tc19_1, s21_math_nan_test);
  suite_add_tcase(s1, tc19_1);
  TCase *tc20_1 = tcase_create("s21_math_log");
  tcase_add_test(tc20_1, s21_math_log_test);
  suite_add_tcase(s1, tc20_1);
  TCase *tc21_1 = tcase_create("s21_math_exp");
  tcase_add_test(tc21_1, s21_math_exp_test);
  suite_add_tcase(s1, tc21_1);
  TCase *tc22_1 = tcase_create("s21_math_exp_stage_2");
  tcase_add_test(tc22_1, s21_math_exp_test_stage_2);
  suite_add_tcase(s1, tc22_1);
  TCase *tc23_1 = tcase_create("s21_math_pow_");
  tcase_add_test(tc23_1, s21_math_pow_test);
  suite_add_tcase(s1, tc23_1);
  TCase *tc24_1 = tcase_create("s21_math_pow_stage_2");
  tcase_add_test(tc24_1, s21_math_pow_test_stage_2);
  suite_add_tcase(s1, tc24_1);
  srunner_run_all(sr, CK_VERBOSE);

  nf = srunner_ntests_failed(sr);
  srunner_free(sr);

  return nf == 0 ? 0 : 1;
}
