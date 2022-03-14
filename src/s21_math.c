#include "s21_math.h"

long int s21_abs(int x) { return x < 0 ? x *= -1 : x; }

long double s21_fabs(double x) {
  long double res = x;
  if (res < 0) res *= -1;
  return res;
}

long double s21_fmod(double x, double y) {
  return (x && y) || y != 0 ? (x / y - (int)(x / y)) * y : s21_NAN;
}

long double s21_exp(double x) {
  double val = 1.0, temp = 1.0;
  double iter = 1.0;
  if (x != x) {
    val = -s21_NAN;
  } else if (x == s21_INFINITY) {
    val = x;
  } else if (x == -s21_INFINITY) {
    val = 0;
  } else if (x > 709) {
    val = s21_INFINITY;
  } else if (x < -709) {
    val = 0;
  } else if (x != 0) {
    if (x > 0) {
      for (; temp >= pr; ++iter) {
        temp *= x / iter;
        val += temp;
      }
    } else {
      x *= -1;
      for (; temp >= pr; ++iter) {
        temp *= x / iter;
        val += temp;
      }
      val = 1 / val;
    }
  } else {
    val = 1;
  }
  return val;
}

long double s21_pow(double base, double exp) {
  long double number;
  if (base < 0) {
    if ((long int)exp == exp) {
      if (exp > 0) {
        number = base;
        for (long int i = 0; i < (long int)exp - 1; i++) {
          number *= base;
        }
      } else if (exp == 0) {
        number = 1;
      } else {
        number = 1 / base;
        for (long int i = 0; i < (long int)exp * (-1) - 1; i++) {
          number /= base;
        }
      }
    } else {
      if (exp == -s21_INFINITY || exp == s21_INFINITY) {
        if (base * (-1) < 1) {
          number = 0;
        } else if (base * (-1) == 1) {
          number = 1;
        } else {
          if (exp == -s21_INFINITY) {
            number = 0;
          } else {
            number = s21_INFINITY;
          }
        }
      } else {
        number = -s21_NAN;
      }
    }
  } else if (base == 0) {
    if (exp == 0) {
      number = 1;
    } else {
      number = 0;
    }
  } else if (base == 1) {
    number = 1;
  } else {
    if ((long int)exp == exp) {
      if (exp > 0) {
        number = base;
        for (long int i = 0; i < (long int)exp - 1; i++) {
          number *= base;
        }
      } else if (exp == 0) {
        number = 1;
      } else {
        number = 1 / base;
        for (long int i = 0; i < (long int)exp * (-1) - 1; i++) {
          number /= base;
        }
      }
    } else {
      number = s21_exp(exp * (double)s21_log(base));
    }
  }
  return number;
}

long double s21_sqrt(double x) {
  long double sum = x;
  int flag_sqrt = 1;
  if (sum < 0) {
    sum = s21_NAN;
    flag_sqrt--;
  }
  if (sum == 0) {
    flag_sqrt = 0;
  }
  if (x != x) {
    flag_sqrt = 0;
  }
  if (flag_sqrt == 1) {
    if (sum > 1) {
      double y = 1;
      double e = 0.000000001;
      while (sum - y >= e) {
        sum = (sum + y) / 2;
        y = x / sum;
      }
    } else if (sum < 1) {
      double var = x - 1;
      double q = 1;
      sum = 1;
      for (int n = 1; n <= 1000; n++) {
        q *= (-1.0) * (2 * n - 1) * 2 * n * var / (n * n * 4);
        sum += q / (1.0 - 2 * n);
      }
    }
  }
  return sum;
}

long double s21_log(double x) {
  long double iter = 0;
  if (x != x) {
    iter = -s21_NAN;
  } else if (x == s21_INFINITY) {
    iter = x;
  } else if (x == -s21_INFINITY || x < 0) {
    iter = -s21_NAN;
  } else if (x == 0) {
    iter = -s21_INFINITY;
  } else {
    if (x >= 0.1 && x <= 10) {
      if (x < 0.3 || x > 0.6) {
        x = (x / (x - 1));
        for (int i = 1; i < 1000; i++) {
          iter += 1.0 / (i * simple_pow(x, i));
        }
      } else {
        iter = s21_log(x * 2) - s21_log(2);
      }
    } else {
      int exp = 0;
      if (x > 10) {
        while (x > 10) {
          x /= 10;
          exp++;
        }
      } else {
        while (x < 0.1) {
          x *= 10;
          exp--;
        }
      }
      iter = s21_log(x) + exp * E_10;
    }
  }
  return iter;
}

long double simple_pow(double x, int y) {
  double result = 1;
  for (int i = 0; i < y; ++i) result *= x;
  return result;
}

long double s21_sin(double x) {
  long double y;
  if (x != x) {
    y = -s21_NAN;
  } else if (x == s21_INFINITY || x == -s21_INFINITY) {
    y = -s21_NAN;
  } else {
    int negative = 1;
    if (x < 0) {
      negative = -1;
      x *= -1;
    }
    while (x > s21_PI_2) {
      x -= s21_PI;
      negative *= -1;
    }
    if (x <= s21_PI_3 && x >= -s21_PI_3) {
      long double P1 = 1.5707963, P2 = 0.64596410, P3 = 0.079692626;
      long double P4 = 0.0046817541, P5 = 0.00016044118, P6 = 0.0000035988432;
      long double d_x3, d_x5, d_x7, d_x9, d_x11;
      x = x / s21_PI_2;
      d_x3 = x * x * x;
      d_x5 = d_x3 * x * x;
      d_x7 = d_x5 * x * x;
      d_x9 = d_x7 * x * x;
      d_x11 = d_x9 * x * x;
      y = P1 * x - P2 * d_x3 + P3 * d_x5 - P4 * d_x7 + P5 * d_x9 + P6 * d_x11;
    } else {
      if (x > s21_PI_3) {
        y = s21_sin(x - s21_PI_6) * s21_SQRT_3_D2 + s21_cos(x - s21_PI_6) * 0.5;
      } else {
        y = (s21_sin(x + s21_PI_6) * s21_SQRT_3_D2 -
             s21_cos(x + s21_PI_6) * 0.5);
      }
    }
    y *= (long double)negative;
    if (y < -1.0) y = -0.99999999999999999;
    if (y > 1.0) y = 0.9999999999999999990;
  }
  return y;
}

long double s21_tan(double x) {
  long double y = 0;
  if (x != x) {
    y = -s21_NAN;
  } else if (x == s21_INFINITY || x == -s21_INFINITY) {
    y = -s21_NAN;
  } else {
    long double sin_x = s21_sin(x), cos_x = s21_cos(x);
    if (cos_x == 0) {
      if (sin_x > 0) {
        y = s21_INFINITY;
      } else if (sin_x < 0) {
        y = -s21_INFINITY;
      }
    } else {
      y = sin_x / cos_x;
    }
  }
  return y;
}

long double s21_cos(double x) {
  long double y;
  if (x != x) {
    y = -s21_NAN;
  } else if (x == s21_INFINITY || x == -s21_INFINITY) {
    y = -s21_NAN;
  } else {
    if (x < 0) {
      x *= -1;
    }
    while (x >= 2 * s21_PI) {
      x -= 2 * s21_PI;
    }
    x = s21_PI_2 - x;
    y = s21_sin(x);
  }
  return y;
}

long double s21_ceil(double x) {
  double arg = x;
  int num = arg;
  long double res;
  if (x == s21_INFINITY || x == -s21_INFINITY) {
    res = x;
  } else if (arg == num) {
    res = num;
  } else if (arg < (long double)num) {
    res = num;
  } else if (arg > (long double)num) {
    res = num + 1;
  }

  return x != x ? s21_NAN : res;
}

long double s21_floor(double x) {
  double arg = x;
  int num = arg;
  long double res;
  if (x == s21_INFINITY || x == -s21_INFINITY) {
    res = x;
  } else if (arg == num) {
    res = num;
  } else if (arg > (long double)num) {
    res = num;
  } else if (arg < (long double)num) {
    res = num - 1;
  }

  return x != x ? s21_NAN : res;
}

long double s21_asin(double x) {
  long double y;
  if (x != x || x > 1.0 || x < -1.0) {
    y = s21_NAN;
  } else if (x == 1) {
    y = s21_PI_2;
  } else if (x == -1) {
    y = -s21_PI_2;
  } else {
    y = s21_atan(x / s21_sqrt(1 - x * x));
  }
  return y;
}

long double s21_acos(double x) {
  long double y;
  if (x != x || x > 1.0 || x < -1.0) {
    y = s21_NAN;
  } else if (x == 1) {
    y = 0;
  } else if (x == -1) {
    y = s21_PI;
  } else {
    y = s21_PI_2 - s21_asin(x);
  }
  return y;
}

long double s21_atan(double x) {
  long double y;
  if (x != x) {
    y = -s21_NAN;
  } else if (x == s21_INFINITY) {
    y = s21_PI_2;
  } else if (x == -s21_INFINITY) {
    y = -s21_PI_2;
  } else {
    long double dx, s21_abss;
    int n = 1, flag = 0, flag_negative = 0;
    if (x < 0) {
      flag_negative = 1;
      x *= -1;
    }
    if (x > 1) {
      x = 1 / x;
      flag = 1;
    }
    dx = x;
    y = x;
    s21_abss = x;
    if (x != 1) {
      while (s21_abss > 0.000000000000001) {
        ++n;
        dx *= -x * x * (2 * n - 3) / (2 * n - 1);
        y += dx;
        if (dx < 0) {
          s21_abss *= -dx;
        } else {
          s21_abss = dx;
        }
      }
    } else {
      y = s21_PI_4;
    }
    if (flag == 1) {
      y = s21_PI_2 - y;
    }
    if (flag_negative == 1) {
      y *= -1;
    }
  }
  return y;
}
