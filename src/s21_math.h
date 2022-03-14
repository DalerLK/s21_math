#ifndef SRC_S21_MATH_H_
#define SRC_S21_MATH_H_

#include <stdio.h>
#include <stdlib.h>

#define s21_PI 3.141592653589793116
#define s21_PI_2 1.570796326794896558
#define s21_PI_3 1.047197551196597853
#define s21_PI_6 0.523598775598298926
#define s21_PI_4 0.785398163397448309
#define SIN_ADD_MINUS -0.00000000000000006123232794L
#define SIN_ADD_PLUS 0.00000000000000006123234309L

#define s21_SQRT_2_D2 0.7071067811865475243818940365159164684L
#define s21_SQRT_3_D2 0.8660254037844386467868626477972782140L
#define s21_SQRT_3_D1 1.7320508075688772935737252955945564281L

#define s21_e 2.718281828459045091
#define logE 1
#define pr 1e-16
#define s21_INFINITY 1.0 / 0.0
#define s21_NAN 0.0 / 0.0
#define E_10 2.302585092994046

long double s21_sin(double x);
long double s21_cos(double x);
long double s21_tan(double x);
long int s21_abs(int x);
long double s21_ln(double x);
long double s21_log(double x);
long double simple_pow(double x, int y);
long double s21_pow(double base, double exp);
long double s21_exp(double x);
long double s21_fmod(double x, double y);
long double s21_sqrt(double x);
long double s21_ceil(double x);
long double s21_floor(double x);
long double s21_fabs(double x);
long double s21_atan(double x);
long double s21_acos(double x);
long double s21_asin(double x);

#endif  // SRC_S21_MATH_H_
