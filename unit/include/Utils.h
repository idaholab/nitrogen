#ifndef UTILS_H
#define UTILS_H

// Relative tolerance to be used when comparing to a value from external fluid
// property packages, which might be using different underlying models
#define REL_TOL_EXTERNAL_VALUE 1e-3

// Relative tolerance to be used when comparing to a value computed directly
// from the code at an earlier date, to ensure that implementations are not
// inadvertently changed
#define REL_TOL_SAVED_VALUE 1e-12

// Relative tolerance to be used for consistency checks - computing properties
// in different ways at the same state
#define REL_TOL_CONSISTENCY 1e-10

// Relative tolerance to be used when comparing a derivative value to a finite
// difference approximation
#define REL_TOL_DERIVATIVE 1e-6

// Macro for computing relative error
#define REL_TEST(value, ref_value, tol)                                                            \
  EXPECT_LE(std::abs((value) - (ref_value)) / (ref_value), (tol))

// Macro for computing absolute error
#define ABS_TEST(value, ref_value, tol) EXPECT_LE(std::abs(((value) - (ref_value))), (tol))

// Relative perturbation for derivative tests
#define REL_PERTURBATION 1e-6

// Macro for performing a derivative test
#define DERIV_TEST(f, a, b, tol)                                                                   \
  {                                                                                                \
    const Real da = REL_PERTURBATION * a;                                                          \
    const Real db = REL_PERTURBATION * b;                                                          \
    const Real df_da_fd = (f(a + da, b) - f(a - da, b)) / (2 * da);                                \
    const Real df_db_fd = (f(a, b + db) - f(a, b - db)) / (2 * db);                                \
    Real f_value, df_da, df_db;                                                                    \
    f(a, b, f_value, df_da, df_db);                                                                \
    REL_TEST(f(a, b), f_value, REL_TOL_CONSISTENCY);                                               \
    REL_TEST(df_da, df_da_fd, tol);                                                                \
    REL_TEST(df_db, df_db_fd, tol);                                                                \
  }

#endif // UTILS_H
