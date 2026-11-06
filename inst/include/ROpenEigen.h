// ROpenEigen.h —— 公共 C++ 头文件
//
// 装了本包之后，任何写了 LinkingTo: ROpenEigen 的包都可以
//
//     #include <ROpenEigen.h>
//
// 直接拿到 Rcpp 与 Eigen 5.0.1，不必再单独依赖 RcppEigen。Eigen 是
// header-only 的，所以除了头文件之外什么都不需要链接。
//
// 注意 R_NO_REMAP：R 的 C API 把 length、error、vector 之类定义成了宏，
// 和 Eigen 的同名成员直接冲突。在包含 Rcpp 之前定义它，R 的符号就只剩
// Rf_ 前缀的那一套，Eigen 才编得过。这一点对下游包同样成立，所以放在
// 公共头里而不是某个 .cpp 里。
//
// 这里只引 Dense 和 Sparse。unsupported/ 整棵子树也随包安装好了
// （SpecialFunctions、MatrixFunctions、FFT、Polynomials、Splines、
// NonLinearOptimization、AutoDiff、NNLS 等），需要时自己
//     #include <unsupported/Eigen/MatrixFunctions>
// 即可——把它们塞进公共头会让每个下游翻译单元多编几十秒，尤其是
// NonLinearOptimization。

#ifndef ROpenEigen__ROpenEigen__h
#define ROpenEigen__ROpenEigen__h

#ifndef R_NO_REMAP
# define R_NO_REMAP
#endif

#include <Rcpp.h>

#include <Eigen/Dense>    // Core, LU, Cholesky, QR, SVD, Geometry, Eigenvalues
#include <Eigen/Sparse>

// 版本自检：下游包若误装了别的 Eigen（比如同时 LinkingTo RcppEigen），
// 编译器会先找到那一个，而本包的封装假设 5.0.1 的语义。早报错比运行期
// 得到另一个库的结果好。
#include <Eigen/Version>
#if !EIGEN_VERSION_AT_LEAST(5, 0, 0)
# error "ROpenEigen requires Eigen 5.0.1; an older Eigen was found on the include path"
#endif

#endif
