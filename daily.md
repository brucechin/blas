#### **第一日**

openblas api文档阅读与测试高频使用的api与裸写的矩阵操作之间速度差距

##### level1 vector * vector

| 函数                 | 解释                                    | size=5W时运行10W次耗时(us) |
| -------------------- | --------------------------------------- | -------------------------- |
| 自己实现的float dot  | 输入为float，输出为float                | 6393417                    |
| 自己实现的double dot | 输入为double，输出为double              | 14511149                   |
| sdot                 | 输入为float，输出为float                | 547730                     |
| ddot                 | 输入为double，输出为double              | 1082742                    |
| dsdot                | 输入float，转为double计算，输出为double | 3209029                    |
| sdsdot               | 输入float，转为double计算，输出为float  | 3273407                    |

1. 观察到vector size从2K->100K变化过程中计算耗时随size线性增长，详细数据不再列出
2. dsdot和sdsdot相对较慢怀疑是float转double的过程额外耗时较多
3. 一般情况下openblas版本要比自己实现的baseline快10-20倍

**level2 vector * matrix**

**level3 matrix * matrix**

直接测的matrix*matrix

| 测试函数\矩阵规模                    | 512*512（耗时单位s） |
| ------------------------------------ | -------------------- |
| 自己实现的未优化版本double精度矩阵乘 | 44.53                |
| cblas_dgemm(输入double，输入double)  | 0.72                 |

1. 没测单精度下的结果，预计类似，耗时与矩阵规模的三次方成正比，具体数据不再上传
2. openblas版本比自己实现的baseline快20-80倍

### 第二日

看完了matrix.java的几个类，用C写了Matrix.h, LogicMatrix.h, matrixFactory.h.

存在疑虑是，在基础类中是否可以用到OpenBLAS库来优化？类似for循环遍历getElement/setElement的操作会很慢的吧？

### 第三日

看了一下[blislab](https://github.com/flame/blislab)优化矩阵计算的一些trick，包括：

> 1. cache friendly pointer access
> 2. loop unrolling
> 3. register variable
> 4. parallelizing with OpenMP

有一个设计上的疑虑在于，矩阵的data是用row-major还是column-major来存？？Eigen有以下的回答：

> which storage order should you use in your program? There is no simple answer to this question; it depends on your application. Here are some points to keep in mind:
>
> - Your users may expect you to use a specific storage order. Alternatively, you may use other libraries than [Eigen](https://eigen.tuxfamily.org/dox-devel/namespaceEigen.html), and these other libraries may expect a certain storage order. In these cases it may be easiest and fastest to use this storage order in your whole program.
> - Algorithms that traverse a matrix row by row will go faster when the matrix is stored in row-major order because of better data locality. Similarly, column-by-column traversal is faster for column-major matrices. It may be worthwhile to experiment a bit to find out what is faster for your particular application.
> - The default in [Eigen](https://eigen.tuxfamily.org/dox-devel/namespaceEigen.html) is column-major. Naturally, most of the development and testing of the [Eigen](https://eigen.tuxfamily.org/dox-devel/namespaceEigen.html) library is thus done with column-major matrices. This means that, even though we aim to support column-major and row-major storage orders transparently, the [Eigen](https://eigen.tuxfamily.org/dox-devel/namespaceEigen.html) library may well work best with column-major matrices.

最后跟mentor聊了一下，定的是row-major。

### 第四日

上午编译了openblas在Windows平台上，但只能编译出lib文件，dll文件出不来很奇怪，后来找了编译好的直接用。。

发现openblas缺乏elementwise操作的各种接口，如两个矩阵每个元素一一比较，round，floor，sqrt，pow等等，但考虑到Intel MKL和openblas接口基本是完全相同的（要符合BLAS约定？），不可能两个库同时使用，暂考虑使用性能稍差一些的Intel MKL

下午实现了如下接口：

- [x] add, sub, div, mul, matrixMul

### 第五日

今日确定了一些设计上的问题：

1. value设为protected属性，只能由MatrixCalculator类访问
2. elementwise操作时不用get/setElement接口，直接指针访存

今日实现了以下接口：

- [x] max, min, bigger, smaller, equal, between
- [x] and, or, not, condition
- [x] rank, round, floor, abs, minus, sqrt, log, exp, sign, inverse, signedpow

### 第六日

讨论了时序数据处理的一些接口中传入参数的形式，java代码中是传入一个一维数组，我最开始在C中实现的是传一个double指针，但是发现难以获取指针指向内存空间的大小，后来把传入参数统一为Matrix，时序数据这些接口传入一维的Matrix

今日实现以下接口

- [x] shift, delay, delta, ratio, sum, product
- [ ] tsMax, tsMin, tsArgmax/min, tsRank, tsMean, tsStd, tsSkewness, tsKurtosis, tsCov, tsCorr, tsCountNaN, tsCountTrue
- [ ] decayLinear, decayExponential, smoothByDecayLinear
- [ ] activate, normalize, neutralize, unify, evalValidPct, evalAbsSum, evalMean, evalVariance, evalCorr, evalCovariance
- [ ] Det, Inverse, inv, treat, diag, inverseDiag, evalBeta
- [x] summaryMean/Variance/Skewness/Kurtosis/Covariance/Corr/Sum/Gini

### 7.15

1. 已基本实现所有接口，但超过一半的是没有调用外部库实现的，可能效率不高？
2. 不改变输入参数Matrix的时候是pass by value，为了降低调用开销是否改成pass by reference to const?(不用调用Matrix的拷贝构造和析构函数等)
3. 开始准备写测试代码，原java文件竟然没有配套测试代码？？？？

### 7.16

1. 测试需要生成size不同的matrix然后把整个对象以二进制存为文件，记得一些corner case的生成

2. openblas windows安装遇到的坑：
   1. mingw32下编译会报缺少某obj文件的error
   2. mingw64下编译一切顺利，之前用visual studio编译cmake生成的项目没有生成libopenblas
   
3. vscode c++ compile and run配置好了，但是有些时候不能自动跳转错误对应行。

4. debug遇到的坑：

   1. ###### Static function: a storage class may not be specified here，static的函数在头文件里声明了之后，在cpp文件里不需要了。
   
   2. 因为类里重载过max min floor abs round，在类里使用std的这些函数的时候没加std::导致找不到对应的函数了。。
   
   3. 把所有的函数参数都变成matrix的指针了，传值的话拷贝构造额外开销怕是要崩
   
   4. const static成员变量不能in class initialized
   
   5. redefinition of class xxx :  在h文件中直接实现接口的话，要在文件开头注明
   
      1. > #ifndef XXX_H #define XXX_H 
         >
         > 代码区
         >
         > #endif
   
   6. matrix类里应该是重载and or not，不能直接定义一个函数，因为他们是保留关键字

### 7.17

1. 开始配置google test环境写单元测试，include的时候遇到的坑是要把include文件夹内的东西原封不动挪到系统的include文件夹内，不然可能会破坏原来的结构
2. vscode下gtest integration失败，不过在visual studio里面很顺利，目前规划在VS里写单元测试，先测试正确性，再在开发机上测速度。在visual studio出的问题又是，没法把openblas库的lib文件link到项目里，出现了经典的LINK2019的bug，找不到MatrixCalculator里调用的openblas函数的定义。后来在Windows下面开了linux subsystem，装好了环境打算先在这上面搞，因为开发机上没权限配环境很麻烦。
3. 此外由于本机内存大小的限制，测试中matrix最大size不打算超过5k*5k
4. 需要系统学习一下g++, makefile，编译原理等
5. gtest一些总结
   1. ASSERT_* 系列的断言，当检查点失败时，退出当前函数（并不会停止后续的TEST()）
   2. EXPECT_* 系列的断言，当检查点失败时，继续往下执行。

### 7.18

1. 写完了主要函数的单元测试，基本都能跑出来，除了Det和Rank会core dump，之后去修一下。速度慢的主要是timeseries的函数，因为是n^3的复杂度，之后打算用sliding window的trick改一改
2. 类似tsmax, tsmin, tssum之类的用sliding window改了，tsvar, tscorr, tscov这种还是算了吧。。
3. 发现老版本的tsmax tsmin都被我抄错了，有毒啊，还是得每个函数以较小的输入时打印出来结果比对来找出实现中的小错误。
4. 老版本的tsargmax和tsargmin写错了，记录的位置应该是j-k而不是k。同时count consecutive true也不太他对，java版本是遇到一个false就停止了，而且是倒着数的，记得汇报
5. 总结一下，目前除了Det，Rank会core dump以外，其他函数都可以运行，但结果正确性无法保证，ts有optimize过的几个函数比较过正确性了，一些简单的加减乘除保证正确性。