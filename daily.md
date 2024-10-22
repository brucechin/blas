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
6. java代码没有注释导致有些函数的实现无法判断是否正确，timeseries部分都是从pivot倒着计算的，很神奇。

### 7.19

1. tsargmax tsargmin是思路不一致的问题，考虑到业务代码上已经采用旧逻辑写完了，我不能修改，但按他的旧逻辑很难实现N^2的优化，ts count consecutive true函数放弃优化了
2. 今天需要用JNI来在java中调用C++，java坑还是挺多的，在CLASSPATH环境变量里需要设项目的root directory，不然java找不到要执行的对应的class，与此同时，java.library.path里包含了java运行时寻找动态链接库的路径，需要把gcc编译出来的so文件弄过去。java里loadlibrary的时候比如load “Hello”，那么lib文件名为“libHello.so”/"libHello.dll"
3. linker出bug的时候可以修改/etc/ld.so.conf然后/sbin/ldconfig，增加寻找的文件夹路径，gcc/g++ 编译的时候要带上-I include路径，-L 库文件路径
4. JNI里的JMatrix里再弄一个专门持有C++内存地址的对象句柄，当JMatrix析构的时候，这个对象调用回调函数释放C++对应地址上的内存，JMatrxi内部的各种计算也通过这个对象操作。
5. 还是要自己查一下如何用java管理调用C++的对象和内存，在C++这边写一个demo class，构造和析构函数输出一些东西，然后用javah生成的对应的C++的函数调用这个demo class对象，看看java这边对象析构时的回调函数能不能把C++这边对应的析构掉。

### 7.20

1. 发现java finalize不会让对象在main()函数结束的时候主动调用进而析构C++对象，要不要主动析构？？？
2. 针对问题1打算是C++这边Matrix在构造的时候通过JNI在JVM里申请内存然后用，这样java那边需要析构的时候jvm会判断内存够不够用需不需要gc
3. 发现JNI生成的native层C++实现的时候需要按java那边生成的函数顺序来实现，否则在link 动态链接库的时候会出现找不到函数定义的error
4. 目前需要考虑MatrixCalculator这个类里如果每个函数接口都要与原java接口大部分一样的话，那么我需要一个Matrix.java?然后问题2中C++通过JNI分配java Matrix类还是C++ matrix类实例化的内存？如果不再搞一个Matrix.java的话，那么MatrixCalculator.java里native的接口就需要每个接口传入Matrix对象的指针long型参数，而且Matrix的初始化也会有问题。

### 7.22

1. 发现JNI的接口NewDoubleArray分配的内存块如果太大了(比如大于30*sizeof(double))在用C++去修改内存上的值的时候会报错  	**libjvm.so+0x6f1dc8(这个地址不固定) JNIHandles::make**

   用env->SetDoubleArrayElements(target, 0, size, buffer)用buffer里的内容写入到target上时，target得是jdoublearray类型，buffer得是jdouble*类型

2. 但是用env->GetDoubleArrayElements(target, *isCopy)获取目前来看这样分配的内存也能在jvm内存不够用时被正确析构，但这会拷贝一下target，太耗内存了。

3. 直接在C++里实例化对象的话，java程序能一直运行不报错out of memory，但之前的对象也没有被析构很难受，估计是用虚拟内存在撑着。

### 7.23

1. 今天主要测了hftdev上单线程运行java和c++版本之间速度差异，最夸张的是matrixMul函数，C++利用openblas优化后快了一百倍，但没利用openblas写的函数大多性能相近，难有一个数量级左右的提升。

| （10k * 10k matrix） | java(s) | c++(s)                 |
| -------------------- | ------- | ---------------------- |
| smaller              | 0.8     | 0.19                   |
| between              | 0.16    | 0.39                   |
| condition            | 1.2     | 0.81                   |
| round                | 0.8     | 0.0003                 |
| tsstd                | 39.045  | 38.427                 |
| tsCountTrue          | 21      | 10.29 (optimized0.82)  |
| tsMax                | 48      | 13.57(optimized 0.057) |
| tsSkewness           | 25      | 25.07                  |
| tsMean               | 22      | 23.01                  |

结论：没用openblas的c++接口跟java性能差异不是很大，想要提升还需要算法上降低复杂度

### 7.24

1. 今天被安排了新的活，要写一个java的threadpool来管理机器上的CPU和内存资源，尽可能的提升任务并行度，提升资源利用率
2. 说要试试sun.misc.Unsafe类在jvm之外分配内存然后在上面搞事情，但这个跟在native层用C++实例化对象的效果差不多？
3. 配了一会clion+WSL失败了，IDEA很好用但sigar这个库一直报错，不知道还有啥获取系统cpu mem使用信息比较好用
4. JNI里的局部引用、全局引用、弱全局引用是native层代码引用JAVA里的实例，影响jvm对这些被native层代码引用的实例GC的时机的，jvm没法自动GC native层分配的内存
5. 使用java里的弱引用试图督促jvm释放过了生命周期的matrix，调用finalize后clear掉native里的空间，但失败了。目前唯一成功的是System.gc()强制回收
6. 在JVM运行空间中，对象的整个生命周期大致可以分为7个阶段：创建阶段（Creation）、应用阶段（Using）、不可视阶段（Invisible）、不可到达阶段（Unreachable）、可收集阶段（Collected）、终结阶段（Finalized）与释放阶段（Free）。

### 7.25

1. 今天学习使用了ThreadPoolExecutor和Timer类，做了一个定时收集CPU使用情况来更新ThreadPoolExecutor的pool size，感觉这种资源管理跟k8s很像啦，不过k8s集群还得调试几个月才能上线。

### 7.29

1. 打算看一下MXNet里java封装的时候如何处理C++那边的内存的
2. mxnet里java接口的调用demo程序中，ndarray在new出来之后没有手动析构，需要看一下怎么做到的。
3. 之前也有查过用java的weak reference或者phantom reference，但考虑到这样业务代码又要重写就很麻烦。

### 7.30

1. 今天看tensorflow的java wrapper发现他们讨论把tensor分配在java还是native(JNI/C++)这边的时候也折腾了很久，让人头秃

   > That sounds good! After looking further into JNI, I suggest the following:
   >
   > - Tensors are always backed by a direct ByteBuffer that is created on the Scala/Java side.
   > - They have an "asNative" method that (through JNI) calls "TF_NewTensor" in the C API. It also the JNI "NewGlobalRef" method to create a global reference to the direct ByteBuffer object. This guarantees that as long as TensorFlow native uses that byte buffer, it won't be garbage collected by the JVM. As a deallocator function we provide a function that the JNI "DeleteGlobalRef" method on the previously created reference.
   >
   > And that's the only interaction between the native TensorFlow tensors API and the Scala/Java API. The rest of it (i.e., elements access, slicing, etc.) is all handled on the Scala side. The same could apply to the current Java API and void the copying.
   
2. tf java里的tensor要么是自己调用close()，要么因为它是继承autocloseable类可以在try语句中new出来，语句块结束后会自动析构。

3. 又回到C++和java性能对比的问题上，突然意识到之前benchmark的错误在于分配10000乘以10000个double的连续内存是非常耗时的，可能比计算还要耗时，那么两者对比的时间差异不大可能是由此导致的。我又修改成128 128的矩阵连续运算5000次发现差距依然不大，那么说明小块内存分配相比小矩阵计算也很耗时。现在需要测在矩阵不同size下整个计算过程中为result矩阵分配空间消耗了多少时间。

4. 然而经测试，内存分配并不怎么消耗时间：）总用时占比小于1%.

5. 进一步直接测C++调用openblas和裸写的速度差异。

6. 查漏补缺，把c++和java版本大部分的接口计算出来的结果二进制保存到文件上再读出来比较结果差异是否超过error_bound，主要是dirty work。

### 7.31

1. 昨天java和c++ matrix计算结果的正确性测试代码已经跑通了，目前有以下函数结果不正确：这几个都是因为本来想调用intel mkl的库但是没法同时link openblas和mkl所以把mkl的实现注释掉了，内部函数为空导致的，打算先裸写一版下来。（✔为修改完毕的）
   - [x] mul(0,0)
   - [x] max(0,0)
   - [x] min(0,0)
   - [x] rank(0,0)
   - [x] tsskewness, tscorr, tscov
   - [x] abs, log, minus, exp, inverse, sqrt, floor, round
2. 发现有些bug是因为访存炫技玩脱了导致的，不该作死的。
3. 另外一个问题是C++里的NAN二进制写到文件里后java再读出来转换成0.0了

### 8.1

1. C++和java库作比较时，如果相对误差控制在10^-9以内也都满足条件
2. rank的也解决了。。居然是因为冒泡排序没写对。。zs了

### 8.2

1. 开始利用AutoCloseable类写Matrix.java + JNI，目前readMatrix/writeMatrix还有构造和析构函数等都可用了，迈出了坚实的第一步（说明在C++这一层的内存访问正常且正常释放）
2. 进一步实现了Matrix的getRow/Col Vector，set/getElement等API，目前Matrix.java类已经基本实现完成，下一步要实现MatrixCalculator类

### 8.5

1. 周末把MatrixCalculator.java和native层都写掉且跑通了
2. 今天早上把MatrixCorrectenessTest.java跑通且都对了，除了tsCorr还有点问题会core dump。因为之前已经验证c++和原java的库结果是一致的，所以这次是JNI的和c++的做比较，由此可见提升代码复用性可以减少很多活。
3. 刚才除了一个小插曲，在1000*1000的规模下，matrixMul和mul的用时相差不到一半，我以为出了问题， 后来发现提升到十倍规模后差距才出来，interesting。
4. **之前写的dynamicthreadpool需要增加一个功能，强行停止现在所有的任务，最好再加一个graceful stop。**
5. 需要再实现一个demo，向MatrixCalculator传入空的hashmap，然后按照一定规律填满key value对
6. 提升IO性能
7. 测试INF NAN处理的如何。

速度测试结果如下，时间单位为毫秒，每个任务均执行50次

| API\SIZE     | 1k*1k      | 3k*3k        |
| ------------ | ---------- | ------------ |
| add\sub      | 152/236    | 3228/2088    |
| div\mul      | 65/147     | 2140/1873    |
| matrixMul    | 251/22894  | 6096/inf     |
| max/min      | 65/587     | 2161/4611    |
| \><=         | 64/124     | 815/979      |
| between      | 79/114     | 2123/971     |
| and          | 230/333    | 2151/2971    |
| or           | 234/294    | 2152/2856    |
| not          | 41/31      | 403/273      |
| condition    | 221/373    | 4407/4148    |
| rank         | 45290/7643 | 125735/inf   |
| round        | 146/193    | 3370/1764    |
| floor        | 98/192     | 2939/1940    |
| abs          | 38/108     | 2327/1374    |
| minus        | 38/113     | 2321/1168    |
| sqrt         | 323/295    | 5012/2759    |
| log          | 691/1744   | 8437/15162   |
| exp          | 719/9263   | 9239/72449   |
| sign         | 39/146     | 2375/2045    |
| inverse      | 38/115     | 2404/999     |
| signedpow    | 347/162    | 5237/1344    |
| shift/delay  | 53/172     | 2597/1201    |
| delta/ratio  | 53/148     | 2602/1292    |
| sum/product  | 2416/3779  | 24945/29868  |
| tsmax/tsmin  | 2643/5800  | 27532/54590  |
| tsargmax/min | 4000/4613  | 39214/42580  |
| tsrank       | 687/1862   | 66914/155220 |
| tsmean       | 516/5665   | 51012/49497  |
| tsstd        | 587/6019   | 57414/54488  |
| tsskewness   | 624/6719   | 61384/61910  |
| tskurtosis   | 664/6737   | 65123/66868  |
| tscov        | 7905/9838  | 72612/92589  |
| tscounttrue  | 2421/3033  | 23658/29409  |
| decayLinear  | 3632/3851  | 34831        |
| neuralize    | 654/1204   | 11283        |

发现由于jvm的懒惰释放内存政策，且以C++为后端的BLAS库能在try代码块结束后立刻释放不用的内存，相同程序下，利用以C++为后端BLAS在用内存量往往远小于java裸写版本.

**此外，应该针对ts函数继续优化复杂度，性能还有可提升的空间。**



### 8.6

1. 对随机生成的矩阵中含有NAN或者INF时进行正确性测试
2. 支持自定义增加新接口？如java段传入一个hashmap，是在c++端计算出结果后传一个array回java再在java这边填满，还是在c++端填满后传回hashmap回java？？？
3. 大矩阵下C++加速效果不明显。

### 8.7

1. 之前的关于native层写一个java hashmap然后传回java端已经实现，其中的坑主要是java hashmap的key和value必须是object，不能直接用jint和jdouble jboolean之类的。
2. 之前的线程池shutdown之后会接着往下执行代码，增加了新功能是，必须要等queue里的tasks执行完之后才shutdown成功，接着往下执行其他代码。
3. 接了一个新需求，需要把现实中的时间转换成trading time，这个是完全连续的，因为现实时间中不是一直在做交易，转换成连续的方便他们策略组跑模型吧
4. 晚上把几个optimized timeseries的函数fix了，之后可以把更复杂的几个ts函数优化一下复杂度

### 8.8

1. 实现trading time和现实时间转换的C++库函数，要编译成.a文件
2. 继续优化几个ts的函数，能提升一到两个数量级的性能很舒服
3. 今天在测正确性的时候发现了一个bug是，如果最后结果都是很小的数，比如e-10级别的，那么他们的相对误差会很大，这个时候还是用绝对误差靠谱一些。不过一般只会在matrix超级大的情况下才会有这种极端结果算出来。最后定位用绝对误差不超过0.001，偶尔在mul, matrixMul, product这种会出e10数量级以上的大数时才会超过这个error_bound，然而按相对误差来算是低于e-12的。
4. 今天查出来的一些低级错误包括for循环里k--写成k++，sumxx写漏一个x，还有就是copy一些想重复使用的代码块的时候，在新的函数里面有些变量其实换了名字，不能直接copy过来
5. 下一步计划把每隔num个数执行一次的ts函数也优化一下，这块难度不大

### 8.9

1. 发现在函数参数带num的API实现中有重大错误是在写res的时候越界了。主要是res和mat1 mat2的size是不一样的，没调用setElement裸写内存访问的时候写错了。。
2. 接了一个新需求要写一个matrix pool，减少在复杂大矩阵预算过程中，内存反复allocate、release的overhead（因为MatrixCalculator里面的API都是先为res分配好内存再往上填）

### 8.12

1. 找到了合适的C++ CSV文件处理库，准备开始写TradingCalendar了
2. 把之前写的matrix pool添加了多线程Lock()，并撰写了多线程的测试代码。

### 8.15

1. 前两天calendar写好了，在debug遇到以下问题：

   1. # [error: passing xxx as 'this' argument of xxx discards qualifiers](https://stackoverflow.com/questions/5973427/error-passing-xxx-as-this-argument-of-xxx-discards-qualifiers)
   
2. 今天完全debug完成了，明天主要写一些单元测试，测一下各个函数的情况

### 8.16

1. 今天跑单元测试，发现生成大的随机数也是一件需要考量的事情
2. 无法跟py代码作比较测正确性，估计要手算以下了。。

### 8.21

1. 摸了两天🐟之后终于又安排了新的活，主要是XML文件的解析，要有三个功能，一个是分析多个xml文件的公共部分，一个是多个xml文件随机挑一个base文件，比较base与其他文件的diff，最后一个是将delta文件作为base文件update的参照，输出新的output.xml
2. 忙了半个下午发现这活还挺棘手的，因为xml文件的tag tree中每个tag都可能有attribute，不方便把tag tree转换成一个唯一识别符

### 8.22

1. 把xml文件里每一个line沿途的path上的tag和attribute转变成他unique key，然后用这个单层的词典进行diff和求相同line的操作，但是merge的时候发现难以把这个四不像dict转回xml了
2. 最后解决办法是用ElementTree把xml文件转成一个tree结构，然后对两个tree进行merge然后再转成xml格式写回到文件

### 8.26

1. 优化了项目结构，增加一些代码注释


### 8.38

1. 项目从本机WSL上移植到hftdev上还是有各种小问题，自己在开发机上啥权限都没有非常烦躁，装个pip包不行，把so文件或者.h文件放在lib或者include文件夹下也不行，难受的一批.
2. 发现旧MatrixCalculator版本有个问题是当matrix size过大时，里面的int2Double等转换函数就出界了.
3. 发现一个问题：C++ matrix speedtest时相比Java(backend c++)版本要快的(用时少10-60%不等)，这是native层带来的额外开销导致的，但是在开发机上测试native层带来的overhead的时候发现又没有那么多了，很是奇怪。
4. 在开发机上同时跑两种实现的JAVA MatrixCalculator的时候，C++ backend的还是总体快一些，尤其是在优化过的TS系列函数上，想优化最好的办法还是优化算法啊。

在1000*1000size的矩阵测试中，C++后端的Java native实现与原来的普通java实现耗时对比(每个函数均测试100次求总用时)

![](C:\Users\10261\Desktop\blas\perf-compare-1k-1k-100times.PNG)

在2000*2000size的矩阵测试中，测试结果如下：

![](C:\Users\10261\Desktop\blas\perf-compare-2k-2k-100times.PNG)