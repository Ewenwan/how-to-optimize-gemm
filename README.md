# How To Optimize Gemm wiki pages
https://github.com/flame/how-to-optimize-gemm/wiki


假设矩阵C = 矩阵A * 矩阵B； 矩阵A的shape为(M, K)，矩阵B的shape为(K, N)，矩阵C的shape为(M，N)。

Copyright by Prof. Robert van de Geijn (rvdg@cs.utexas.edu).

Adapted to Github Markdown Wiki by Jianyu Huang (jianyu@cs.utexas.edu).

# Table of contents

  * [The GotoBLAS/BLIS Approach to Optimizing Matrix-Matrix Multiplication - Step-by-Step](../../wiki#the-gotoblasblis-approach-to-optimizing-matrix-matrix-multiplication---step-by-step)
  * [NOTICE ON ACADEMIC HONESTY](../../wiki#notice-on-academic-honesty)
  * [References](../../wiki#references)
  * [Set Up](../../wiki#set-up)
  * [Step-by-step optimizations](../../wiki#step-by-step-optimizations)
  * [Computing four elements of C at a time](../../wiki#computing-four-elements-of-c-at-a-time)
    * [Hiding computation in a subroutine](../../wiki#hiding-computation-in-a-subroutine)
    * [Computing four elements at a time](../../wiki#computing-four-elements-at-a-time)
    * [Further optimizing](../../wiki#further-optimizing)
  * [Computing a 4 x 4 block of C at a time](../../wiki#computing-a-4-x-4-block-of-c-at-a-time)
    * [Repeating the same optimizations](../../wiki#repeating-the-same-optimizations)
    * [Further optimizing](../../wiki#further-optimizing-1)
    * [Blocking to maintain performance](../../wiki#blocking-to-maintain-performance)
    * [Packing into contiguous memory](../../wiki#packing-into-contiguous-memory)
  * [Acknowledgement](../../wiki#acknowledgement)

# Related Links
* [BLISlab: A Sandbox for Optimizing GEMM](https://github.com/flame/blislab)
* [GEMM: From Pure C to SSE Optimized Micro Kernels](http://apfel.mathematik.uni-ulm.de/~lehn/sghpc/gemm/)

# Acknowledgement
This material was partially sponsored by grants from the National Science Foundation (Awards ACI-1148125/1340293 and ACI-1550493).

_Any opinions, findings and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation (NSF)._
#  c/c++/python基本上是以行存储优先的，本文将以行存储优先作为基础进行优化分析。
    
考虑两种情况：
          
（1）当AB矩阵较小时，根据计算机结构可知，当从RAM中读取AB矩阵内存，根据局部性原理可以将AB矩阵放到cache中，因为cpu访问cache比访问主存的快。

（2）当AB矩阵较大时，超过cache大小时，根据矩阵乘的普通方法，由于访问“行优先存储的B矩阵”的时候内存不连续（读取B矩阵的一列），造成缓存cache频繁的换入换出，从RAM读取内存的次数大于AB矩阵的大小。

因此第一种优化方法：
       
1. 向量化（SIMD）

向量化可以使一条指令并行的使多个相同操作数执行相同的操作，减少每次循环迭代时评估操作类型的开销。
     
2. 内存对齐
         
内存对齐的原则：任何K字节的基本对象的地址必须都是K的倍数。
             
“假设 cache line 为 32B。待访问数据大小为 64B，地址在 0x80000001，则需要占用 3 条 cache 映射表项；若地址在 0x80000000 则只需要 2 条。内存对齐变相地提高了 cache 命中率。” 假定kernel一次计算执行4*4 大小的block, 根据MMult_4x4_7.c 和 MMult_4x4_8.c 代码，可以看出MMult_4x4_8.c使用了偏移量完成内存对齐。 
    
3. 分块计算，需要内存重排
   
   假设矩阵C = 矩阵A * 矩阵B; 矩阵A的shape为(M, K)，矩阵B的shape为(K, N)，矩阵C的shape为(M，N)。
   
   结果矩阵 C[ M, N ] = A[M, K]   *  B[ K, N]

最朴素版本:


1行x1列:
```c
for (int m = 0; m < M; m++) {    // 结果矩阵c的一行 
  for (int n = 0; n < N; n++) {  // 结果矩阵c的列行 
    C[m][n] = 0;
    for (int k = 0; k < K; k++) {
      C[m][n] += A[m][k] * B[k][n]; // A 矩阵的一行 * B 矩阵的一列
    }
  }
}

```

1行x4列： 循环展开
```c
for (int m = 0; m < M; m++) {    // 结果矩阵c的一行 
  for (int n = 0; n < N; n += 4) {  // 结果矩阵c的列行   一次计算四个结果
   
    // 结果储值位0
    C[m][n + 0] = 0;
    C[m][n + 1] = 0;
    C[m][n + 2] = 0;
    C[m][n + 3] = 0;
       
    for (int k = 0; k < K; k++) {
      C[m][n + 0] += A[m][k] * B[k][n + 0]; // A 矩阵的一行 * B 矩阵的四列
      C[m][n + 1] += A[m][k] * B[k][n + 1];
      C[m][n + 2] += A[m][k] * B[k][n + 2];
      C[m][n + 3] += A[m][k] * B[k][n + 3];
    }
  }
}

```
简单的观察即可发现，上述伪代码的最内侧计算使用的矩阵 A的元素是一致的。因此可以将 
A[m][k]读取到寄存器中，从而实现 4 次数据复用（这里不再给出示例）。一般将最内侧循环称作计算核（micro kernel）。进行这样的优化后，内存访问操作数量变为 
(2+1+1/4)MNK,其中 1/4是对 A优化的效果。

 

4行x4列： 循环展开
```c
for (int m = 0; m < M;   m += 4) {   // 结果矩阵c的4行    一次计算16个结果
  for (int n = 0; n < N; n += 4) {   // 结果矩阵c的4行
   
    // 结果储值位0
    C[m][n + 0] = 0;
    C[m][n + 1] = 0;
    C[m][n + 2] = 0;
    C[m][n + 3] = 0;
    // [m + 1] = 0
    // [m + 2] = 0
    // [m + 3] = 0
    
    for (int k = 0; k < K; k++) {
      C[m][n + 0] += A[m][k] * B[k][n + 0]; // A 矩阵的一行 * B 矩阵的四列
      C[m][n + 1] += A[m][k] * B[k][n + 1];
      C[m][n + 2] += A[m][k] * B[k][n + 2];
      C[m][n + 3] += A[m][k] * B[k][n + 3];
      // 
       C[m + 1][n + 0] += A[m + 1][k] * B[k][n + 0];
       //...
      // 
      C[m + 2][n + 0] += A[m + 2][k] * B[k][n + 0];
      //
      C[m + 3][n + 0] += A[m + 3][k] * B[k][n + 0];
    }
  }
}

```
同样地，将计算核心展开，可以得到下面的伪代码。这里我们将 1×4中展示过的 N 维度的计算简化表示。这种拆分可看成是 4×1×4，这样 
A 和 B的访存均可复用四次。由于乘数效应，4×4 的拆分可以将对输入数据的访存缩减到 2MNK+1/4
2MNK+1/4MNK+1/4MNK=5/2MNK。这相对于最开始的 4MNK 已经得到了 1.6X 的改进，这些改进都是通过 **展开循环后利用寄存器存储数据减少访存** 得到的。

